package MitoScape

/*
 * Copyright 2021 Larry N. Singh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import org.apache.spark.broadcast.Broadcast
import org.apache.spark.sql.{DataFrame, SparkSession}
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.adam.rdd.read.AlignmentDataset
import org.bdgenomics.adam.sql.Alignment
import org.seqdoop.hadoop_bam.SAMFormat

import scala.annotation.tailrec

/** This module classifies reads based on heuristics and a machine learning classifier.
  */

object MTClassify {
  /** main function
    *
    * @param args - command line arguments
    */
  def main(args: Array[String]): Unit = {
    val usage =
      """
        Usage: mtclassify
                          --prefix <prefix of bam file>
                          --out <output file>
                          --prob <probability threshold>
                          --threads <number of threads>
      """.stripMargin

    type OptionMap = Map[String, String]

    /* Recursively parse the command line options.
     *
     * @param optMap  - resultant map containing parameters and corresponding values
     * @param argxList - the command line arguments to parse.
     * @return optMap
     */
    @tailrec def optionMap(optMap: OptionMap, argList: List[String]): OptionMap = {
      // Recursive matching of options.
      argList match {
        // Check if we have all the required options and we're done.
        case Nil =>
          if (optMap.isEmpty) {
            println(usage)
            sys.exit(-1)
          } else optMap
        case "--prefix" :: value :: rest =>
          optionMap(optMap ++ Map("prefix" -> value), rest)
        case "--out" :: value :: rest =>
          optionMap(optMap ++ Map("out" -> value), rest)
        case "--prob" :: value :: rest =>
          optionMap(optMap ++ Map("maxProb" -> value), rest)
        case "--threads" :: value :: rest =>
          optionMap(optMap ++ Map("threads" -> value), rest)
        case "--ld" :: value :: rest =>
          optionMap(optMap ++ Map("ld" -> value), rest)
        case "--numt" :: value :: rest =>
          optionMap(optMap ++ Map("numt" -> value), rest)
        case "--classifier" :: value :: rest =>
          optionMap(optMap ++ Map("class" -> value), rest)
        case option :: _ => println("Unknown option: " + option)
          sys.exit(1)
      }
    }
    val options: OptionMap = optionMap(Map(), args.toList)

    /* Helper function to get options.
     *
     * @param opt - option we are looking for
     * @return String value of the option
     */
    def getOpt(opt: String): String = {
      options.get(opt) match {
        case Some (value) => value
        case None =>
          println ("ERROR: " + opt + " must be specified.")
          println(usage)
          sys.exit(-1)
      }
    }

    // Print error and exit
    def printAndExit(msg: String) = {
      Console.err.println(msg)
      sys.exit(-1)
    }

    // Utility variables and objects
    // My understanding is that since the resources below are unmmanaged and and not stored relative to the class's
    // package, I need to specify an absolute path with "/".
    val ldFile: String = options.getOrElse("ld", "mitomap.ld")
    val numtFile: String = options.getOrElse("numt", "NUMTs_hg38.txt")
    val classifierPath: String = options.getOrElse("class", "MTClassifierModel.RF")

    val prefix = getOpt("prefix")
    val outFile = getOpt("out")

    val maxProb = try {
        options.getOrElse("maxProb", "0.5").toDouble
    } catch {
      case _: Exception => printAndExit("ERROR: Invalid probability specified.")
    }
    if(maxProb < 0.0 || maxProb > 1.0) {
      printAndExit("ERROR: Probability (" + maxProb + ") option must be in [0.0, 1.0].")
    }

    val threads = try {
      options.getOrElse("threads", "1").toInt
    } catch {
      case _: Exception => printAndExit("ERROR: Invalid number of threads specified.")
    }

    // Create a spark session

    val spark: SparkSession = SparkSession
      .builder
      .master("local[" + threads + "]")
      .appName("MTClassify")
      .config("spark.sql.autoBroadcastJoinThreshold ", -1)
      .config("spark.driver.bindAddress", "127.0.0.1")
      .getOrCreate

    // Load the NUMTs from file.
    val NUMTs = MTClassifierModel.loadNUMTs(numtFile)

    /* Create the MitoMap dataframe. It's important that the executors can't have
     * spark sessions, so you have to create everything that uses a sparksession here in the
     * main object. We need the mitomap lookup in the executors though, so we need to broadcast
     * it.
     */
    val ldScores: LD = new LD(ldFile, spark)
    val ldScoresBV: Broadcast[LD] = spark.sparkContext.broadcast(ldScores)

    // Load the dataframe with the bans
    val bamDF: DataFrame = MTClassifierModel.DF(
      prefix,
      spark,
      NUMTs,
      ldScoresBV,
      0.5   // Initial label
    )

    def filterAlignments(mtAlignments: AlignmentDataset): AlignmentDataset =
      // If there are no reads in bamDF, then no need to do any filtering.
      if(!bamDF.isEmpty && maxProb > 0.0) {
        // Load the RF model
        val model = MTClassifierModel.loadRF(classifierPath)

        // Create an R model and perform predictions
        val readsDF: DataFrame = model
          .transform(MTClassifierModel prepareDF bamDF)

        // Filter the reads that have low probability
        import spark.sqlContext.implicits._

        val probFilterReads = MTClassifierModel
          .getMaxProb(readsDF, spark).where($"MaxProb" >=  maxProb)

        // Get the reads
        val MTReads = readsDF
          .where("prediction = " + MTClassifierModel.MT_LABEL)
          .select("Read")

        // load the alignments again from file.
        mtAlignments
          .transformDataset(
            ds => {
              /* Filter:
               *
               * 1. Join with the mito reads.
               * 2. Filter low probability
               */
              ds.join(MTReads, MTReads.col("Read") === ds.col("readName"), "leftsemi")
                .join(probFilterReads, probFilterReads.col("Read") === ds.col("readName"), "leftsemi")
                .as[Alignment]
            }
          )
      } else {
        Console.err.println("No filtering done.")
        mtAlignments
      }

    // Load the alignments again from file
    val mtAlignments: AlignmentDataset = new ADAMContext(spark.sparkContext)
      .loadAlignments(prefix + "_MT.bam") // AlignmentRecordRDD

    filterAlignments(mtAlignments)
      .saveAsSam(outFile, asType=Some(SAMFormat.BAM), asSingleFile=true)

    // Done, wrap-up
    ldScoresBV.unpersist(true)
    ldScoresBV.destroy()

    spark.close()

    // We require a return exit status.
    System.exit(0)
  }
}
