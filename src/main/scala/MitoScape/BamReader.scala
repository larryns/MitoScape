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

package MitoScape

/** Provides methods for classifying NUMTs from bam files.
  *
  * LOG NOTES: I dropped the MappingScore (Alignment score), since we're taking the best
  * alignment, this score doesn't change much.
  */

import org.apache.spark.broadcast.Broadcast
import org.apache.spark.rdd.RDD
import org.apache.spark.sql._
import org.apache.spark.sql.types._
import org.apache.spark.sql.functions.{collect_list, sum}
import org.bdgenomics.adam.rdd.ADAMContext

/** Classifies a single bam file. Each bam file is treated independently.
 */

/** Classifies a single bam file.
  *
  * @param spark the spark session we are working in.
  */
abstract class BamReader(spark: SparkSession) extends Serializable {
  // I don't want to use the implicit so I'll make my own adam context.
  val ac: ADAMContext = new ADAMContext(spark.sparkContext)

  /** Summarizes/aggregates the dataframe (df) bby removing redundant rows.
    *
    * @param df - input dataframe with redundat row.
    * @return Dataframe where each row corresponds to a unique read.
    */
  protected def summarize(df: DataFrame): DataFrame

  /** Converts the features to a dataset.
    *
    * @return Dataframe of features.
    */
  def DF: DataFrame
}

/** Read a nuclear bam file and extract the features for the classifier.
  *
  * @param prefix prefix of bam file to process (including path)
  * @param spark the spark session we are working in.
  */
class NucReader(prefix: String, spark: SparkSession, NUMTs: Map[String, List[(Int, Int, Float)]]
) extends BamReader(spark) {
    /** Summarizes the dataframe, by aggregating the rows per read.
    *
    * @param df - input dataframe with redundant row.
    * @return Dataframe where each row corresponds to a unique read.
    */
  protected def summarize(df: DataFrame): DataFrame = df
      .groupBy("Read")
      .agg(
        sum("NTMapQ").alias("NTMapQ"),
        sum("NTNumAlignments").alias("NTNumAlignments"),
        sum("NTEditDist").alias("NTEditDist"),
        sum("NTScore").alias("NTScore")
      )

  /** Gets the Nuclear Reads in a dataframe format to be used by the classifier.
    *
    * @return Dataframe of features.
    */
  def DF: DataFrame = {
    def NUMTOverlaps(feature: NucFeature): Int =
      NUMTs
        .getOrElse(feature.chromosome, Nil)
        .map(_._3)
        .sum
        .toInt

    val featureSchema = StructType(
      Seq(
        StructField(name = "Read", dataType = StringType, nullable = false),
        StructField(name = "NTMapQ", dataType = IntegerType, nullable = false),
        StructField(name = "NTNumAlignments", dataType = IntegerType, nullable = false),
        StructField(name = "NTEditDist", dataType = IntegerType, nullable = false),
        StructField(name = "NTScore", dataType = IntegerType, nullable = false)
      )
    )

    val featureRDD: RDD[Row] = ac
        .loadAlignments(prefix + "_NT.bam")
      .rdd.map(new NucFeature(_))
      .filter(_.isValid)
      .map {
        feature =>
          Row(
            feature.readName,
            feature.mapQ,
            feature.mappingScore,
            feature.numAlignments,
            feature.editDistance,
            NUMTOverlaps(feature)
          )
    }

    // At this point, there are multiple rows with the same reads, so we need to aggregate them.
    summarize(spark.createDataFrame(featureRDD, featureSchema))
  }

}


/** Reads a MT bam file
  *
  * @param prefix prefix of bam file to process (including path)
  * @param spark the spark session we are working in.
  * @param ldScoresBV - LD file for retrieving the linkage scores for mito reads; must be a broadcast variable
  */
class MTReader(prefix: String, spark: SparkSession, ldScoresBV: Broadcast[LD]) extends BamReader(spark) {
  // For the toDF function call below, we need an encoder
  import spark.implicits._

  /** Removes the repeated rows for every row by aggregating them.
    *
    * @param df - input dataframe to summarize/aggregate
    * @return Summarized dataframe, i.e. a dataframe where each row is a unique read name.
    */
  protected def summarize(df: DataFrame): DataFrame = df
      .groupBy(
        "Read"
      ) // Group by read name to get the pairs
      .agg(
      sum("MTMapQ"),
      sum("MTNumAlignments"),
      sum("MTEditDist"),
      collect_list("variants").alias("VariantList")
    ).map(row => (
        row.getString(0), // read name
        row.getLong(1), // mapping quality
        row.getLong(2), // number of alignments
        row.getLong(3), // edit distance
        computeR(row.getSeq[List[String]](4).toList) // Get Linkage Disequilibrium scores
      ))
      .toDF("Read", "MTMapQ", "MTNumAlignments", "MTEditDist", "LD")

  /** Get the features in dataframe format.
    *
    * @return Dataframe of features.
    */
  def DF: DataFrame = {
    // Get a list of variants from a read. First group the reads to get the pairs together.
    val featureSchema = StructType(
      Seq(
        StructField(name = "Read", dataType = StringType, nullable = false),
        StructField(name = "MTMapQ", dataType = IntegerType, nullable = false),
        StructField(name = "MTNumAlignments", dataType = IntegerType, nullable = false),
        StructField(name = "MTEditDist", dataType = IntegerType, nullable = false),
        StructField(name = "Variants", dataType = ArrayType(StringType), nullable = true)
      )
    )

    val featureRDD: RDD[Row] = ac
      .loadAlignments(prefix + "_MT_MD.bam")
      .rdd.map(new MTFeature(_))
      .filter(_.isValid)
      .map {
        feature =>
          Row(
            feature.readName,
            feature.mapQ,
            feature.numAlignments,
            feature.editDistance,
            feature.variants
          )
      }

    /* At this point, we have a data frame with the components of the read that we need
     * minus the LD scores. We need to take the sequences and now convert them to LD scores.
     */
    summarize(spark.createDataFrame(featureRDD, featureSchema))
  }

  /** Compute R/LD score from a list of variants stored as a list of string.
    *
    * @param variants - Variants on the same paired-read
    * @return - Aggregate LD score for all the variants on the paired read.
    */
  private def computeR(variants: List[List[String]]): Int = variants
      .flatten
      .combinations(2)
      .map(varPair => ldScoresBV.value.getLD(varPair.head, varPair(1)))
      .sum

}
