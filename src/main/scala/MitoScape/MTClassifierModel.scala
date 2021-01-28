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

/** This object creates the machine learning (Random Forest) model for classifying the ambiguous
  * mtDNA reads. Use this object to create and save the model.
  */

import org.apache.spark.broadcast.Broadcast
import org.apache.spark.ml.Pipeline
import org.apache.spark.ml.classification.{GBTClassifier, RandomForestClassificationModel, RandomForestClassifier}
import org.apache.spark.ml.evaluation.MulticlassClassificationEvaluator
import org.apache.spark.ml.feature.RFormula
import org.apache.spark.ml.linalg.Vector
import org.apache.spark.sql.{DataFrame, Row, SparkSession}
import org.apache.spark.sql.functions.{avg, lit, stddev_samp}
import org.apache.spark.sql.types.DoubleType

import scala.io.Source

object MTClassifierModel {
  val MT_LABEL = 0.0
  val RHO0_LABEL = 1.0
  val RF_NUM_TREES = 128   // Number of trees to use for the RandomForest Classifier

  /* I removed MT/NTNumAlignments because these variables don't seem to have much effect.
   * IMPORTANT: I originally had MTMapQ and NTMapQ as features, but it was obvious that the sequencing
   * performed for the MT reads had a different distribution from the mapq scores from the WGS (rho0) cells.
   * Whether the difference in distribution is due to the difference in sequencing platforms and preparation or
   * there is a legitimate difference in distributions is unclear. I cannot resolve this issue, so I
   * removed mapq scores entirely. In hindsight, the amplified mtDNA reads and rho0 reads should have
   * been barcoded, multiplexed and sequenced together. I hope this approach can be taken in the future.
   *
   * An alternative is to perform a Z-score transformation on all the Mapping qualities so that they can be compared.
   * At this time, the algorithm is working very well, so I don't think this is necessary, but it may be worth
   * considering in the future.
   */
  val R_FORMULA = "label ~ MTEditDist + LD  + NTEditDist + NTScore + MTNumAlignments + NTNumAlignments"

  def time[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block    // call-by-name
    val t1 = System.nanoTime()
    println("Elapsed time: " + (t1 - t0) + "ns")
    result
  }

  /** Load the mitochondrial and nuclear DNA reads and join them.
    *
    * @param prefix of the MT and Nuclear DNA sequencing data
    * @param spark session
    * @param NUMTs map storing the NUMT information.
    * @param ldScoresBV linkage disequilibrium scores for the mitochondrial data
    * @param label label either case (1) or control (0)
    * @return Dataframe with the merged MT and Nuclear scores for RF classificaiton
    */
  def DF(prefix: String, spark: SparkSession, NUMTs: Map[String, List[(Int, Int, Float)]], ldScoresBV: Broadcast[LD], label: Double): DataFrame = {
    // Get all the reads from the bam aligned to just the nuclear chromosomes.
    val nucBamFeaturesDF: DataFrame = new NucReader(prefix, spark, NUMTs).DF

    // Get all the reads from the bam aligned to the mito chromosomes.
    val mtBamFeaturesDF: DataFrame = new MTReader(prefix, spark, ldScoresBV).DF
    if(mtBamFeaturesDF.isEmpty) {
      Console.err.println(s"WARNING: MT File: $prefix is empty!")
      sys.exit(0)
    }

    // Left-join the nuclear and MT features.
    val features: DataFrame = mtBamFeaturesDF.join(nucBamFeaturesDF, "Read")

    // add the label and normalize the mapping qualities
    NormalizeMapQ(features.withColumn("label", lit(label)))
  }

  /** Convert a dataframe into a 3 column vector with ReadName, MaxProb, and Prediction, where MaxProb
    * is the maximum probability assigned to the label.
    *
    * @param df - input dataframe
    * @return dataframe with ReadName, MaxProb and Prediction
    */
  def getMaxProb(df: DataFrame, spark: SparkSession): DataFrame = {
    import spark.implicits._

    // Get the indices of the probabilities we need.
    val probColIndex = df.columns.indexOf("probability")
    val predictionColIndex = df.columns.indexOf("prediction")

    // Create a dataframe with just the read, maximum probability of the classifier and the prediction
    df.map {
      row => {
        val probMax: Double = row(probColIndex).asInstanceOf[Vector].toArray.max
        (row.getString(0), probMax, row.getDouble(predictionColIndex))
      }
    }.toDF("Read", "MaxProb", "Prediction")
  }

  /** The mapping qualities for the mito sequencing and rho-zero sequencing are
    * different. If you want to compare the mapping quality scores, then you must
    * normalize both the MT and NT mapping qualities together for each of the
    * mito sequencing and rho-zero sequencing, independently.
    *
    * @param df - unnormalized dataframe
    * @return DataFrame with normalized MT MappingQ and NT MappingQ scores
    */
  def NormalizeMapQ(df: DataFrame): DataFrame = {
    // It's possible that the datframe is empty.
    if(df.isEmpty)
      df
        .withColumn("MTMapQ", lit(null).cast(DoubleType))
        .withColumn("NTMapQ", lit(null).cast(DoubleType))
    else {
      // The MT and NT mapping qualities are in two separate columns. We need the mean
      // and stddev of both columns to do a Z-transform.
      val summaryStats: Row = df
        .select(df.col("MTMapQ").alias("MapQ"))
        .union(df.select(df.col("NTMapQ").alias("MapQ")))
        .select(
          avg("MapQ").alias("Avg"),
          stddev_samp("MapQ").alias("SD")
        ).first()

      val mean = summaryStats.getDouble(0)
      val sd = summaryStats.getDouble(1)

      // Now perform the Z transformation
      df
        .withColumn("MTMapQ", (df.col("MTMapQ") - mean) / sd)
        .withColumn("NTMapQ", (df.col("NTMapQ") - mean) / sd)
    }
  }

  /** Convert a line delimited by tabs into a paired tuple with the first element being the first element.
    *
    * @param line input line
    * @return A paired tuple with the first element being the chromosome, and the second element the rest.
    */
  def toPair(line: String): Array[String] = line.trim().split("\t")

  /** Stores the NUMTs read from file. The format of the map is:
    * chromosome -> (start, end, score)
    *
    */
  def loadNUMTs(numtFile: String): Map[String, List[(Int, Int, Float)]] = {
    val numtReader = Source.fromFile(numtFile)

    val numtMap = numtReader
      .getLines()
      .map(toPair)
      .toList
      .groupBy(_.head)
      .map {
        case (chrom, positions) =>
          (chrom, positions.map(pos => (pos(1).toInt, pos(2).toInt, pos(3).toFloat)))
      }

    numtReader.close()
    numtMap
  }

  /** Prepare/transform the data to be used by the classifier.
    *
    * @param df - dataframe to be transformed
    * @return transformed dataframe
    */
  def prepareDF(df: DataFrame): DataFrame = {
    val formula = new RFormula().setFormula(R_FORMULA)
    val fittedRF = formula.fit(df)

    fittedRF.transform(df)
  }

  // Create the RF classifier, tweak the parameters and try them.
  def testRFModel(trainingProb: Float, trainingDF: DataFrame, spark: SparkSession): Double = {
    // First convert the data frame to a form to be used for the classifiers.

    // Create a random split of the data into training and test data.
    val preparedDF = prepareDF(trainingDF)
    val splits = preparedDF.randomSplit(Array(trainingProb, 1.0-trainingProb))
    val (trainingData, testData) = (splits(0), splits(1))

    // Train the model
    val rf = new RandomForestClassifier()
      .setNumTrees(RF_NUM_TREES)

    // Chain indexers and forest in a Pipeline.
    val pipeline = new Pipeline()
      .setStages(Array(rf))

    // Train model
    val model = pipeline.fit(trainingData)

    // Make predictions.  Change max prob to tune the accuracy.
    val predictionsOrig = model.transform(testData)
    val predictionsMaxProb = getMaxProb(predictionsOrig, spark).where("MaxProb >= 0.0")
    val predictions = predictionsOrig
      .join(
        predictionsMaxProb,
        predictionsOrig.col("Read") === predictionsMaxProb.col("Read"),
        "left_semi"
      )

    // Select (prediction, true label) and compute test error.
    val evaluator = new MulticlassClassificationEvaluator()
      .setMetricName("accuracy")

    // Return the accuracy
    evaluator.evaluate(predictions)
  }

  // Create the Gradient Boosting classifier, tweak the parameters and try them.
  def testGBModel(trainingProb: Float, trainingDF: DataFrame): Double = {
    // First convert the data frame to a form to be used for the classifiers.
    val formula = new RFormula()
      .setFormula(R_FORMULA)
    val fittedRF = formula.fit(trainingDF)
    val preparedDF = fittedRF.transform(trainingDF)

    // Create a random split of the data into training and test data.
    val splits = preparedDF.randomSplit(Array(trainingProb, 1.0-trainingProb))
    val (trainingData, testData) = (splits(0), splits(1))

    // Train a GradientBoostedTrees model.
    val gbtClassifier = new GBTClassifier()
      .setMaxDepth(8)
      .setMaxIter(50)

    // Chain indexers and forest in a Pipeline.
    val pipeline = new Pipeline()
      .setStages(Array(gbtClassifier))

    // Train model
    val model = pipeline.fit(trainingData)

    // Make predictions.
    val predictions = model.transform(testData)

    // Select (prediction, true label) and compute test error.
    val evaluator = new MulticlassClassificationEvaluator()
      .setMetricName("precision")

    // Return the accuracy
    evaluator.evaluate(predictions)
  }

  /** Train the Random Forest (RF) classifier
    *
    * @param trainingDF - Training Data to build the Random Forest (RF) classifier
    * @param RFPath - path to store the built RF model
    * @return the trained model
    */
  def trainRF(trainingDF: DataFrame, RFPath: String): RandomForestClassificationModel = {
    def trainRF: RandomForestClassificationModel = {
      val rf = new RandomForestClassifier().setNumTrees(RF_NUM_TREES)
      val formula = new RFormula().setFormula(R_FORMULA)
      val fittedRF = formula.fit(trainingDF)
      val preparedDF = fittedRF.transform(trainingDF)

      rf.fit(preparedDF)
    }

    val model = trainRF
    model
        .write
        .overwrite()
        .save(RFPath)
    model
  }

  /** Load the Random forest classifier model from disk.
    *
    * @param RFPath - path where the RF is stored
    * @return RandomForestClassificationModel trained on MT and Rho0 data.
    */
  def loadRF(RFPath: String): RandomForestClassificationModel = RandomForestClassificationModel.load(RFPath)
}
