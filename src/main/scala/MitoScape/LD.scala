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

import org.apache.spark.sql.SparkSession

/** Stores linkage disequilibrium (LD) scores obtained from mitomap for the mitochondrial
  * variants.
  */
class LD(ldFile: String, spark: SparkSession) extends Serializable {
  private val RMultiplier: Int = 100000

  /* ldFile contains the raw data information. We can't use the csv method because a delimiter
   * has to be a single character not a regex. We get an unparsed array of lines containing
   * the linkage information.
   */
  import spark.implicits._

  /** Creates a map of the linkage disequilibrium (LD) scores.
    *
    * @return Map of (Variant1, Variant2) --> R LD score.
    */
  type LDMap = Map[(String, String), Int]

  Console.err.println(s"Loading LD file: $ldFile")
  private val ldMap: LDMap = spark.read
    .option("header", "true")
    .option("delimiter", "\t")
    .option("inferSchema", "true")
    .csv(ldFile)
    .toDF("Variant1", "Variant2", "R")
    .map {
      row => ((row.getString(0), row.getString(1)), (row.getDouble(2) * RMultiplier).toInt)
    }
    .filter(_._2 != 0)
    .collect()
    .toMap

  /** Get the linkage score for two variants.
    *
    * @param var1 - variant 1
    * @param var2 - variant 2
    * @return LD score corresponding to both variants
    */
  def getLD(var1: String, var2: String): Int = {
    ldMap.get((var1, var2)) match {
      case Some(ld: Int) => ld
      case _ => ldMap.get((var2, var1)) match {
        case Some(ld: Int) => ld
        case _ => 0
      }
    }
  }
}
