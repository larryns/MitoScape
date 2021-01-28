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

import org.bdgenomics.formats.avro.Alignment

/* Optional SAM tags for reference.
 *
 * XM: Cigar string
 * RG: Read group
 * MD: String for mismatching positions. Used to try to get the SNP/indel calls without looking at the reference.
 * NH: * Number of reported alignments that contains the query in the current alignment.
 * HI: Query hit index, indicating the alignment record is the i-th one stored in the current record.
 * NM: Edit distance to the reference, including ambiguous bases but excluding clipping.
 * SM: Template-independent mapping quality.
 * XQ: Non-normalized mapping quality score.
 * X2: * Second best XQ score among all multimapping alignments for a given read. X2 = 0 if there is only one alignment.
 * XO: Output type:
 *      NM - nomapping, the entire read could not be aligned.
 *      CU - * concordant unique, both ends of a paired-end read could be aligned concordantly aligned to a single
 *           position in the genome. The order of optimal alignments are concordant, paired and then unpaired. CU is
 *           the ideal alignment.
 *      CM - * concordant multiple.
 *      CX - * concordant multiple excess.
 *      CT - concordant translocation.
 *      CC - * concordant circular.
 *      PI - paired unique inversion.
 *      PS - paired unique scramble.
 *      PL - paired unique long.
 *      PC - paired unique circular.
 *      PM - paired multiple.
 *      PX - paired multiple excess.
 *      HU, HM, HT, HC - halfmapping unique, halfmapping multiple, halfmapping translocation, halfmapping circular
 *      HX - halfmapping multiple excess.
 *      UU, UM, UT, UC - unpaired unique, unpaired multiple, unpaired translocation, unpaired circular.
 * XG: Method that generated the alignment:
 *      A: suffix array method.
 *      B: GMAP alignment produced from suffix array.
 *      M: hash table alignment
 *      T: terminal alignment
 *      O: merging of overlaps
 */

/** Converts a single sam record to a feature to be used by the machine learning classifier.
  * IMPORTANT: This class is designed to operate on a single read, and return per read, NOT
  * paired read.
  *
  * @param read The read to convert to a feature.
  */
abstract class Feature(read: Alignment) {
  /** Parse the attributes of the read. The attributes we are interested in are:
    *
    * XQ: ignored
    * XO: ignored
    * XM: ignored
    * SM: ignored
    * NM: Edit distance to the reference
    * HI: ignored
    * NH: Number of alignments
    * XG: ignored
    * RG: ignored
    * X2: second best mapping score
    *
    * @return A map: attribute name --> value of attribute.
    */
  private def parseAttributes(): Map[String, String] = {
    // takes a token and returns a tuple with the code and value
    def parseToken(token: String): (String, String) = {
      val tokens = token.split(":")

      // We must have at least 3 tokens
      if(tokens.size < 3) (tokens(0), "")
      else (tokens(0), tokens(2))
    }

    read
      .getAttributes        // Get the attributes
      .split("\\s+")  // Then tokenize the attributes
      .map(parseToken).toMap  // Parse the tokens and convert to a map
  }

  protected val attribs: Map[String, String] = parseAttributes()

  /** Converts the MD tag and read sequence to variants for a single read.
    *
    */
  protected def getVariants: List[String] = {
    val md: String = read.getMismatchingPositions
    val parser = new MDParser(read.getSequence)

    parser.parseAll(parser.md, md) match {
      case parser.Success(result: List[String], _) => result
      case _ => throw new Exception("Parse error in MD tag: " + md)
    }
  }

  /** Check if the read is a valid/legitimate paired read.
    *
    * @return  True if read is legitimate.
    */
  def isValid: Boolean = read.getPrimaryAlignment &&
    read.getReadPaired &&
    read.getProperPair &&
    read.getMateMapped &&
    !read.getSupplementaryAlignment

  // The following are all the attributes/features we are going to use in the machine learning algorithm.
  val editDistance: Int = attribs.getOrElse("NM", "0").toInt
  val numAlignments: Int = attribs.getOrElse("NH", "1").toInt
  val mappingScore: Int = attribs.getOrElse("XQ", "0").toInt
  val readName: String = read.getReadName
  val mapQ: Int = if(isValid) read.getMappingQuality else 0
  val chromosome: String = read.getReferenceName

  // For printing
  override def toString: String = readName + " " + attribs.mkString("; ")
}

/** A nuclear chromosome feature.
  *
  * @param read The read to convert to a feature.
  */
class NucFeature(read: Alignment) extends Feature(read) {

}

/** A MT chromosome feature. This class stores variants.
  *
  * @param read The read to convert to a feature.
  */
class MTFeature(read: Alignment) extends Feature(read) {
  val variants: List[String] = getVariants
}
