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

import scala.annotation.tailrec
import scala.util.parsing.combinator.RegexParsers

/** Parser for MD tags from SAM files. Convert MD tags to substitution variants.
  *
  * @param seq The dna sequence from a read to parse and convert to a list of variants.
  */
class MDParser(seq: String) extends RegexParsers {
  /* There are 3 possible variants:
   *
   * 1. Substitutions - consumes both the read and the reference
   * 2. Insertions - consumes only the read
   * 3. Deletions - consumes only the reference
   */

  // These are state variables: reference and read indices/pointers.
  private var referenceLoc: Int = 0
  private var readLoc: Int = 0

  /************************* AUXILIARY FUNCTIONS ******************************/

  /** Find the position of the nth = in a sequence.
    *
    * @param sequence The sequence to search.
    * @param n Number of ='s to match.
    * @return Position of the nth =
    */
  def nthMatch(sequence: String, n: Int): Int = {

    /* Internal tail recursive function for matching.
     *
     * @param sequence - The sequence in the read.
     * @param n - number of ='s to match
     * @param i - read index
     * @return - position of the nth "="
     */
    @tailrec def nthMatchRec(sequence: String, n: Int, i: Int): Int = {
      if (n <= 0) i // Base case, we're done
      else if(i >= sequence.length)
        // Sanity check. We have a badly formed MD tag.
        throw new Exception(s"Parse error in: $sequence, number of matches ($n) exceeds sequence length (${sequence.length}).")

      else if (sequence.charAt(i) == '=') nthMatchRec(sequence, n - 1, i + 1)  // Match an =
      else nthMatchRec(sequence, n, i + 1)                              // Did not match an =
    }

    nthMatchRec(sequence, n, 0)
  }

  /*************************** PRODUCTIONS  ***********************************/

  /* A match could be a match or insertion. Matches advance both read and reference
   * indices. Insertions only advance the read index. The reference is always advanced by
   * number of matches. The read is advanced until we get number of "="'s equal to number
   * of matches.
   */
  def matches: Parser[String] =
  """\d+""".r ^^ {
    m => {
      // We always increase reference by number of matches
      val numMatches = m.toInt
      referenceLoc = referenceLoc + numMatches

      // Read loc is trickier because of insertions
      readLoc = readLoc + nthMatch(seq.substring(readLoc), numMatches)

      "" // No variant in a match, since we don't care about insertions.
    }
  }

  def md: Parser[List[String]] = matches ~ rep(subOrDel ~ matches) ^^ {
    case matchesHead ~ matchesTail => matchesTail.foldLeft(List(matchesHead)) {
      case (m, v ~ ms) => m ++ List(v) ++ List(ms)
    }.filter(_ != "")
  }

  def subOrDel: Parser[String] = substitution | deletion

  def substitution: Parser[String] = """[acgtnACGTN]""".r ^^ {
    // A substitution can only be a single base according to the SAM format guidelines.
    base => {
      val variant: String =
        if(base == "N") ""
        else (referenceLoc+1).toString + seq.substring(readLoc, readLoc+1)

      // Advance the indices for reference and read
      referenceLoc = referenceLoc + 1
      readLoc = readLoc + 1

      variant
    }
  }

  def deletion: Parser[String] = """\^[acgtnACGTN]+""".r ^^ {
    deletion => {
      referenceLoc = referenceLoc + deletion.length   // Deletions only affect reference index
      ""    // No variant, we don't care about deletions.
    }
  }
}
