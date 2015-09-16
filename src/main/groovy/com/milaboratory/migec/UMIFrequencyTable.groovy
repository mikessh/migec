/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.milaboratory.migec

//========================
//          CLI
//========================
def cli = new CliBuilder(usage: 'UmiFrequencyTable [options] input.fastq[.gz] output')
def scriptName = getClass().canonicalName
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}
def inputFileName = opt.arguments()[0],
    outputFileName = opt.arguments()[1]

if (new File(outputFileName).parentFile)
    new File(outputFileName).parentFile.mkdirs()

//========================
//      MISC UTILS
//========================


def getUmiEntry = { String header ->
    def splitHeader = header.split(" ")
    def umiEntry = splitHeader.find { it.startsWith("UMI:") }
    if (umiEntry == null) {
        println "[${new Date()} $scriptName] Error: no UMI header in input. Terminating"
        System.exit(2)
    }
    umiEntry.split(":")
}

//========================
//         BODY
//========================
println "[${new Date()} $scriptName] Reading data from $inputFileName, collecting UMI headers"
def header
def reader = Util.getReader(inputFileName)
def umiMap = new HashMap<String, Integer>()
def n = 0
while ((header = reader.readLine()) != null) {
    if (!header.startsWith("@")) {
        println "[${new Date()} $scriptName] Not a FASTQ!"
        System.exit(2)
    }

    def umi = getUmiEntry(header)[1]
    umiMap.put(umi, (umiMap.get(umi) ?: 0) + 1)

    reader.readLine()
    reader.readLine()
    reader.readLine()

    if (++n % 500000 == 0)
        println "[${new Date()} $scriptName] $n reads processed"
}
println "[${new Date()} $scriptName] Finished processing $n reads"

new File(outputFileName).withPrintWriter { pw ->
    umiMap.entrySet().collect().sort { -it.value }.each { pw.println(it.key + "\t" + it.value) }
}

