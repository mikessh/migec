/**
 Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.migec

import java.util.zip.GZIPInputStream

//========================
//          CLI
//========================
def cli = new CliBuilder(usage: 'UMIFrequencyTable [options] input.fastq[.gz] output')
def scriptName = getClass().canonicalName
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(0)
}
def inputFileName = opt.arguments()[0],
    outputFileName = opt.arguments()[1]

if (new File(outputFileName).parentFile)
    new File(outputFileName).parentFile.mkdirs()

//========================
//      MISC UTILS
//========================
def getReader = { String fname ->
    new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fname)) :
            new FileInputStream(fname)))
}

def getUmiEntry = { String header ->
    def splitHeader = header.split(" ")
    def umiEntry = splitHeader.find { it.startsWith("UMI:") }
    if (umiEntry == null) {
        println "[${new Date()} $scriptName] Error: no UMI header in input. Terminating"
        System.exit(-1)
    }
    umiEntry.split(":")
}

//========================
//         BODY
//========================
println "[${new Date()} $scriptName] Reading data from $inputFileName, collecting UMI headers"
def header
def reader = getReader(inputFileName)
def umiMap = new HashMap<String, Integer>()
def n = 0
while ((header = reader.readLine()) != null) {
    if (!header.startsWith("@")) {
        println "[${new Date()} $scriptName] Not a FASTQ!"
        System.exit(-1)
    }

    def umi = getUmiEntry(header)[1]
    umiMap.put(umi, (umiMap.get(umi) ?: 0) + 1)

    reader.readLine()
    reader.readLine()
    reader.readLine()

    if (++n % 500000 == 0)
        println "[${new Date()} $scriptName] $n reads processed"
}

new File(outputFileName).withPrintWriter { pw ->
    umiMap.entrySet().collect().sort { -it.value }.each { pw.println(it.key + "\t" + it.value) }
}