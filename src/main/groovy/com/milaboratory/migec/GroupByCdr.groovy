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

package com.milaboratory.migec

import java.util.zip.GZIPInputStream


def W_S = "5", W_M_Q = "15"
def cli = new CliBuilder(usage: 'GroupByCdr [options] input.fastq.gz ' +
        '(with CDR3 field in header, generatedby CdrBlast with --cdr3-fastq-file set) ' +
        '[cdr3_filter_table from FilterCdrBlastResults or -] output')
cli.f('Quality filter reads with sliding window')
cli._(longOpt: 'window-min-qual', args: 1, "Min average quality in window. Default: $W_M_Q")
cli._(longOpt: 'window-size', args: 1, "Size of sliding window. Default: $W_S")
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(-1)
}

boolean qualFilter = opt.f
int windowSize = Integer.parseInt(opt.'window-size' ?: W_S),
    windowMinQual = Integer.parseInt(opt.'window-min-qual' ?: W_M_Q) * windowSize

def isGoodQual = { String qual ->
    def qualArray = qual.toCharArray().collect { (int) it - 33 }
    int currentQual = 0
    for (int i = 0; i < windowSize; i++)
        currentQual += qualArray[i]

    if (currentQual < windowMinQual)
        return false

    for (int i = windowSize; i < qual.size(); i++) {
        currentQual -= qualArray[i - windowSize]
        currentQual += qualArray[i]

        if (currentQual < windowMinQual)
            return false
    }
    return true
}

def scriptName = getClass().canonicalName
def inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[2]
def cdr3Set = new HashSet<String>()
boolean onlyListedCdrs = opt.arguments()[1] != "-"
if (onlyListedCdrs) {
    new File(opt.arguments()[1]).splitEachLine("\t") { line ->
        if (line.size() > 1)
            cdr3Set.add(line[0])
        else
            cdr3Set.add(line[2])
    }
}

def getReader = { String fname ->
    new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fname)) :
            new FileInputStream(fname)))
}

class CdrReadData {
    def seqs = new ArrayList<String>()
    def cdrFroms = new ArrayList<Integer>()
    int minX = Integer.MAX_VALUE, minY = Integer.MAX_VALUE
}

def cdrDataMap = new HashMap<String, CdrReadData>()
def reader = getReader(inputFileName)
def header
println "[${new Date()} $scriptName] Reading sequences from $inputFileName"
int n = 0
while ((header = reader.readLine()) != null) {
    def splitHeader = header.split(" ")
    def cdrEntry = splitHeader.find { it.startsWith("CDR3:") }
    if (cdrEntry == null) {
        println "[${new Date()} $scriptName] Error: no CDR3 entry in header: $header. Terminating"
        System.exit(-1)
    }
    cdrEntry = cdrEntry.split(":")
    def cdrSeq = cdrEntry[1]
    int cdrFrom = Integer.parseInt(cdrEntry[2])

    if (!onlyListedCdrs || cdr3Set.contains(cdrSeq)) {
        def fullSeq = reader.readLine()
        reader.readLine()
        def qual = reader.readLine()
        if (!qualFilter || isGoodQual(qual)) {
            def cdrReadData = cdrDataMap.get(cdrSeq)
            if (cdrReadData == null)
                cdrDataMap.put(cdrSeq, cdrReadData = new CdrReadData())

            cdrReadData.seqs.add(fullSeq)
            cdrReadData.cdrFroms.add(cdrFrom)
            cdrReadData.minX = Math.min(cdrReadData.minX, cdrFrom)
            cdrReadData.minY = Math.min(cdrReadData.minY, fullSeq.size() - cdrFrom)
        }
    } else {
        reader.readLine()
        reader.readLine()
        reader.readLine()
    }
    n++
}
println "[${new Date()} $scriptName] Finished, $n reads with ${cdrDataMap.size()} unique CDR3 regions."

new File(outputFileName).withPrintWriter { pw ->
    cdrDataMap.each {
        def cdrReadData = it.value
        for (int i = 0; i < it.value.seqs.size(); i++) {
            pw.println(cdrReadData.seqs.size() + "\t" + it.key + "\t" + cdrReadData.seqs[i].substring(cdrReadData.cdrFroms[i] - cdrReadData.minX,
                    cdrReadData.cdrFroms[i] + cdrReadData.minY))
        }
    }
}