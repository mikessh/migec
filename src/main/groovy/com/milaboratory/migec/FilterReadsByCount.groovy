package com.milaboratory.migec

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

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

def DEFAULT_THRESHOLD = "10"

def cli = new CliBuilder(usage:
        'FilterReadsByCount [options] input.fastq[.gz] output.fastq[.gz]\n[for benchmarking purposes]')
cli.t(args: 1, "Count threshold. $DEFAULT_THRESHOLD")

def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def threshold = (opt.t ?: DEFAULT_THRESHOLD).toInteger()
def inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

def reader = Util.getReader(inputFileName)

class ReadData {
    int count = 0
    final long[] qualArr

    ReadData(String seq) {
        this.qualArr = new long[seq.length()]
    }

    void append(String qual) {
        count++
        for (int i = 0; i < qual.length(); i++)
            qualArr[i] += Util.qualFromSymbol(qual.charAt(i))
    }

    String finalizeQual() {
        char[] qual = new char[qualArr.length]

        for (int i = 0; i < qualArr.length; i++)
            qual[i] = Util.symbolFromQual((byte) (qualArr[i] / count))

        new String(qual)
    }
}

def readDataMap = new HashMap<String, ReadData>()

while (reader.readLine() != null) {
    def seq = reader.readLine()
    reader.readLine()
    def qual = reader.readLine()

    def readData = readDataMap[seq]

    if (readData == null)
        readDataMap.put(seq, readData = new ReadData(seq))

    readData.append(qual)
}

def writer = Util.getWriter(outputFileName)

readDataMap.each {
    if (it.value.count >= threshold)
        writer.writeLine("@GroupedRead Count:$it.value.count\n" + it.key + "\n+\n" + it.value.finalizeQual())
}

writer.close()