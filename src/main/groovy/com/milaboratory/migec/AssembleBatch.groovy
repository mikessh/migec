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

def DEFAULT_MODE = "0:1"
def cli = new CliBuilder(usage: 'AssembleBatch [options] checkout.filelist.txt histogram.estimates.txt')
cli.p(args: 1, 'number of threads to use')
cli.m(longOpt: "assembly-mode", args: 1, argName: 'assembly mode in format X:X',
        "Identifier(s) of read(s) to assemble. Default: \"$DEFAULT_MODE\".")
cli.c("compressed output")

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

boolean compressed = opt.c
String mode = opt.m
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
String filelistFileName = opt.arguments()[0], estimatesFileName = opt.arguments()[1]

def sampleFileNamesMap = filelistFileName.readLines()[1..-1].collectEntries {
    def splitLine = it.split('\t')
    [(splitLine[0..1].join('\t')): splitLine.size() == 4 ? splitLine[2..3] : [splitLine[2], '-']]
}

estimatesFileName.readLines()[1..-1].each {
    def splitLine = it.split('\t')
    def sampleName = splitLine[0..1].join('\t')
    def sampleFileNames = sampleFileNamesMap[sampleName]

    ESTIMATES_HEADER = "SAMPLE_ID\tFILE_TYPE\t" +
            "TOTAL_READS\tTOTAL_MIGS\t" +
            "OVERSEQ_THRESHOLD\tCOLLISION_THRESHOLD\t" +
            "UMI_QUAL_THRESHOLD\tUMI_LEN"

    def totalUmis = splitLine[3].toInteger(), overseqThreshold = splitLine[4].toInteger(),
        collThreshold = splitLine[5].toInteger(), umiQualThreshold = Byte.parseByte(splitLine[6]),
        umiLen = splitLine[7].toInteger(), filterCollisions = false

    if (collThreshold > overseqThreshold && totalUmis < Math.pow(4, umiLen))
        filterCollisions = true    // safe to filter collisions

    def args = [sampleFileNames].flatten().join(" ")

}

//int minCount

//class SampleData {
//    final String path
//    final int[] overSeq, overSeqHeader, collisions, collisionsHeader
//}
//boolean collisionFiltering

// Collision filtering if umi_size ? #UMIs
