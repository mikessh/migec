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

import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicIntegerArray

def cli = new CliBuilder(usage: 'Histogram [options] checkout.filelist.txt output_prefix')
cli.q(args: 1, argName: 'read quality (phred)', 'barcode region quality threshold. Default: 15')
cli.p(args: 1, 'number of threads to use')
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
byte umiQualThreshold = opt.q ? Byte.parseByte(opt.q) : Util.DEFAULT_UMI_QUAL_THRESHOLD

def scriptName = getClass().canonicalName
int nBins = 17

def scale = { Integer value ->
    (int) (Math.min(Math.log((double) value) / Math.log(2.0), nBins - 1))
}

def samples = new File(opt.arguments()[0]).readLines().collect {
    it.split("\t")
}

def outputFilePrefix = opt.arguments()[1]

if (new File(outputFilePrefix).parentFile)
    new File(outputFilePrefix).parentFile.mkdirs()

def HEADER = "SAMPLE_ID\tFILE_TYPE\t" + (0..<nBins).collect { (int) Math.pow(2, it) }.join("\t")

new File("${outputFilePrefix}.overseq.txt").withPrintWriter { oWriter ->
    oWriter.println(HEADER)

    new File("${outputFilePrefix}.overseq-units.txt").withPrintWriter { ouWriter ->
        ouWriter.println(HEADER)

        new File("${outputFilePrefix}.collision1.txt").withPrintWriter { cWriter ->
            cWriter.println(HEADER)

            new File("${outputFilePrefix}.collision1-units.txt").withPrintWriter { cuWriter ->
                cuWriter.println(HEADER)

                samples.each { sampleEntry ->
                    println "[${new Date()} $scriptName] Processing ${sampleEntry[0]} (${sampleEntry[1]})"

                    // Accumulate UMIs
                    def umiCountMap = new HashMap<String, Integer>()
                    def reader = Util.getReader(sampleEntry[2])

                    String header

                    int nReads = 0
                    while ((header = reader.readLine()) != null) {
                        def umi = Util.getUmi(header, umiQualThreshold)

                        if (umi != null)
                            umiCountMap.put(umi, (umiCountMap.get(umi) ?: 0) + 1)

                        reader.readLine()
                        reader.readLine()
                        reader.readLine()

                        if (++nReads % 1000000 == 0)
                            println "[${new Date()} $scriptName] Processed $nReads, ${umiCountMap.size()} UMIs so far"
                    }

                    println "[${new Date()} $scriptName] Processed $nReads, ${umiCountMap.size()} UMIs total"

                    def overseqHist = new AtomicIntegerArray(nBins), overseqHistUnits = new AtomicIntegerArray(nBins),
                        collisionHist = new AtomicIntegerArray(nBins), collisionHistUnits = new AtomicIntegerArray(nBins)

                    def nUmis = new AtomicInteger()

                    GParsPool.withPool THREADS, {
                        umiCountMap.eachParallel { Map.Entry<String, Integer> umiEntry ->
                            int thisCount = umiEntry.value, bin = scale(thisCount)

                            // Append to cumulative overseq
                            overseqHist.addAndGet(bin, thisCount)
                            overseqHistUnits.incrementAndGet(bin)

                            // Calculate 1-mm collisions
                            char[] umi = umiEntry.key.toCharArray()

                            for (int i = 0; i < umi.length; i++) {
                                for (int j = 0; j < 4; j++) {
                                    char prevChar, nt = Util.NTS[j]
                                    if (umi[i] != nt) {
                                        prevChar = umi[i]
                                        umi[i] = nt
                                        def otherCount = umiCountMap.get(new String(umi))
                                        if (otherCount != null && thisCount < otherCount) {
                                            collisionHist.addAndGet(bin, thisCount)
                                            collisionHistUnits.incrementAndGet(bin)
                                        }
                                        umi[i] = prevChar
                                    }
                                }
                            }

                            nUmisCurrent = nUmis.incrementAndGet()

                            if (nUmisCurrent % 10000 == 0)
                                println "[${new Date()} $scriptName] Collecting stats, $nUmisCurrent UMIs processed"
                        }
                    }

                    def row = sampleEntry[0..1].join("\t")

                    oWriter.println(row + "\t" + Util.toString(overseqHist))
                    cWriter.println(row + "\t" + Util.toString(collisionHist))
                    ouWriter.println(row + "\t" + Util.toString(overseqHistUnits))
                    cuWriter.println(row + "\t" + Util.toString(collisionHistUnits))
                }
            }
        }
    }
}