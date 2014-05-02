/**
 Copyright 2013 Mikhail Shugay (mikhail.shugay@gmail.com)

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

package migec.control

@Grab(group = 'org.codehaus.gpars', module = 'gpars', version = '1.0.0')

import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger
import java.util.zip.GZIPInputStream

def cli = new CliBuilder(usage: 'groovy Histogram [options] checkout.filelist.txt output_prefix')
cli.q(args: 1, argName: 'read quality (phred)', 'barcode region quality threshold. Default: 15')
cli.p(args: 1, 'number of threads to use')
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    cli.usage()
    return
}
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
int umiQualThreshold = opt.q ?: 15

def scriptName = getClass().canonicalName
int nBins = 17
def nts = ['A', 'T', 'G', 'C']
def getReader = { String fname ->
    new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fname)) :
            new FileInputStream(fname)))
}
def scale = { Integer value ->
    (int) (Math.min(Math.log((double) value) / Math.log(2.0), nBins - 1))
}
def qualFromSymbol = { char symbol ->
    (int) symbol - 33
}
def getUmi = { String header ->
    def splitHeader = header.split(" ")

    def umiEntry = splitHeader.find { it.startsWith("UMI:") }

    String umi = umiEntry.split(":")[1]
    for (int i = umiEntry.length() - umi.length(); i < umiEntry.length(); i++)
        if (qualFromSymbol(umiEntry.charAt(i)) < umiQualThreshold)   // quality can contain :
            return null
    umi
}

def fileNameList = new File(opt.arguments()[0]).readLines().collect { it.split("\t")[0] },
    outputFilePrefix = opt.arguments()[1]

new File(outputFilePrefix).mkdirs()

// Load UMIs  if (!(new File(logFileName).exists())) {
boolean oHeader = !new File("${outputFilePrefix}.overseq.txt").exists(),
        ouHeader = !new File("${outputFilePrefix}.overseq-units.txt").exists(),
        cHeader = !new File("${outputFilePrefix}.collision1.txt").exists(),
        cuHeader = !new File("${outputFilePrefix}.collision1.txt").exists()

new File("${outputFilePrefix}.overseq.txt").withWriterAppend { oWriter ->
    if (oHeader)
        oWriter.println("SAMPLE_ID\t" + (0..<nBins).collect { (int) Math.pow(2, it) }.join("\t"))
    new File("${outputFilePrefix}.overseq-units.txt").withWriterAppend { ouWriter ->
        if (ouHeader)
            ouWriter.println("SAMPLE_ID\t" + (0..<nBins).collect { (int) Math.pow(2, it) }.join("\t"))

        new File("${outputFilePrefix}.collision1.txt").withWriterAppend { cWriter ->
            if (cHeader)
                cWriter.println("SAMPLE_ID\t" + (0..<nBins).collect { (int) Math.pow(2, it) }.join("\t"))

            new File("${outputFilePrefix}.collision1-units.txt").withWriterAppend { cuWriter ->
                if (cuHeader)
                    cuWriter.println("SAMPLE_ID\t" + (0..<nBins).collect { (int) Math.pow(2, it) }.join("\t"))
                fileNameList.each { fileName ->
                    println "[${new Date()} $scriptName] Processing $fileName"

                    // Accumulate UMIs
                    def umiCountMap = new HashMap<String, Integer>()
                    def reader = getReader(fileName)
                    String header
                    int nReads = 0
                    while ((header = reader.readLine()) != null) {
                        def umi = getUmi(header)
                        if (umi != null)
                            umiCountMap.put(umi, (umiCountMap.get(umi) ?: 0) + 1)
                        reader.readLine()
                        reader.readLine()
                        reader.readLine()
                        if (++nReads % 1000000 == 0)
                            println "[${new Date()} $scriptName] Processed $nReads, ${umiCountMap.size()} UMIs so far"
                    }
                    println "[${new Date()} $scriptName] Processed $nReads, ${umiCountMap.size()} UMIs total"

                    def overseqHist = new AtomicInteger[nBins], overseqHistUnits = new AtomicInteger[nBins],
                        collisionHist = new AtomicInteger[nBins], collisionHistUnits = new AtomicInteger[nBins]
                    for (int i = 0; i < nBins; i++) {
                        overseqHist[i] = new AtomicInteger()
                        collisionHist[i] = new AtomicInteger()
                        overseqHistUnits[i] = new AtomicInteger()
                        collisionHistUnits[i] = new AtomicInteger()
                    }
                    def nUmis = new AtomicInteger()
                    GParsPool.withPool THREADS, {
                        umiCountMap.eachParallel { Map.Entry<String, Integer> umiEntry ->
                            int thisCount = umiEntry.value

                            // Append to cumulative overseq
                            overseqHist[scale(thisCount)].addAndGet(thisCount)
                            overseqHistUnits[scale(thisCount)].incrementAndGet()

                            // Calculate 1-mm collisions
                            char[] umi = umiEntry.key.toCharArray()
                            for (int i = 0; i < umi.length; i++) {
                                for (int j = 0; j < 4; j++) {
                                    char prevChar, nt = nts[j]
                                    if (umi[i] != nt) {
                                        prevChar = umi[i]
                                        umi[i] = nt
                                        def otherCount = umiCountMap.get(new String(umi))
                                        if (otherCount != null && thisCount < otherCount) {
                                            collisionHist[scale(thisCount)].addAndGet(thisCount)
                                            collisionHistUnits[scale(thisCount)].incrementAndGet()
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

                    oWriter.println(fileName + "\t" + overseqHist.collect().join("\t"))
                    cWriter.println(fileName + "\t" + collisionHist.collect().join("\t"))
                    ouWriter.println(fileName + "\t" + overseqHistUnits.collect().join("\t"))
                    cuWriter.println(fileName + "\t" + collisionHistUnits.collect().join("\t"))
                }
            }
        }
    }
}