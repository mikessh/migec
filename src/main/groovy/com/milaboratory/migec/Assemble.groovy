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

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicInteger

//========================
//          CLI
//========================
def DEFAULT_MODE = "1:1", DEFAULT_MIN_COUNT = "10", DEFAULT_PARENT_CHILD_RATIO = "0.1"
def cli = new CliBuilder(usage: 'Assemble [options] R1.fastq[.gz] [R2.fastq[.gz] or -] output_prefix ' +
        '[assembly_log, to append]')
cli.q(args: 1, argName: 'read quality (phred)',
        "barcode region quality threshold. Default: $Util.DEFAULT_UMI_QUAL_THRESHOLD")
cli._(longOpt: 'assembly-mode', args: 1, argName: 'X:Y, X=0/1, Y=0/1',
        "Mask for read(s) in pair that should be assembled. " +
                "0:0 indicates single read from overlapped paired-end. Default: \"$DEFAULT_MODE\".")
cli.p(args: 1,
        "number of threads to use. Default: all available processors")
cli.c("compressed output")

cli._(longOpt: 'alignment-file-prefix', args: 1, argName: 'string',
        "File name prefix to output multiple alignments generated during assembly, for \"com.milaboratory.migec.control.BacktrackSequence\"")
cli.m(longOpt: 'min-count', args: 1, argName: 'integer',
        "Min number of reads in MIG. Should be set according to 'Histogram.groovy' output. Default: $DEFAULT_MIN_COUNT")
cli._(longOpt: 'filter-collisions',
        "Collision filtering. Should be set if collisions (1-mismatch erroneous UMI sequence variants) " +
                "are observed in 'Histogram.groovy' output")
cli._(longOpt: 'collision-ratio', args: 1, argName: 'double, < 1.0',
        "Min parent-to-child MIG size ratio for collision filtering. Default value: $DEFAULT_PARENT_CHILD_RATIO")
cli._(longOpt: 'assembly-offset', args: 1, argName: 'integer',
        'Assembly offset range. Default: 3')
cli._(longOpt: 'assembly-mismatches', args: 1, argName: 'integer',
        'Assembly max mismatches. Default: 5')
cli._(longOpt: 'assembly-anchor', args: 1, argName: 'integer',
        'Assembly anchor region half size. Default: 10')

def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 3) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(-1)
}

//========================
//         PARAMS
//========================
def scriptName = getClass().canonicalName

// Parameters
boolean compressed = opt.c, filterCollisions = opt.'filter-collisions'
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
byte umiQualThreshold = opt.q ? Byte.parseByte(opt.q) : Util.DEFAULT_UMI_QUAL_THRESHOLD
int minCount = Integer.parseInt(opt."min-count" ?: DEFAULT_MIN_COUNT)
double collisionRatioThreshold = Double.parseDouble(opt.'collision-ratio' ?: DEFAULT_PARENT_CHILD_RATIO)

// I/O
def fastq1 = opt.arguments()[0],
    fastq2 = opt.arguments()[1],
    outputFilePrefix = opt.arguments()[2]

// Assembly parameters
def offsetRange = Integer.parseInt(opt.'assembly-offset' ?: '5'),
    maxMMs = Integer.parseInt(opt.'assembly-mismatches' ?: '5'),
    anchorRegion = Integer.parseInt(opt.'assembly-anchor' ?: '10')

// I/O parameters
boolean paired = fastq2 != "-"
def assemblyIndices = paired ? (opt.'assembly-mode' ?: DEFAULT_MODE).split(":").collect { Integer.parseInt(it) > 0 } :
        [true, false]
boolean bothReads = assemblyIndices[0] && assemblyIndices[1], overlapped = false
if (!assemblyIndices.any()) {
    assemblyIndices[0] = true
    overlapped = true
}

// Misc output
String alignmentFilePrefix = opt.'alignment-file-prefix'
String logFileName = opt.arguments().size() > 3 ? opt.arguments()[3] : null

if (new File(outputFilePrefix).parentFile)
    new File(outputFilePrefix).parentFile.mkdirs()
if (opt.'alignment-file-prefix' && new File(alignmentFilePrefix).parentFile)
    new File(alignmentFilePrefix).parentFile.mkdirs()
if (logFileName && new File(logFileName).parentFile)
    new File(logFileName).parentFile.mkdirs()

//=================================
//   PRE-LOAD DATA FOR ASSEMBLY
//=================================
def migData = new HashMap<String, Map<String, Integer>>[2]
migData[0] = new HashMap<String, Map<String, Integer>>(1000000)
if (assemblyIndices[0] && assemblyIndices[1])
    migData[1] = new HashMap<String, Map<String, Integer>>(1000000)
def putData = { int readId, String umi, String seq ->
    def seqCountMap = migData[readId].get(umi)
    if (seqCountMap == null)
        migData[readId].put(umi, seqCountMap = new HashMap<String, Integer>())
    seqCountMap.put(seq, (seqCountMap.get(seq) ?: 0) + 1)
}

int nReads = 0, nGoodReads = 0
println "[${new Date()} $scriptName] Pre-loading data for $fastq1, $fastq2.."
def reader1 = Util.getReader(fastq1), reader2 = paired ? Util.getReader(fastq2) : null
String header1, seq1
String seq2 = ""
int MIN_READ_SZ = 2 * anchorRegion + 1 + offsetRange
while ((header1 = reader1.readLine()) != null) {
    seq1 = reader1.readLine()
    reader1.readLine()
    reader1.readLine()

    if (paired) {
        reader2.readLine() // Skip header of read 2
        seq2 = reader2.readLine()
        reader2.readLine()
        reader2.readLine()
    }

    if (seq1.length() > MIN_READ_SZ && (!paired || seq2.length() > MIN_READ_SZ)) {
        def umi = Util.getUmi(header1, umiQualThreshold)

        if (umi != null && seq1.length() > 0) {
            if (assemblyIndices[0])
                putData(0, umi, seq1)
            if (bothReads)
                putData(1, umi, seq2)
            else if (assemblyIndices[1])
                putData(0, umi, seq2) // put all to 1st file if mode=0,1
            nGoodReads++
        }
        if (++nReads % 500000 == 0)
            println "[${new Date()} $scriptName] Processed $nReads reads, " +
                    "unique UMIs so far ${migData[0].size()}"
    }
}
println "[${new Date()} $scriptName] Processed $nReads reads, " +
        "unique UMIs ${migData[0].size()}"

//=================================
//   PERFORM ASSEMBLY
//=================================
def writeQueue = new LinkedBlockingQueue<String[]>(2048)
def suffix = assemblyIndices[0] ? (overlapped ? "R12" : "R1") : "R2"
def outputFileName1 = "${outputFilePrefix}_${suffix}.fastq",
    outputFileName2 = bothReads ? "${outputFilePrefix}_R2.fastq" : "-"
def writer1 = Util.getWriter(outputFileName1, compressed),
    writer2 = bothReads ? Util.getWriter(outputFileName2, compressed) : null
def detailsWriter1 = alignmentFilePrefix ? Util.getWriter("${alignmentFilePrefix}_${suffix}.asm", compressed) : null,
    detailsWriter2 = bothReads && alignmentFilePrefix ? Util.getWriter("${alignmentFilePrefix}_${suffix}.asm", compressed) : null

println "[${new Date()} $scriptName] Starting assembly.."
def writeThread = new Thread({  // Writing thread, listening to queue
    String[] result
    while (true) {
        result = writeQueue.take()

        if (result.length == 0)
            break

        writer1.writeLine(result[0])

        if (bothReads)
            writer2.writeLine(result[1])

        if (alignmentFilePrefix)
            detailsWriter1.writeLine(result[2])

        if (alignmentFilePrefix && bothReads)
            detailsWriter2.writeLine(result[3])
    }

    writer1.close()

    if (bothReads)
        writer2.close()

    if (alignmentFilePrefix)
        detailsWriter1.close()

    if (alignmentFilePrefix && bothReads)
        detailsWriter2.close()
} as Runnable)
writeThread.start()

def nMigs = new AtomicInteger(),
    nReadsInMigs = new AtomicInteger(), nCollisions = new AtomicInteger()
def nGoodMigs = new AtomicInteger[3], nReadsInGoodMigs = new AtomicInteger[3]
nGoodMigs[0] = new AtomicInteger()
nGoodMigs[1] = new AtomicInteger()
nGoodMigs[2] = new AtomicInteger()
nReadsInGoodMigs[0] = new AtomicInteger()
nReadsInGoodMigs[1] = new AtomicInteger()
nReadsInGoodMigs[2] = new AtomicInteger()

def getCoreSeq = { String seq, int offset ->
    int mid = seq.length() / 2
    seq.substring(mid - anchorRegion - offset, mid + anchorRegion + 1 - offset)
}

GParsPool.withPool THREADS, {
    migData[0].eachParallel { migEntry ->
        String umi = migEntry.key
        Map<String, Integer> reads1 = migEntry.value
        int count = (int) reads1.values().sum()

        def migsToAssemble = [reads1]
        if (bothReads) // only for 1,1
            migsToAssemble.add(migData[1].get(umi))
        def assembledReads = new String[4]

        int nMigsCurrent = nMigs.incrementAndGet()
        int nCollisionsCurrent = nCollisions.get()
        int nReadsInMigsCurrent = nReadsInMigs.addAndGet(count)
        int[] nGoodMigsCurrent = new int[3], nReadsInGoodMigsCurrent = new int[3]
        nGoodMigsCurrent[0] = nGoodMigs[0].get()
        nReadsInGoodMigsCurrent[0] = nReadsInGoodMigs[0].get()
        nGoodMigsCurrent[1] = nGoodMigs[1].get()
        nReadsInGoodMigsCurrent[1] = nReadsInGoodMigs[1].get()
        nGoodMigsCurrent[2] = nGoodMigs[2].get()
        nReadsInGoodMigsCurrent[2] = nReadsInGoodMigs[2].get()

        // Search for collisions
        boolean noCollision = true
        if (filterCollisions) {
            // A standard hash-based 1-loop single-mm search..
            char[] umiCharArray = umi.toCharArray()
            char oldChar
            for (int i = 0; i < umiCharArray.length; i++) {
                oldChar = umiCharArray[i]
                for (int j = 0; j < 4; j++) {
                    char nt = Util.code2nt(j)
                    if (nt != oldChar) {
                        umiCharArray[i] = nt
                        String otherUmi = new String(umiCharArray)

                        Map<String, Integer> otherReads1 = migData[0].get(otherUmi)
                        if (otherReads1 != null) {
                            int otherCount = (int) otherReads1.values().sum()
                            if (count / (double) otherCount < collisionRatioThreshold) {
                                noCollision = false
                                nCollisionsCurrent = nCollisions.incrementAndGet()
                                break
                            }
                        }
                    }
                }
                umiCharArray[i] = oldChar
            }
        }

        // Do assembly
        if (noCollision && count >= minCount) {
            migsToAssemble.eachWithIndex { Map<String, Integer> mig, int ind ->
                // Step 1: collect core regions with different offsets to determine most frequent one
                def coreSeqMap = new HashMap<String, int[]>()
                mig.each { Map.Entry<String, Integer> read ->
                    for (int offset = -offsetRange; offset <= offsetRange; offset++) {
                        String coreSeq = getCoreSeq(read.key, offset)
                        int[] coreSeqData = coreSeqMap.get(coreSeq)
                        if (coreSeqData == null)
                            coreSeqMap.put(coreSeq, coreSeqData = new int[2])
                        coreSeqData[0] += read.value
                        coreSeqData[1] += Math.abs(offset)
                    }
                }

                String bestCoreSeq = ""
                int[] bestCoreData = new int[2]

                coreSeqMap.each {
                    if (it.value[0] > bestCoreData[0] ||
                            (it.value[0] == bestCoreData[0] && it.value[1] < bestCoreData[1])) {
                        bestCoreSeq = it.key
                        bestCoreData = it.value
                    }
                }

                int maxX = Integer.MIN_VALUE, maxY = Integer.MIN_VALUE

                def migAlignmentData = new HashMap<String, int[]>()
                // Step 2: For all reads find optimal position against the core & append to pwm
                mig.each { Map.Entry<String, Integer> read ->
                    // 2.1 Determine best offset vs core
                    def bestOffset = 0, bestOffsetMMs = anchorRegion
                    for (int offset = -offsetRange; offset <= offsetRange; offset++) {
                        int offsetMMs = 0
                        String coreSeq = getCoreSeq(read.key, offset)

                        if (coreSeq == bestCoreSeq) {
                            bestOffset = offset
                            bestOffsetMMs = 0
                            break       // keep match
                        } else {
                            for (int j = 0; j < coreSeq.length(); j++)
                                if (coreSeq.charAt(j) != bestCoreSeq.charAt(j))
                                    offsetMMs++

                            if (offsetMMs < bestOffsetMMs) {
                                bestOffsetMMs = offsetMMs
                                bestOffset = offset
                            }
                        }
                    }

                    // 2.2 Keep if more than 'maxMMs' per 'anchorRegion'
                    if (bestOffsetMMs <= maxMMs) {
                        int l = read.key.length(), mid = l / 2
                        int x = mid - bestOffset, y = l - x
                        maxX = Math.max(maxX, x)
                        maxY = Math.max(maxY, y)

                        int[] data = new int[3]
                        data[0] = read.value
                        data[1] = x
                        data[2] = y
                        migAlignmentData.put(read.key, data)
                    } else
                        count -= read.value // drop it
                }

                // Still good?
                if (count >= minCount) {
                    // Step 3.1: Select region to construct PWM, append reads to PWM
                    int pwmLen = maxY + maxX
                    int[][] pwm = new int[pwmLen][4]

                    String detailInfo = "@" + umi + "\n" // details header

                    migAlignmentData.each {
                        int redundancyCount = it.value[0], x = it.value[1], y = it.value[2]
                        String seqRegion = 'N' * (maxX - x) + it.key + 'N' * (maxY - y)

                        for (int i = 0; i < pwmLen; i++) {
                            if (seqRegion.charAt(i) == 'N')
                                for (int j = 0; j < 4; j++)
                                    pwm[i][j] += redundancyCount / 4
                            else
                                pwm[i][Util.nt2code(seqRegion.charAt(i))] += redundancyCount
                        }

                        if (alignmentFilePrefix) // detailes - aligned read
                            detailInfo += seqRegion + "\t" + redundancyCount + "\n"
                    }

                    // Step 3.2: Calculate new quality
                    def consensus = new StringBuilder(), qual = new StringBuilder()
                    for (int i = 0; i < pwmLen; i++) {
                        int mostFreqLetter = 0, maxLetterFreq = 0
                        for (int j = 0; j < 4; j++) {
                            if (maxLetterFreq < pwm[i][j]) {
                                maxLetterFreq = pwm[i][j]
                                mostFreqLetter = j
                            }
                        }
                        consensus.append(Util.code2nt(mostFreqLetter))
                        qual.append(Util.symbolFromQual(Math.max(2, (int) ((maxLetterFreq / count - 0.25) / 0.75 * 40.0))))
                    }
                    assembledReads[ind] = "@MIG UMI:$umi:$count\n${consensus.toString()}\n+\n${qual.toString()}".toString()

                    if (alignmentFilePrefix)
                        assembledReads[ind + 2] = detailInfo

                    nGoodMigsCurrent[ind] = nGoodMigs[ind].incrementAndGet()
                    nReadsInGoodMigsCurrent[ind] = nReadsInGoodMigs[ind].addAndGet(count)
                }
            }
        }

        if ((!assemblyIndices[0] || assembledReads[0] != null) &&
                (!assemblyIndices[1] || assembledReads[0] != null) &&
                (!(assemblyIndices[0] && assemblyIndices[1]) || assembledReads[1] != null)) {
            writeQueue.put(assembledReads)
            nGoodMigsCurrent[2] = nGoodMigs[2].incrementAndGet()
            nReadsInGoodMigsCurrent[2] = nReadsInGoodMigs[2].addAndGet(count)
        }


        if (nMigsCurrent % 10000 == 0)
            println "[${new Date()} $scriptName] Processed $nMigsCurrent MIGs, $nReadsInMigsCurrent reads total, " +
                    "$nCollisionsCurrent collisions detected, assembled so far: " +
                    "$suffix ${nGoodMigsCurrent[0]} MIGs, ${nReadsInGoodMigsCurrent[0]} reads" +
                    (bothReads ? "; R2 ${nGoodMigsCurrent[1]} MIGs, ${nReadsInGoodMigsCurrent[1]} reads" : "") +
                    (bothReads ? "; Overall ${nGoodMigsCurrent[2]} MIGs, ${nReadsInGoodMigsCurrent[2]} reads" : "")
    }

}

println "[${new Date()} $scriptName] Processed ${nMigs.get()} MIGs, ${nReadsInMigs.get()} reads total, " +
        "${nCollisions.get()} collisions detected, assembled so far: " +
        "$suffix ${nGoodMigs[0].get()} MIGs, ${nReadsInGoodMigs[0].get()} reads" +
        (bothReads ? "; R2 ${nGoodMigs[1].get()} MIGs, ${nReadsInGoodMigs[1].get()} reads" : "") +
        (bothReads ? "; Overall ${nGoodMigs[2].get()} MIGs, ${nReadsInGoodMigs[2].get()} reads" : "")

writeQueue.put(new String[0])
writeThread.join()

if (logFileName != null) {
    if (!(new File(logFileName).exists())) {
        new File(logFileName).withPrintWriter { pw ->
            pw.println("#SAMPLE_ID\tSAMPLE_TYPE\tASSEMBLE_FASTQ1\tASSEMBLE_FASTQ2\tPAIRED_MASK\t" +
                    "MIG_COUNT_THRESHOLD\tMIGS_TOTAL\tREADS_TOTAL\tMIGS_GOOD\tREADS_GOOD")
        }
    }
    new File(logFileName).withWriterAppend { writer ->
        writer.println([outputFilePrefix =~ /[^\/]*$/, paired ? "paired" : (overlapped ? "overlapped" : "unpaired"),
                        outputFileName1, outputFileName2,
                        assemblyIndices.collect { it ? 0 : 1 }.join(":"),
                        minCount, nMigs.get(), nReadsInMigs.get(),
                        nGoodMigs[2].get(), nReadsInGoodMigs[2].get()].join("\t"))
    }
}
println "[${new Date()} $scriptName] Finished"