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
def DEFAULT_ASSEMBLE_MASK = "1:1", DEFAULT_MIN_COUNT = "10", DEFAULT_PARENT_CHILD_RATIO = "0.1"
def cli = new CliBuilder(usage:
        'Assemble [options] R1.fastq[.gz] [R2.fastq[.gz] or -] output_dir/')
cli.q(args: 1, argName: 'read quality (phred)',
        "barcode region quality threshold. Default: $Util.DEFAULT_UMI_QUAL_THRESHOLD")
cli._(longOpt: 'assembly-mask', args: 1, argName: 'X:Y, X=0/1, Y=0/1',
        "Mask for read(s) in pair that should be assembled. Default: \"$DEFAULT_ASSEMBLE_MASK\".")
cli.p(args: 1,
        "number of threads to use. Default: all available processors")
cli.c("compressed output")
cli._(longOpt: 'log-file', args: 1, argName: 'fileName', "File to output assembly log")
cli._(longOpt: 'log-overwrite', "Overwrites provided log file")
cli._(longOpt: 'log-sample-name', "Sample name to use in log [default = N/A]")
cli._(longOpt: 'log-sample-type', "Sample type to use in log, i.e. unpaired, paired and overlapped [default = N/A]")
cli._(longOpt: 'alignment-details',
        "Output multiple alignments generated during assembly as .asm files, " +
                "for \"BacktrackSequence\"")
cli.m(longOpt: 'min-count', args: 1, argName: 'integer',
        "Minimal number of reads in MIG. Should be set according to 'Histogram.groovy' output. " +
                "Default: $DEFAULT_MIN_COUNT")
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
boolean compressed = opt.c, filterCollisions = opt.'filter-collisions', overwriteLog = opt.'log-overwrite'
String sampleName = opt.'log-sample-name' ?: "N/A", sampleType = opt.'log-sample-type' ?: "N/A"
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
byte umiQualThreshold = opt.q ? Byte.parseByte(opt.q) : Util.DEFAULT_UMI_QUAL_THRESHOLD
int minMigSize = Integer.parseInt(opt.'min-count' ?: DEFAULT_MIN_COUNT)
double collisionRatioThreshold = Double.parseDouble(opt.'collision-ratio' ?: DEFAULT_PARENT_CHILD_RATIO)

// I/O
def inputFileName1 = opt.arguments()[0],
    inputFileName2 = opt.arguments()[1],
    outputDir = opt.arguments()[2],
    outputFilePrefix1 = '-', outputFilePrefix2 = '-'

String logFileName = opt.'log-file' ?: null

if (!(inputFileName1.endsWith(".fastq") || inputFileName1.endsWith(".fastq.gz"))) {
    println "[ERROR] Bad file extension $inputFileName1. Either .fastq or .fastq.gz should be provided as R1 file."
    System.exit(-1)
} else {
    outputFilePrefix1 = Util.getFastqPrefix(inputFileName1) + ".t" + minMigSize + (filterCollisions ? ".cf" : "")
}

if (inputFileName2 != "-") {
    if (!(inputFileName2.endsWith(".fastq") || inputFileName2.endsWith(".fastq.gz"))) {
        println "[ERROR] Bad file extension $inputFileName2. Either .fastq, .fastq.gz or \'-\' should be provided as R2 file."
        System.exit(-1)
    } else {
        outputFilePrefix2 = Util.getFastqPrefix(inputFileName2) + ".t" + minMigSize + (filterCollisions ? ".cf" : "")
    }
}

// Assembly parameters
def offsetRange = Integer.parseInt(opt.'assembly-offset' ?: '5'),
    maxMMs = Integer.parseInt(opt.'assembly-mismatches' ?: '5'),
    anchorRegion = Integer.parseInt(opt.'assembly-anchor' ?: '10')

// I/O parameters
boolean paired = inputFileName2 != "-"
def assemblyIndices = [true, false]
boolean bothReads = false
if (paired) {
    def assemblyMask = (opt.'assembly-mask' ?: DEFAULT_ASSEMBLE_MASK).toString()
    if (!Util.MASKS.any { it == assemblyMask }) {
        println "[ERROR] Bad mask $assemblyMask. Allowed masks for paired-end mode are ${Util.MASKS.join(", ")}"
        System.exit(-1)
    }
    if (assemblyMask == "0:0") {
        println "[WARNING] Blank mask specified for paired-end data, skipping"
        System.exit(0)
    }
    assemblyIndices = (opt.'assembly-mask' ?: DEFAULT_ASSEMBLE_MASK).split(":").collect { Integer.parseInt(it) > 0 }
    bothReads = assemblyIndices[0] && assemblyIndices[1]
}

String outputFileNameNoExt1 = (!paired || assemblyIndices[0]) ? outputFilePrefix1 : outputFilePrefix2,
       outputFileNameNoExt2 = outputFilePrefix2

// Misc output
boolean alignmentDetails = opt.'alignment-details'

new File(outputDir).mkdirs()

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
println "[${new Date()} $scriptName] Pre-loading data for $inputFileName1, $inputFileName2.."
def reader1 = Util.getReader(inputFileName1), reader2 = paired ? Util.getReader(inputFileName2) : null
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
                putData(0, umi, seq2) // put all to 1st file if mask=0,1
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

def writer1 = Util.getWriter(outputDir + "/" + outputFileNameNoExt1 + ".fastq", compressed),
    writer2 = bothReads ? Util.getWriter(outputDir + "/" + outputFileNameNoExt2 + ".fastq", compressed) : null
def detailsWriter1 = alignmentDetails ? Util.getWriter(outputDir + "/" + outputFileNameNoExt1 + ".asm", compressed) : null,
    detailsWriter2 = bothReads && alignmentDetails ? Util.getWriter(outputDir + "/" + outputFileNameNoExt2 + ".asm", compressed) : null

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

        if (alignmentDetails)
            detailsWriter1.writeLine(result[2])

        if (alignmentDetails && bothReads)
            detailsWriter2.writeLine(result[3])
    }

    writer1.close()

    if (bothReads)
        writer2.close()

    if (alignmentDetails)
        detailsWriter1.close()

    if (alignmentDetails && bothReads)
        detailsWriter2.close()
} as Runnable)
writeThread.start()

def nMigs = new AtomicInteger(),
        nReadsInMigs = new AtomicInteger(), nDroppedReads = new AtomicInteger(),
    nCollisions = new AtomicInteger()
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
        if (noCollision && count >= minMigSize) {
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
                    } else {
                        count -= read.value // drop it
                        nDroppedReads.incrementAndGet()
                    }
                }

                // Still good?
                if (count >= minMigSize) {
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

                        if (alignmentDetails) // detailes - aligned read
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

                    if (alignmentDetails)
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
                    "$outputFileNameNoExt1 ${nGoodMigsCurrent[0]} MIGs, ${nReadsInGoodMigsCurrent[0]} reads" +
                    (bothReads ? "; $outputFileNameNoExt2 ${nGoodMigsCurrent[1]} MIGs, ${nReadsInGoodMigsCurrent[1]} reads" : "") +
                    (bothReads ? "; Overall ${nGoodMigsCurrent[2]} MIGs, ${nReadsInGoodMigsCurrent[2]} reads" : "")
    }

}

println "[${new Date()} $scriptName] Processed ${nMigs.get()} MIGs, ${nReadsInMigs.get()} reads total, " +
        "${nCollisions.get()} collisions detected, assembled so far: " +
        "$outputFileNameNoExt1 ${nGoodMigs[0].get()} MIGs, ${nReadsInGoodMigs[0].get()} reads" +
        (bothReads ? "; $outputFileNameNoExt2 ${nGoodMigs[1].get()} MIGs, ${nReadsInGoodMigs[1].get()} reads" : "") +
        (bothReads ? "; Overall ${nGoodMigs[2].get()} MIGs, ${nReadsInGoodMigs[2].get()} reads" : "")

writeQueue.put(new String[0])
writeThread.join()

println "[${new Date()} $scriptName] Finished"

def logLine = [assemblyIndices[0] ? new File(inputFileName1).absolutePath : '-',
               assemblyIndices[1] ? new File(inputFileName2).absolutePath : '-',
               assemblyIndices[0] ? new File(outputDir).absolutePath + '/' +
                       outputFilePrefix1 + ".fastq${compressed ? ".gz" : ""}" : '-',
               assemblyIndices[1] ? new File(outputDir).absolutePath + '/' +
                       outputFilePrefix2 + ".fastq${compressed ? ".gz" : ""}" : '-',

               minMigSize,

               nGoodMigs[0].get(), nGoodMigs[1].get(), nGoodMigs[2].get(), nMigs.get(),

               nReadsInGoodMigs[0].get(), nReadsInGoodMigs[1].get(), nReadsInGoodMigs[2].get(), nReadsInMigs.get(),
               nDroppedReads.get()].join("\t")


if (logFileName) {
    def logFile = new File(logFileName)
    if (logFile.exists()) {
        if (overwriteLog)
            logFile.delete()
    } else {
        logFile.absoluteFile.parentFile.mkdirs()
        logFile.withPrintWriter { pw ->
            pw.println(Util.ASSEMBLE_LOG_HEADER)
        }
    }

    logFile.withWriterAppend { writer ->
        writer.println("$sampleName\t$sampleType\t" + logLine)
    }
}

return logLine