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

package migec.pipeline

@Grab(group = 'org.codehaus.gpars', module = 'gpars', version = '1.0.0')

import groovyx.gpars.GParsPool

import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicInteger
import java.util.regex.Pattern
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream


//========================
//          CLI
//========================
def mode = "0:1", mc = "10", p_c_r = "0.1"
def cli = new CliBuilder(usage: 'groovy [options] Assemble R1.fastq[.gz] [R2.fastq[.gz] or -] output_prefix ' +
        '[assembly_log, to append]')
cli.q(args: 1, argName: 'read quality (phred)',
        "barcode region quality threshold. Default: 15")
cli.m(longOpt: "assembly-mode", args: 1, argName: 'assembly mode in format X:X',
        "Identifier(s) of read(s) to assemble. Default: \"$mode\". In case of \"0:0\" will try to overlap reads.")
cli.p(args: 1,
        "number of threads to use. Default: all available processors")
cli.c("compressed output")

cli._(longOpt: 'alignment-file-prefix', args: 1, argName: 'string',
        "File name prefix to output multiple alignments generated during assembly, for \"migec.control.BacktrackSequence\"")
cli._(longOpt: 'min-count', args: 1, argName: 'integer',
        "Min number of reads in MIG. Should be set according to 'Histogram.groovy' output. Default: $mc")
cli.f(longOpt: 'filter-collisions',
        "Collision filtering. Should be set if collisions (1-mismatch erroneous UMI sequence variants) are observed in 'Histogram.groovy' output")
cli._(longOpt: 'collision-ratio', args: 1, argName: 'double, < 1.0',
        "Min parent-to-child MIG size ratio for collision filtering. Default value: $p_c_r")
cli._(longOpt: 'assembly-offset', args: 1, argName: 'integer',
        'Assembly offset range. Default: 3')
cli._(longOpt: 'assembly-mismatches', args: 1, argName: 'integer',
        'Assembly max mismatches. Default: 5')
cli._(longOpt: 'assembly-anchor', args: 1, argName: 'integer',
        'Assembly anchor region half size. Default: 10')

def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 3) {
    cli.usage()
    return
}

//========================
//         PARAMS
//========================
def scriptName = getClass().canonicalName

// Parameters
boolean compressed = opt.c, filterCollisions = opt.f
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
int umiQualThreshold = opt.q ?: 15, minCount = Integer.parseInt(opt."min-count" ?: mc)
double collisionRatioThreshold = Double.parseDouble(opt.'collision-ratio' ?: p_c_r)

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
def assemblyIndices = paired ? (opt.m ?: mode).split(":").collect { Integer.parseInt(it) > 0 } : [true, false]
boolean overlap = false, bothReads = assemblyIndices[0] && assemblyIndices[1]
if (!assemblyIndices.any()) {
    assemblyIndices = [true, false]
    overlap = true
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

//========================
//      MISC UTILS
//========================
def qualFromSymbol = { char symbol ->
    (int) symbol - 33
}

def symbolFromQual = { int qual ->
    qual = qual < 2 ? 2 : qual
    qual = qual > 40 ? 40 : qual
    (char) (qual + 33)
}

def code2nt = { int code ->
    switch (code) {
        case 0:
            return 'A'
        case 1:
            return 'T'
        case 2:
            return 'G'
        case 3:
            return 'C'
    }
}

def nt2code = { char symbol ->
    switch (symbol) {
        case 'A':
            return 0
        case 'T':
            return 1
        case 'G':
            return 2
        case 'C':
            return 3
    }
}

def getReader = { String fname ->
    new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fname)) :
            new FileInputStream(fname)))
}

def getWriter = { String outfile ->
    if (compressed)
        outfile += ".gz"
    new BufferedWriter(new OutputStreamWriter(compressed ?
            new GZIPOutputStream(new FileOutputStream(outfile)) : new FileOutputStream(outfile)))
}

def getUmi = { String header ->
    def splitHeader = header.split(" ")
    def umiEntry = splitHeader.find { it.startsWith("UMI:") }
    if (umiEntry == null) {
        println "[${new Date()} $scriptName] Error: no UMI header in input. Terminating"
        System.exit(-1)
    }
    String umi = umiEntry.split(":")[1]
    for (int i = umiEntry.length() - umi.length(); i < umiEntry.length(); i++)
        if (qualFromSymbol(umiEntry.charAt(i)) < umiQualThreshold)   // quality can contain :
            return null
    umi
}

// Overlap, select top quality nts for overlap zone
int maxOffset = 5, mmOverlapSz = 10
int k = 5
def overlapReads = { String r1, String r2, String q1, String q2 ->
    for (int i = 0; i < maxOffset; i++) {
        def kmer = r2.substring(i, i + k)
        def pattern = Pattern.compile(kmer)
        def matcher = pattern.matcher(r1)
        // Find last match
        int position
        while (matcher.find()) {
            position = matcher.start()
            if (position >= 0) {
                // Start fuzzy align
                int nmm = 0
                for (int j = 0; j < mmOverlapSz; j++) {
                    def posInR1 = position + k + j, posInR2 = i + k + j
                    if (posInR1 + 1 > r1.length())
                        break  // went to end of r1, all fine
                    if (r1.charAt(posInR1) != r2.charAt(posInR2)) {
                        if (++nmm > 1)
                            break  // two consequent mismatches
                    } else {
                        nmm = 0 // zero counter
                    }
                }
                if (nmm < 2) {
                    // take best qual nts
                    def seq = new StringBuilder(r1.substring(0, position))
                    int pos2 = i - 1
                    for (int j = position; j < r1.length(); j++) {
                        pos2++
                        if (pos2 == r2.length())
                            break // should not happen
                        seq.append(qualFromSymbol(q1.charAt(j)) > qualFromSymbol(q2.charAt(pos2)) ?
                                r1.charAt(j) : r2.charAt(pos2))
                    }
                    for (int j = pos2 + 1; j < r2.length(); j++)
                        seq.append(r2.charAt(j)) // fill the remainder
                    return seq.toString() // passed test
                }
            }
        }
    }
    return "" // failed
}

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
println "[${new Date()} $scriptName] Pre-loading data fo $fastq1, $fastq2.."
def reader1 = getReader(fastq1), reader2 = paired ? getReader(fastq2) : null
String header1, seq1, qual1
String seq2 = "", qual2 = ""
int MIN_READ_SZ = Math.max(k, 2 * anchorRegion + 1 + offsetRange)
while ((header1 = reader1.readLine()) != null) {
    seq1 = reader1.readLine()
    reader1.readLine()
    qual1 = reader1.readLine()

    if (paired) {
        reader2.readLine() // Skip header of read 2
        seq2 = reader2.readLine()
        reader2.readLine()
        qual2 = reader2.readLine()
    }

    if (seq1.length() > MIN_READ_SZ && (!paired || seq2.length() > MIN_READ_SZ)) {
        def umi = getUmi(header1)
        if (overlap)
            seq1 = overlapReads(seq1, seq2, qual1, qual2)

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
                    (overlap ? "successfully overlapped and " : "") + "found good UMI header for $nGoodReads reads, " +
                    "unique UMIs so far ${migData[0].size()}"
    }
}
println "[${new Date()} $scriptName] Processed $nReads reads, " +
        (overlap ? "successfully overlapped and " : "") + "found good UMI header for $nGoodReads reads, " +
        "unique UMIs ${migData[0].size()}"

//=================================
//   PERFORM ASSEMBLY
//=================================
def writeQueue = new LinkedBlockingQueue<String[]>(2048)
def suffix = overlap ? "RO" : (assemblyIndices[0] ? "R1" : "R2")
def writer1 = getWriter("${outputFilePrefix}_${suffix}.fastq"),
    writer2 = bothReads ? getWriter("${outputFilePrefix}_R2.fastq") : null
def detailsWriter1 = alignmentFilePrefix ? getWriter("${alignmentFilePrefix}_${suffix}.asm") : null,
    detailsWriter2 = bothReads && alignmentFilePrefix ? getWriter("${alignmentFilePrefix}_${suffix}.asm") : null

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
                    char nt = code2nt(j)
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
                                pwm[i][nt2code(seqRegion.charAt(i))] += redundancyCount
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
                        consensus.append(code2nt(mostFreqLetter))
                        qual.append(symbolFromQual(Math.max(2, (int) ((maxLetterFreq / count - 0.25) / 0.75 * 40.0))))
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
            pw.println("SAMPLE_ID\tMIG_COUNT_THRESHOLD\tMIGS_TOTAL\tREADS_TOTAL\tMIGS_GOOD\tREADS_GOOD")
        }
    }
    new File(logFileName).withWriterAppend { writer ->
        writer.println([outputFilePrefix, minCount,
                        nMigs.get(), nReadsInMigs.get(),
                        nGoodMigs[2].get(), nReadsInGoodMigs[2].get()].collect().join("\t"))
    }
}
println "[${new Date()} $scriptName] Finished"