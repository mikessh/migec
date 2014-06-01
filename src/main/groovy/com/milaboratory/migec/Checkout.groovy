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

import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicLong
import java.util.regex.Pattern
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

def mm = "15:0.2:0.05", rcm = "0:1", mtrim = "10"
def cli = new CliBuilder(usage:
        'Checkout [options] barcode_file R1.fastq[.gz] [R2.fastq[.gz] or -] [path/to/out/]')
cli.o('Oriented reads, so master barcode has to be in R1. ' +
        'Default: scans both reads (if pair-end) for master barcode')
cli.u('Save UMI region specified by capital N\'s in barcode sequence to the header')
cli.t('Trim barcode sequences. Also trims boundary sequence between barcode and read start/end if it is short.')
cli.e('Remove template-switching trace, ^T{0,3}G{3,7}. Only works in compination with -t')
cli._(longOpt: 'max-trim-nts', 'Maximal number of nucleotides between barcode and read start/end that could be trimmed')
cli.r(args: 1, "RC mask to apply after master-slave is determined, e.g. will RC slave read if set to 0:1. " +
        "Default: $rcm")
cli.p(args: 1, 'Number of threads. Default: all available processors.')
cli._(args: 1, longOpt: 'first', 'Number of reads to try. If <0 will process all reads. [default = -1]')
cli._(longOpt: 'overlap', 'Will try to overlap paired reads.')
cli._(longOpt: 'rc-barcodes', 'Also search for reverse-complement barcodes.')
cli.m(args: 1,
        argName: 'LQ:E0:E1', "Low quality threshold : Ratio of low-quality errors : Ratio of high-quality errors. " +
        "Used to calculate the maximum number of allowed mismatches when performing alignment to full barcode, " +
        "ratios are scaled to full barcode length. Default: $mm")
cli.c('Compressed output')
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def scriptName = getClass().canonicalName
boolean compressed = opt.c, oriented = opt.o, extractUMI = opt.u, addRevComplBc = opt."rc-barcodes",
        trimBc = opt.t, removeTS = opt.e, overlap = opt.'overlap'
int trimSizeThreshold = (opt.'max-trim-nts' ?: mtrim).toInteger()
int firstReadsToTake = (opt.'first' ?: "-1").toInteger()
def rcReadMask = (opt.r ?: rcm).split(":").collect { Integer.parseInt(it) > 0 }
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
def bcfile = opt.arguments()[0],
    fastq1 = opt.arguments()[1],
    fastq2 = opt.arguments()[2],
    out = opt.arguments()[3] ?: "."
def mmData = (opt.m ?: mm).split(":").collect { Double.parseDouble(it) }
def paired = fastq2 != '-'

if (overlap && !paired) {
    println "ERROR Overlap requested for unpaired reads"
    System.exit(-1)
}

//========================
//          MISC
//========================
def complements = ['A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'y': 'r',
                   'r': 'y', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n']

def redundancy = ['A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': '[ATGC]', 'R': '[AG]', 'Y': '[CT]', 'M': '[AC]',
                  'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'V': '[ACG]', 'D': '[AGT]', 'H': '[ACT]', 'B': '[CGT]']

def redundancy2 = new HashMap<Character, Set<Character>>([
        ((Character) 'A'): new HashSet<Character>([(Character) 'A']),
        ((Character) 'T'): new HashSet<Character>([(Character) 'T']),
        ((Character) 'G'): new HashSet<Character>([(Character) 'G']),
        ((Character) 'C'): new HashSet<Character>([(Character) 'C']),
        ((Character) 'N'): new HashSet<Character>([(Character) 'A', (Character) 'T', (Character) 'G', (Character) 'C']),
        ((Character) 'R'): new HashSet<Character>([(Character) 'A', (Character) 'G']),
        ((Character) 'Y'): new HashSet<Character>([(Character) 'C', (Character) 'T']),
        ((Character) 'M'): new HashSet<Character>([(Character) 'A', (Character) 'C']),
        ((Character) 'S'): new HashSet<Character>([(Character) 'G', (Character) 'C']),
        ((Character) 'W'): new HashSet<Character>([(Character) 'A', (Character) 'T']),
        ((Character) 'K'): new HashSet<Character>([(Character) 'G', (Character) 'T']),
        ((Character) 'V'): new HashSet<Character>([(Character) 'A', (Character) 'C', (Character) 'G']),
        ((Character) 'D'): new HashSet<Character>([(Character) 'A', (Character) 'G', (Character) 'T']),
        ((Character) 'H'): new HashSet<Character>([(Character) 'A', (Character) 'C', (Character) 'T']),
        ((Character) 'B'): new HashSet<Character>([(Character) 'C', (Character) 'G', (Character) 'T'])])

def qualFromSymbol = { char symbol ->
    (int) symbol - 33
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

def revCompl = { String seq ->
    def chars = seq.reverse().toCharArray()
    for (int i = 0; i < chars.length; i++) {
        if (chars[i] == (char) 'A')
            chars[i] = (char) 'T'
        else if (chars[i] == (char) 'T')
            chars[i] = (char) 'A'
        else if (chars[i] == (char) 'G')
            chars[i] = (char) 'C'
        else if (chars[i] == (char) 'C')
            chars[i] = (char) 'G'
        else if (chars[i] == (char) 'N')
            chars[i] = (char) 'N'
        else if (chars[i] == (char) 'a')
            chars[i] = (char) 't'
        else if (chars[i] == (char) 't')
            chars[i] = (char) 'a'
        else if (chars[i] == (char) 'g')
            chars[i] = (char) 'c'
        else if (chars[i] == (char) 'c')
            chars[i] = (char) 'g'
        else
            chars[i] = (char) 'N'
    }
    return new String(chars)
}

//========================
//       BARCODE DATA
//========================
def sampleIds = new ArrayList<String>()
def barcodeList = new List<String>[2]
def seeds = new List<Pattern>[2]
def umiPositions = new List<List<Integer>>[2]
def maxGoodMMs = new List<Integer>[2], maxBadMMs = new List<Integer>[2]

// Init
barcodeList[0] = new ArrayList<String>()
barcodeList[1] = new ArrayList<String>()
seeds[0] = new ArrayList<Pattern>()
seeds[1] = new ArrayList<Pattern>()
umiPositions[0] = new ArrayList<List<Integer>>()
umiPositions[1] = new ArrayList<List<Integer>>()
maxBadMMs[0] = new ArrayList<Integer>()
maxBadMMs[1] = new ArrayList<Integer>()
maxGoodMMs[0] = new ArrayList<Integer>()
maxGoodMMs[1] = new ArrayList<Integer>()

// Subs for parsing barcode data
def getSeed = { String barcode ->
    Pattern.compile(barcode.replaceAll(/./) { redundancy."$it" ?: '[ATGC]' }) // redundancy has only uppecase
}

def parseUmiPositions = { String barcode ->
    def posList = new ArrayList<Integer>()
    for (int i = 0; i < barcode.length(); i++)
        if (barcode.charAt(i) == 'N')
            posList.add(i)
    posList
}

def addBarcode = { String barcode, int slave ->
    barcodeList[slave].add(barcode.toUpperCase()) // important - no lcase masking later
    seeds[slave].add(getSeed(barcode))
    umiPositions[slave].add(parseUmiPositions(barcode))
    maxBadMMs[slave].add((int) (barcode.length() * mmData[1]))
    maxGoodMMs[slave].add((int) (barcode.length() * mmData[2]))
}

// Load barcode data from file
new File(bcfile).splitEachLine("[\t ]") { sl ->
    if (!sl[0].startsWith("#")) {
        sampleIds.add(sl[0])
        if (addRevComplBc)
            sampleIds.add(sl[0])
        if (sl[1] != null && (sl[1] = sl[1].trim()).length() > 0 &&
                sl[1].toCharArray().every { complements.keySet().contains((String) it) }) {
            addBarcode(sl[1], 0)
            if (addRevComplBc)
                addBarcode(revCompl(sl[1]), 0)
        } else {
            println "ERROR: Bad barcode ${sl[1]}. Terminating"
            System.exit(-1)
        }
        if (sl[2] != null && (sl[2] = sl[2].trim()).length() > 0 &&
                sl[2].toCharArray().every { complements.keySet().contains((String) it) }) {
            addBarcode(sl[2], 1)
            if (addRevComplBc)
                addBarcode(revCompl(sl[2]), 1)
        } else {
            addBarcode('n', 1)
            if (addRevComplBc)
                addBarcode('n', 1)
        }
    }
}

println "[${new Date()} $scriptName] Loaded barcode data for " +
        "${addRevComplBc ? (sampleIds.size() / 2) : sampleIds.size()} sample(s)"

//========================
//  BARCODE SEARCH UTILS
//========================
def findAllMatches = { String seq, Pattern seed ->  // all matches for seed
    def positions = new ArrayList<Integer>()
    def matcher = seed.matcher(seq)
    while (matcher.find()) {
        positions.add(matcher.start())
    }
    positions
}
def lq = mmData[0]
def hasFuzzyMatch = { String barcode, String seq, String qual, int from, int bcIndex, int slave -> // fuzzy match from seed match
    int nBadMMs = 0, nGoodMMs = 0
    int goodThreshold = maxGoodMMs[slave][bcIndex], badThreshold = maxBadMMs[slave][bcIndex]
    for (int i = 0; i < barcode.size(); i++) {
        if (!(redundancy2[(Character) barcode.charAt(i)].contains((Character) seq.charAt(from + i)))) {
            if (qualFromSymbol(qual.charAt(i)) > lq)
                if (++nGoodMMs > goodThreshold)
                    break
                else if (++nBadMMs > badThreshold)
                    break
        }
    }
    nGoodMMs <= goodThreshold && nBadMMs <= badThreshold
}
def findMatch = { String barcode, Pattern seed, String seq, String qual, int bcIndex, int slave ->
    def seedOccurences = findAllMatches(seq, seed) // All seed occurences

    for (int i = 0; i < seedOccurences.size(); i++) {
        if (hasFuzzyMatch(barcode, seq, qual, seedOccurences[i], bcIndex, slave))  // Exhaustive search
            return seedOccurences[i] // till first best match
    }
    return -1
}

//========================
//  READ PROCESSING UTILS
//========================
// Overlap, select top quality nts for overlap zone
int maxOffset = 5, mmOverlapSz = 10
int K = 5
def overlapReads = { String r1, String r2, String q1, String q2 ->
    String[] result = null

    for (int i = 0; i < maxOffset; i++) {
        def kmer = r2.substring(i, i + K)
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
                    def posInR1 = position + K + j, posInR2 = i + K + j
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
                    def seq = new StringBuilder(r1.substring(0, position)),
                        qual = new StringBuilder(q1.substring(0, position))

                    int pos2 = i - 1
                    for (int j = position; j < r1.length(); j++) {
                        pos2++
                        if (pos2 == r2.length())
                            break // should not happen

                        seq.append(qualFromSymbol(q1.charAt(j)) > qualFromSymbol(q2.charAt(pos2)) ?
                                r1.charAt(j) : r2.charAt(pos2))

                        qual.append(qualFromSymbol(q1.charAt(j)) > qualFromSymbol(q2.charAt(pos2)) ?
                                q1.charAt(j) : q2.charAt(pos2))
                    }
                    for (int j = pos2 + 1; j < r2.length(); j++) {
                        // fill the remainder
                        seq.append(r2.charAt(j))

                        qual.append(q2.charAt(j))
                    }

                    // report overlap
                    result = new String[2]
                    result[0] = seq.toString()
                    result[1] = qual.toString()

                    return result
                }
            }
        }
    }

    result // failed
}

def flip = { String[] x, int i, int j ->
    String tmp
    tmp = x[i]
    x[i] = x[j]
    x[j] = tmp
}
def counters = new HashMap<String, AtomicLong[]>()
def readCounter = new AtomicLong(), goodReadCounter = new AtomicLong(), masterFirstCounter = new AtomicLong()
def overlapCounter = new AtomicLong()
def wrapRead = { String[] readData, StringBuilder[] umiData, int readIndex, String sampleId, String masterSampleId,
                 boolean good ->
    String umiHeader = umiData[0].size() > 0 ? " UMI:${umiData[0].toString()}:${umiData[1].toString()}" : ""
    readData[0] = readData[0] + umiHeader

    if (paired) {
        readData[3] = readData[3] + umiHeader

        if (readIndex > 0) { // flip if master is R2
            flip(readData, 0, 3)
            flip(readData, 1, 4)
            flip(readData, 2, 5)
        }

        if (rcReadMask[1]) { // rc slave if needed
            readData[4] = revCompl(readData[4])
            readData[5] = readData[5].reverse()
        }
        readData[6] = sampleId
    } else readData[3] = sampleId

    if (rcReadMask[0]) {  // rc master if needed
        readData[1] = revCompl(readData[1])
        readData[2] = readData[2].reverse()
    }

    // Record stats and try overlap
    long nReads = readCounter.incrementAndGet(),
         nGoodReads = good ? goodReadCounter.incrementAndGet() : goodReadCounter.get(),
         nMasterFirst = (good && (readIndex > 0)) ? masterFirstCounter.incrementAndGet() : masterFirstCounter.get(),
         overlapCount = overlapCounter.get()

    if (good) {
        counters.get(sampleId)[0].incrementAndGet()
        if (paired) {
            counters.get(sampleId)[1].incrementAndGet()

            def overlapResult = overlapReads(readData[1], readData[4], readData[2], readData[5])

            if (overlapResult) {
                overlapCounter.incrementAndGet()
                readData[1] = overlapResult[0]
                readData[2] = overlapResult[1]
                readData[7] = "yes"

                counters.get(sampleId)[2].incrementAndGet()
            }
        }
    } else {
        counters.get(masterSampleId)[0].incrementAndGet() // where master was found
        if (paired)
            counters.get(sampleId)[1].incrementAndGet()
    }

    if (nReads % 25000 == 0)
        println "[${new Date()} $scriptName] Processed $nReads, " +
                "identified $nGoodReads (${((int) (10000 * (double) nGoodReads / (double) nReads)) / 100}%), " +
                (overlap ? "overlapped ${((int) (10000 * (double) overlapCount / (double) nGoodReads)) / 100}% of them, " : "") +
                "assymetry (master first): $nMasterFirst (${((int) (10000 * (double) nMasterFirst / (double) nGoodReads)) / 100}%)"

    readData
}

//========================
//          BODY
//========================
def readQueue = new LinkedBlockingQueue<String[]>(4096)  // SINGLE: h1 r1 q1 ''       PAIRED: h1 r1 q1 h2 r2 q2 ''       ''
def writeQueue = new LinkedBlockingQueue<String[]>(4096) // SINGLE: h1 r1 q1 sampleId PAIRED: h1 r1 q1 h2 r2 q2 sampleId overlapped

def reader1 = getReader(fastq1), reader2 = (fastq2 == '-') ? null : getReader(fastq2)
def writers = new HashMap<String, BufferedWriter[]>()
println "[${new Date()} $scriptName] Started processing for $fastq1, $fastq2"
new File(out).mkdir()
[sampleIds, 'undef-s', 'undef-m'].flatten().each { String sampleId ->
    def writerTrio = new BufferedWriter[3]

    writerTrio[0] = getWriter("$out/${sampleId}_R1.fastq")
    if (paired) {
        writerTrio[1] = getWriter("$out/${sampleId}_R2.fastq")
        if (overlap)
            writerTrio[2] = getWriter("$out/${sampleId}_R12.fastq")
    }

    writers.put(sampleId, writerTrio)

    def counterTrio = new AtomicLong[3]

    counterTrio[0] = new AtomicLong()
    if (paired) {
        counterTrio[1] = new AtomicLong()
        if (overlap)
            counterTrio[2] = new AtomicLong()
    }

    counters.put(sampleId, counterTrio)
}

int nProcessors = 2 * THREADS
def rDataSize = paired ? 8 : 4
int readsTaken = 0
def readThread = new Thread({  // Reading thread
    String header1
    while ((header1 = reader1.readLine()) != null) {
        def readData = new String[rDataSize]
        readData[0] = header1
        readData[1] = reader1.readLine()
        reader1.readLine()
        readData[2] = reader1.readLine()

        if (paired) {
            readData[3] = reader2.readLine()
            readData[4] = reader2.readLine()
            reader2.readLine()
            readData[5] = reader2.readLine()
        }
        readQueue.put(readData)

        if (firstReadsToTake >= 0 && ++readsTaken > firstReadsToTake)
            break
    }
    for (int k = 0; k < nProcessors; k++)
        readQueue.put(new String[0]) // empty read - finish. tell ALL processors to stop
} as Runnable)
readThread.start()

def pool = Executors.newFixedThreadPool(THREADS)
def futures = new ArrayList<Future>()

def tgPattern = /^T{0,3}G{3,7}/, gtPattern = /G{3,7}T{0,3}$/

for (int k = 0; k < nProcessors; k++) {
    futures.add(pool.submit(new Runnable() { // Processors
        @Override
        void run() {
            while (true) {
                def readData = readQueue.take()
                if (readData.length == 0)
                    break
                def umiData = new StringBuilder[2]
                umiData[0] = new StringBuilder()
                umiData[1] = new StringBuilder()

                // 1. Search master barcode
                int sampleIndex = -1, readIndex = -1, from = -1
                for (int i = 0; i < barcodeList[0].size(); i++) {
                    // R1
                    if ((from = findMatch(barcodeList[0][i], seeds[0][i],
                            readData[1], readData[2], i, 0)) >= 0) {
                        sampleIndex = i
                        readIndex = 0
                        break
                    }
                    // R2
                    if (!oriented && paired &&
                            ((from = findMatch(barcodeList[0][i], seeds[0][i],
                                    readData[4], readData[5], i, 0)) >= 0)) {
                        sampleIndex = i
                        readIndex = 3
                        break
                    }
                }

                if (from >= 0) {
                    // 1.2 Extract UMI from master
                    if (extractUMI) {
                        def seq = readData[readIndex + 1], qual = readData[readIndex + 2]
                        umiPositions[0][sampleIndex].each { int pos ->
                            umiData[0].append(seq.charAt(from + pos))
                            umiData[1].append(qual.charAt(from + pos))
                        }
                    }

                    // trim adapter
                    if (trimBc) {
                        int to = barcodeList[0][sampleIndex].size() + from
                        def seq = readData[readIndex + 1], qual = readData[readIndex + 2]
                        readData[readIndex + 1] = (from > trimSizeThreshold ? seq.substring(0, from) : "") +
                                ((seq.length() - to > trimSizeThreshold) ? seq.substring(to) : "")
                        readData[readIndex + 2] = (from > trimSizeThreshold ? qual.substring(0, from) : "") +
                                ((seq.length() - to > trimSizeThreshold) ? qual.substring(to) : "")
                        if (removeTS) {
                            seq = readData[readIndex + 1]
                            qual = readData[readIndex + 2]
                            def tsHit = (seq =~ tgPattern)
                            if (tsHit) {
                                readData[readIndex + 1] = seq.substring(tsHit.end())
                                readData[readIndex + 2] = qual.substring(tsHit.end())
                            }
                        }
                    }

                    // 2. Look for slave
                    if (paired) {
                        boolean slaveFound = barcodeList[1][sampleIndex] == 'N' // slave could be skipped
                        if (!slaveFound) {
                            // Look in other read
                            from = findMatch(barcodeList[1][sampleIndex], seeds[1][sampleIndex],
                                    readData[4 - readIndex], readData[5 - readIndex], sampleIndex, 1)

                            slaveFound = from >= 0

                            if (slaveFound) {
                                def seq = readData[4 - readIndex], qual = readData[5 - readIndex]

                                // 2.1 Extract UMI from slave
                                if (extractUMI) {
                                    umiPositions[1][sampleIndex].each { int pos ->
                                        umiData[0].append(seq.charAt(from + pos))
                                        umiData[1].append(qual.charAt(from + pos))
                                    }
                                }

                                // trim adapter
                                if (trimBc) {
                                    int to = barcodeList[1][sampleIndex].size() + from
                                    readData[4 - readIndex] = (from > trimSizeThreshold ? seq.substring(0, from) : "") +
                                            ((seq.length() - to > trimSizeThreshold) ? seq.substring(to) : "")
                                    readData[5 - readIndex] = (from > trimSizeThreshold ? qual.substring(0, from) : "") +
                                            ((seq.length() - to > trimSizeThreshold) ? qual.substring(to) : "")

                                    if (removeTS) {
                                        seq = readData[4 - readIndex]
                                        qual = readData[5 - readIndex]
                                        def tsHit = (seq =~ gtPattern)
                                        if (tsHit) {
                                            readData[readIndex + 1] = seq.substring(0, tsHit.start())
                                            readData[readIndex + 2] = qual.substring(0, tsHit.start())
                                        }
                                    }
                                }
                            }
                        }
                        if (slaveFound) {
                            // 2.2 Report full match
                            writeQueue.put(wrapRead(readData, umiData, readIndex,
                                    sampleIds[sampleIndex], sampleIds[sampleIndex], true))
                        } else {
                            // 2.3 Report failed slave
                            writeQueue.put(wrapRead(readData, umiData, readIndex,
                                    'undef-s', sampleIds[sampleIndex], false))
                        }
                    } else { // Single-end, master match
                        writeQueue.put(wrapRead(readData, umiData, readIndex,
                                sampleIds[sampleIndex], sampleIds[sampleIndex], true))
                    }
                } else {
                    // 1.3 Report failed master
                    writeQueue.put(wrapRead(readData, umiData, readIndex,
                            'undef-m', 'undef-m', false))
                }
            }
        }
    } as Runnable))
}

def writeThread = new Thread({  // Writing thread
    String[] result
    while (true) {
        result = writeQueue.take()

        if (result.length == 0)
            break

        def sampleId = paired ? result[6] : result[3]

        def writerTrio = writers.get(sampleId)

        if (!paired)
            writerTrio[0].writeLine(result[0] + "\n" + result[1] + "\n+\n" + result[2])
        else {
            if (result[7]) {
                writerTrio[2].writeLine(result[0] + "\n" + result[1] + "\n+\n" + result[2])
            } else {
                writerTrio[0].writeLine(result[0] + "\n" + result[1] + "\n+\n" + result[2])
                writerTrio[1].writeLine(result[3] + "\n" + result[4] + "\n+\n" + result[5])
            }
        }
    }
} as Runnable)
writeThread.start()

readThread.join() // wait for read to finish

futures.each { it.get() } // wait for processing to finish
writeQueue.put(new String[0]) // tell writers this is last one
pool.shutdown()

writeThread.join()  // wait for write to finish

writers.values().each { // don't forget to flush
    it[0].close()
    if (paired) {
        it[1].close()
        if (overlap)
            it[2].close()
    }
}

new File("$out/checkout.filelist.txt").withPrintWriter { pw ->
    pw.println("\tSAMPLE_ID\tSAMPLE_TYPE\tCHECKOUT_FASTQ1\tCHECKOUT_FASTQ2")
    new HashSet(sampleIds).each { String sampleId -> // only unique
        pw.println(sampleId + (paired ? "\tpaired\t" : "\tunpaired\t") +
                new File("$out/${sampleId}_R1.fastq.gz").absolutePath + "\t" +
                (paired ? new File("$out/${sampleId}_R2.fastq.gz").absolutePath : "-"))
        if (overlap)
            pw.println(sampleId + "\toverlapped\t" +
                    new File("$out/${sampleId}_R12.fastq.gz").absolutePath) + "\t-"
    }
}

new File("$out/checkout.log.txt").withPrintWriter { pw ->
    pw.println("\t" + counters.collect { it.key }.join("\t"))
    pw.println("Master\t" + counters.collect { it.value[0] }.join("\t"))
    if (paired) {
        pw.println("Master+slave\t" + counters.collect { it.value[1] }.join("\t"))
        if (overlap)
            pw.println("Overlapped\t" + counters.collect { it.value[2] }.join("\t"))
    }
}
println "[${new Date()} $scriptName] Finished"