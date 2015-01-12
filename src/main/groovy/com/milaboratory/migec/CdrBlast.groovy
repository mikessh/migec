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

import org.codehaus.groovy.runtime.ResourceGroovyMethods

//////////////
//   CLI   //
////////////
def cli = new CliBuilder(usage: "CdrBlast [options] reads1.fastq[.gz] [reads2.fastq[.gz] ...] output_file\n" +
        "NOTE: NCBI-BLAST+ package required")
cli.h("usage")
cli.R(args: 1, argName: "gene",
        "Receptor gene: 'TRA', 'TRB', 'TRG', 'TRD',  'IGL', 'IGK' or 'IGH' [required]. " +
                "Use --print-library to print the list of allowed species-gene combinations.")
cli.S(args: 1, argName: "species",
        "Species: 'HomoSapiens'[default], 'MusMusculus', ... " +
                "Use --print-library to print the list of allowed species-gene combinations.")
cli.a("Input data is assembled consensuses set, not reads. " +
        "The read header must contain 'UMI:NNNNN:COUNT' entry.")
cli.q(args: 1, argName: "Phred",
        "minimum quality of read to pass filter, will be " +
                "applied only to CDR3 region [default=30 for -a or 25, out of 40]")
cli.p(args: 1,
        "number of threads to use [default = all available processors]")
cli.N(args: 1,
        "Number of reads to take. For downsampling purposes")
cli._(longOpt: "no-sort",
        "Do not sort output by clonotype count, can save some time")
cli._(longOpt: "blast-path", args: 1, argName: "directory",
        "Path to blast executable.")
cli._(longOpt: "same-sample", "Assembled data only (-a). Input files come from the same sample so " +
        "each UMI will be counted once.")
cli._(longOpt: "all-segments",
        "Use full V/D/J segment library (including pseudogens, etc).")
cli._(longOpt: "print-library",
        "Prints out allowed species-gene pairs. " +
                "To account non-functional segment data use together with --all-segments")
cli._(longOpt: "strict-nc-handling",
        "Will not discard clonotypes with partial CDR3, " +
                "but rather try to append to corresponding complete clonotypes.")
cli._(longOpt: "debug",
        "Prints out alignment details, stores all temporary alignment files.")
cli._(longOpt: "cdr3-fastq-file", args: 1, argName: "file",
        "Store reads with extracted CDR3, appends CDR3 data in header. " +
                "Needed for 'com.milaboratory.migec.post.GroupByCdr' script.")
cli._(longOpt: "cdr3-umi-table", args: 1, argName: "file",
        "Stores UMIs and their associated clonotypes. Only applicable with -a option.")
cli._(longOpt: "log-file", args: 1, argName: "file",
        "File to output cdr extraction log.")
cli._(longOpt: "log-overwrite",
        "Overwrites provided log file.")
cli._(longOpt: "log-sample-name",
        "Sample name to use in log [default = N/A].")

def opt = cli.parse(args)

if (opt == null) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(-1)
}

if (opt.h) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(0)
}

// SEGMENTS STUFF
boolean includeNonFuncitonal = opt.'all-segments'
if (opt.'print-library') {
    println "CDR3 extraction is possible for the following data (segments include non-functional = $includeNonFuncitonal):"
    Util.listAvailableSegments(includeNonFuncitonal)
    System.exit(0)
}

if (!opt.R) {
    println "[ERROR] Receptor gene not provided"
    cli.usage()
    System.exit(-1)
}

String chain = opt.R, species = opt.S ?: "HomoSapiens"

if (!chain) {
    println "[ERROR] Chain argument is required for CdrBlast"
    System.exit(-1)
}

if (!Util.isAvailable(species, chain, includeNonFuncitonal)) {
    println "[ERROR] Sorry, no analysis could be performed for $species gene $chain " +
            "(include non-functional = $includeNonFuncitonal). " +
            "Possible variants are:\n"
    Util.listAvailableSegments(includeNonFuncitonal)
    System.exit(-1)
}

// SYSTEM
def DEBUG = opt.'debug',
    THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors(),
    SCRIPT_NAME = getClass().canonicalName

def timestamp = {
    "[${new Date()} $SCRIPT_NAME]"
}

def blastPath = opt.'blast-path' ?: ""
blastPath = blastPath.length() > 0 ? blastPath + "/" : ""

// Check for BLAST installation
try {
    ["${blastPath}convert2blastmask", "${blastPath}makeblastdb", "${blastPath}blastn"].each { it.execute().waitFor() }
} catch (IOException e) {
    println "[ERROR] Problems with BLAST installation. " + e.message
    System.exit(-1)
}

// INPUT, OUTPUT AND TEMPORARY FILES
if (opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(-1)
}

def inputFileNames = opt.arguments()[0..-2].collect { it.toString() },
    outputFileName = opt.arguments()[-1].toString()

def TMP_FOLDER = new File(inputFileNames[0]).absolutePath + "-cdrblast-" + UUID.randomUUID().toString()

def TMP_FOLDER_FILE = new File(TMP_FOLDER)
TMP_FOLDER_FILE.mkdirs()
//if (!DEBUG)
//    TMP_FOLDER_FILE.deleteOnExit()

def outputFile = null
if (outputFileName != "-" && (outputFile = new File(outputFileName)).parentFile)
    outputFile.parentFile.mkdirs()

// BLAST SETTINGS
int ALLELE_TAIL_INNER = 10, ALLELE_TAIL_OUTER = 6,
    ALLELE_TAIL_V_MAX = 40, ALLELE_TAIL_J_MAX = 20,
    GAP_OPEN = 5, GAP_EXTEND = 2, WORD_SIZE = 7, REWARD = 2, PENALTY = -3

String BLAST_FLAGS = "-lcase_masking"

// BLAST results filtering
int TOP_SEQS = 1,
    MIN_SEGMENT_IDENT = 7,
    MIN_CDR_LEN = 7, // at least 1 nt + conserved AAs
    MAX_CDR_LEN = 300 // there are some IGH with 35+ AAs in CDR3

// LOGGING
String logFileName = opt.'log-file' ?: null
boolean overwriteLog = opt.'log-overwrite'
String sampleName = opt.'log-sample-name' ?: "N/A"

// PARAMETERS
def cdr3FastqFile = opt.'cdr3-fastq-file', cdr3UmiTable = opt.'cdr3-umi-table'

boolean strictNcHandling = opt.'strict-nc-handling',
        assembledInput = opt.a,
        doSort = !opt.'no-sort',
        sameSample = opt.'same-sample'

int qualThreshold = opt.q ? Integer.parseInt(opt.q) : (opt.a ? 30 : 25),
    nReads = Integer.parseInt(opt.N ?: "-1")

////////////////////////////////////////////
// Stage 0: Pre-load segment information //
//////////////////////////////////////////
println "${timestamp()} Loading ${chain} Variable and Joining segment data"

class Allele {
    String alleleId, segmentId, seq, type
    int refPoint
}

class Segment {
    String id, type
    def alleles = new ArrayList<Allele>()
}

def segments = new HashMap<String, Segment>(), alleles = new HashMap<String, Allele>()
def vAlleles = new ArrayList<Allele>(), jAlleles = new ArrayList<Allele>()
def collapseAlleleMap = new HashMap<String, Allele>()

def resFile = Util.getSegmentsFile(includeNonFuncitonal)

resFile.splitEachLine("\t") {
    if (species.toUpperCase() == it[0].toUpperCase() &&
            chain.toUpperCase() == it[1].toUpperCase()) { // Take only alleles of a given chain and species
        def type = it[2].charAt(0).toUpperCase(), seq = it[5], refPoint = Integer.parseInt(it[4])

        // Mask all except seed region
        seq = seq.toLowerCase()
        seq = seq.toCharArray()
        int tailFrom = (type == "V" ? ALLELE_TAIL_OUTER : ALLELE_TAIL_INNER),
            tailTo = (type == "V" ? ALLELE_TAIL_INNER : ALLELE_TAIL_OUTER)
        for (int i = -tailFrom; i < tailTo; i++)
            if (refPoint + i >= 0 && refPoint + i < seq.length)
                seq[refPoint + i] = seq[refPoint + i].toUpperCase()
        seq = new String(seq)

        // Trim extending tails of alleles which are not covered by read
        if (type == "V") {
            if (ALLELE_TAIL_V_MAX <= refPoint) {
                seq = seq.substring(refPoint - ALLELE_TAIL_V_MAX)
                refPoint = ALLELE_TAIL_V_MAX
            }
        } else
            seq = seq.substring(0, Math.min(refPoint + ALLELE_TAIL_J_MAX + 1, seq.length()))
        def allele = collapseAlleleMap[seq] // Collapse alleles with identical sequences
        boolean newAllele = allele == null
        if (newAllele) {
            // it[0] - species
            allele = new Allele(alleleId: it[3], segmentId: it[3].split("[*]")[0],
                    refPoint: refPoint,
                    seq: seq,
                    type: type)
            collapseAlleleMap.put(seq, allele)
        } else
            allele.alleleId += "," + it[3]

        if (newAllele) {
            // V and J alleles separately - for alignment code simplicity
            if (allele.type == "V")
                vAlleles.add(allele)
            else
                jAlleles.add(allele)

            // Segment info - for allele frequencies
            def segment = segments[allele.segmentId]
            if (segment == null)
                segment = new Segment(id: allele.segmentId, type: allele.type)
            segment.alleles.add(allele)
            segments.put(segment.id, segment)
        }
    }
}

collapseAlleleMap.values().each {
    alleles.put(it.alleleId, it)
}

// Make blast db
["V", "J"].each { seg ->
    println "${timestamp()} Generating BLAST database for $chain $seg"
    new File("$TMP_FOLDER/${chain}_${seg}.fa").withPrintWriter { pw ->
        alleles.each { allele ->
            if (allele.value.type == seg)
                pw.println(">$allele.key\n$allele.value.seq")
        }
    }

    ("${blastPath}convert2blastmask -in $TMP_FOLDER/${chain}_${seg}.fa -out $TMP_FOLDER/${chain}_${seg}.msk " +
            "-masking_algorithm 'CDRBLAST' -masking_options 'NA'").execute().waitFor()

    ("${blastPath}makeblastdb -in $TMP_FOLDER/${chain}_${seg}.fa -mask_data $TMP_FOLDER/${chain}_${seg}.msk " +
            "-dbtype nucl -out $TMP_FOLDER/${chain}_${seg}").execute().waitFor()
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stage I: blast for unique sequences. Here we ONLY determine V and J segment presence and their IDs //
// seqSet (unique sequences) -> blast mapping data (unique sequence + V and J mapping)               //
//////////////////////////////////////////////////////////////////////////////////////////////////////
def seqMap = new HashMap<String, Integer>(), seqList = new ArrayList<String>()
int counter = 0, uniqueCounter = 0
inputFileNames.each { inputFileName ->
    println "${timestamp()} Reading $inputFileName and generating list of unique sequences to map V and J genes"
    def reader = Util.getReader(inputFileName)
    String header
    while (((header = reader.readLine()) != null) && (nReads < 0 || counter < nReads)) {
        if (assembledInput) {
            if (!header.contains("UMI:")) {
                println "[ERROR] Assembled input specified, but no UMI field in header"
                System.exit(-1)
            } else if (!header.split("[@ ]").find { it.startsWith("UMI") }.split(":")[2].isInteger()) {
                println "[ERROR] Assembled input specified, but UMI field in header does not contain count, " +
                        "maybe you're running on raw checkout output?"
                System.exit(-1)
            }
        }

        def seq = reader.readLine()
        if (!seqMap.containsKey(seq)) {
            seqMap.put(seq, uniqueCounter) // build a set of unique sequences to reduce the amount of data
            seqList.add(seq)
            uniqueCounter++
        }
        counter++
        reader.readLine()
        reader.readLine()

        if (counter % 1000000 == 0)
            println "${timestamp()} $counter reads processed, $uniqueCounter unique"
    }
}

def queryFilePrefix = TMP_FOLDER + "/" + Util.getFastqPrefix(inputFileNames[0])
println "[${new Date()} $SCRIPT_NAME] Making a temporary FASTA file for BLAST queries"
int chunk_sz = seqList.size() / THREADS
for (int p = 0; p < THREADS; p++) { // split fasta for blast parallelizaiton
    def file = new File(queryFilePrefix + "_${p}.fa")
    file.withPrintWriter { pw ->
        for (int i = p * chunk_sz; i < Math.min(seqList.size(), (p + 1) * chunk_sz); i++) {
            pw.println(">$i\n${seqList[i]}")
        }
    }
    //if (!DEBUG)
    //    file.deleteOnExit()
}

["V", "J"].each { seg -> // run blast for v and j segments separately
    println "${timestamp()} Pre-aligning $chain $seg segment with BLAST <$THREADS threads>"
    def blastProcs = (0..(THREADS - 1)).collect { p ->
        def blastOutFname = "${queryFilePrefix}_${chain}_${seg}_${p}.blast" // temp
        //if (!DEBUG)
        //    new File(blastOutFname).deleteOnExit()

        // A trick to pass -outfmt argument correctly
        def blastCmd = ["${blastPath}blastn",
                        "-query", "${queryFilePrefix}_${p}.fa",
                        BLAST_FLAGS.split(" "),
                        "-gapopen", "$GAP_OPEN", "-gapextend", "$GAP_EXTEND",
                        "-word_size", "$WORD_SIZE", "-reward", "$REWARD", "-penalty", "$PENALTY",
                        "-db", "$TMP_FOLDER/${chain}_${seg}",
                        "-max_target_seqs", "$TOP_SEQS", "-out", "$blastOutFname"]

        blastCmd = [blastCmd, "-outfmt", "6 qseqid sseqid qstart qend sstart send qseq sseq nident"].flatten()
        blastCmd.execute()
        //"bash $TMP_FOLDER/runblast.sh ${queryFilePrefix}_${p}.fa $TMP_FOLDER/${chain}_${seg} $blastOutFname".execute()
    }

    blastProcs.each { it.waitFor() }  // This is the only way to parallelize blast
    blastProcs.each { it.closeStreams() }
    blastProcs.each { it.destroy() }

    println "[${new Date()} $SCRIPT_NAME] Done"
}

//////////////////////////////////////////////////
// Stage II: Read blast results & extract CDRs //
////////////////////////////////////////////////
class BlastResult {
    int qFrom, qTo, aFrom, aTo
    boolean[] gaps
    Allele allele
    int nIdent

    BlastResult(Allele allele,
                int qFrom, int qTo,
                int aFrom, int aTo,
                String qseq,
                int nIdent) {
        this.allele = allele
        this.qFrom = qFrom
        this.qTo = qTo
        this.aFrom = aFrom
        this.aTo = aTo
        gaps = new boolean[qseq.length()]
        // we dont need any sequences further, only gaps in query
        for (int i = 0; i < qseq.length(); i++)
            gaps[i] = qseq.charAt(i) == '-'
        this.nIdent = nIdent
    }
}

def vMappingData = new HashMap<Integer, BlastResult>(),
    jMappingData = new HashMap<Integer, BlastResult>()

println "${timestamp()} Reading in BLAST results"
["V", "J"].eachWithIndex { segment, ind -> // parse blast output
    (0..(THREADS - 1)).each { p ->
        new File("${queryFilePrefix}_${chain}_${segment}_${p}.blast").splitEachLine("\t") { line ->
            int nIdent = Integer.parseInt(line[8])
            if (nIdent >= MIN_SEGMENT_IDENT) {
                def id = Integer.parseInt(line[0])
                def result = new BlastResult(
                        alleles[line[1]],
                        Integer.parseInt(line[2]) - 1, Integer.parseInt(line[3]) - 1,
                        Integer.parseInt(line[4]) - 1, Integer.parseInt(line[5]) - 1,
                        line[6],
                        nIdent
                )
                def mappingData = segment == 'V' ? vMappingData : jMappingData
                def prevMapping = mappingData[id]

                if (prevMapping == null || prevMapping.nIdent < nIdent)
                    mappingData.put(id, result)   // take result with highest identity
            }
        }
    }
}
println "${timestamp()} Done ${vMappingData.size()} V and ${jMappingData.size()} J mapped, " +
        "out of total ${seqMap.size()} sequence variants"

class ClonotypeData {
    Allele vAllele, jAllele
    String cdr3Seq
    int cdrFrom, cdrTo, vEndRel, jStartRel
    boolean rc
    int[] counts = new int[4] // counters: good events/total events/good reads/total reads with this V-J-CDR3
    String aaSeq

    ClonotypeData(Allele vAllele, Allele jAllele,
                  String cdr3Seq,
                  int cdrFrom, int cdrTo, int vEnd, int jStart,
                  boolean rc) {
        this.vAllele = vAllele
        this.jAllele = jAllele
        this.cdr3Seq = cdr3Seq
        this.cdrFrom = cdrFrom
        this.cdrTo = cdrTo
        this.vEndRel = vEnd - cdrFrom
        this.jStartRel = jStart - cdrFrom
        this.rc = rc
    }

    boolean isCanonical() {
        aaSeq = aaSeq ?: Util.translate(cdr3Seq)
        aaSeq.contains("?") || aaSeq.contains("*") || aaSeq.endsWith("F") || aaSeq.endsWith("W")
    }

    String getSignature() {
        [cdr3Seq, aaSeq = aaSeq ?: Util.translate(cdr3Seq),
         vAllele.alleleId, jAllele.alleleId, ".",
         vEndRel, ".", ".", jStartRel].join("\t") // for tabular output
    }

    String getShortSignature() {
        vAllele.alleleId + "|" + jAllele.alleleId + "|" + cdrFrom + "|" + cdr3Seq // for hash
    }

    @Override
    boolean equals(Object obj) {
        this.shortSignature == (obj as ClonotypeData).shortSignature
    }

    @Override
    int hashCode() {
        this.shortSignature.hashCode()
    }
}

def readId2ClonotypeData = new HashMap<Integer, ClonotypeData>()
def uniqueClonotypes = new HashMap<ClonotypeData, ClonotypeData>() // util
println "${timestamp()} Extracting CDRs"

int vRefNotInAlign = 0, jRefNotInAlign = 0, reverseFound = 0, badOrientation = 0,
    cdrStartOut = 0, cdrEndOut = 0, shortCdr = 0, longCdr = 0
int badExamples = 0, goodExamples = 0

seqMap.each {
    String seq = it.key
    int id = it.value
    def vMapping = vMappingData[id], jMapping = jMappingData[id]
    if (vMapping != null && jMapping != null) {
        Allele vAllele = vMapping.allele,
               jAllele = jMapping.allele
        boolean vRC = vMapping.aFrom > vMapping.aTo, jRC = jMapping.aFrom > jMapping.aTo

        // In case of rc flip so that start<end
        if (vRC) {
            int temp = vMapping.aTo
            vMapping.aTo = vMapping.aFrom
            vMapping.aFrom = temp
        }

        if (jRC) {
            int temp = jMapping.aTo
            jMapping.aTo = jMapping.aFrom
            jMapping.aFrom = temp
        }

        int vRef = vAllele.refPoint,
            jRef = jAllele.refPoint

        boolean vRefInAlign = vRef + 1 >= vMapping.aFrom && vRef - 1 <= vMapping.aTo,
                jRefInAlign = jRef + 1 >= jMapping.aFrom && jRef - 1 <= jMapping.aTo

        boolean bad = false

        int cdrFrom, cdrTo, cdrLen

        if (vRefInAlign && jRefInAlign && vRC == jRC) {
            if (vRC) { // after that step RC is same to non-RC
                reverseFound++
                seq = Util.revCompl(seq)

                int temp = vMapping.qTo
                vMapping.qTo = seq.length() - vMapping.qFrom - 1
                vMapping.qFrom = seq.length() - temp - 1

                temp = jMapping.qTo
                jMapping.qTo = seq.length() - jMapping.qFrom - 1
                jMapping.qFrom = seq.length() - temp - 1
            }

            // Determine v&j ref positions with respect to gaps
            vRef -= vMapping.aFrom
            int fixedRef = vRef
            for (int i = 0; i < fixedRef; i++)
                if (vMapping.gaps[i])
                    vRef--
            vRef += vMapping.qFrom

            jRef -= jMapping.aFrom
            fixedRef = jRef
            for (int i = 0; i < fixedRef; i++)
                if (jMapping.gaps[i])
                    jRef--
            jRef += jMapping.qFrom

            // Extract CDR3
            // Cys & Phe are included
            cdrFrom = vRef - 3
            cdrTo = jRef + 4
            cdrLen = cdrTo - cdrFrom

            if (cdrFrom >= 0 && cdrTo <= seq.length() && cdrLen >= MIN_CDR_LEN && cdrLen <= MAX_CDR_LEN) {
                // Do not create new CDR objects - save memory
                // CDR object has to store sequence string for hashing
                def clonotypeData = new ClonotypeData(vAllele, jAllele,
                        seq.substring(cdrFrom, cdrTo),
                        cdrFrom, cdrTo, vMapping.qTo, jMapping.qFrom, vRC) // to get hash

                if (uniqueClonotypes.containsKey(clonotypeData))
                    clonotypeData = uniqueClonotypes.get(clonotypeData) // already generated
                else
                    uniqueClonotypes.put(clonotypeData, clonotypeData) // make new

                // Unique read identifier -> CDR3 extraction result
                readId2ClonotypeData.put(id, clonotypeData)
            } else {
                if (cdrFrom < 0)
                    cdrStartOut++
                if (cdrTo > seq.length())
                    cdrEndOut++
                if (cdrLen < MIN_CDR_LEN)
                    shortCdr++
                if (cdrLen > MAX_CDR_LEN)
                    longCdr++
            }
        } else {
            if (!vRefInAlign)
                vRefNotInAlign++
            if (!jRefInAlign)
                jRefNotInAlign++
            if (vRC != jRC)
                badOrientation++
            bad = true
        }

        if (DEBUG) {
            if (bad) {
                if (badExamples < 5) {
                    println ">BAD"
                    println((vRC == jRC && vRC) ? Util.revCompl(seq) : seq)
                    if (!vRefInAlign)
                        println vAllele.alleleId + "\t" +
                                (vRC ? Util.revCompl(vAllele.seq.toUpperCase()) :
                                        vAllele.seq.toUpperCase()) + "\t" +
                                "aFrom=$vMapping.aFrom aTo=$vMapping.aTo ref=$vRef"

                    if (!jRefInAlign)
                        println jAllele.alleleId + "\t" +
                                (vRC ? Util.revCompl(jAllele.seq.toUpperCase()) :
                                        jAllele.seq.toUpperCase()) + "\t" +
                                "aFrom=$jMapping.aFrom aTo=$jMapping.aTo ref=$jRef"
                    badExamples++
                }
            } else if (goodExamples < 5) {
                println ">GOOD"
                println seq // already reversed
                println seq.substring(cdrFrom, cdrTo)
                println vAllele.alleleId + "\t" +
                        vAllele.seq.toUpperCase() + "\t" +
                        "qTo=$vMapping.qTo qFrom=$vMapping.qFrom aTo=$vMapping.aTo aFrom=$vMapping.aFrom ref=$vRef"
                println jAllele.alleleId + "\t" +
                        Util.revCompl(jAllele.seq.toUpperCase()) + "\t" +
                        "qTo=$jMapping.qTo qFrom=$jMapping.qFrom aTo=$jMapping.aTo aFrom=$jMapping.aFrom ref=$jRef"
                goodExamples++
            }
        }
    }
}
println "${timestamp()} Done. ${readId2ClonotypeData.size()} CDRs mapped, " +
        "out of total ${seqMap.size()} sequence variants"

if (DEBUG) {
    println "reverseFound=" + reverseFound + " (${(int) (reverseFound / seqMap.size() * 100)}%)"
    println "badOrientation=" + badOrientation
    println "vRefNotInAlign=" + vRefNotInAlign
    println "jRefNotInAlign=" + jRefNotInAlign
    println "cdrStartOut=" + cdrStartOut
    println "cdrEndOut=" + cdrEndOut
    println "shortCdr3=" + shortCdr + " (<$MIN_CDR_LEN nt)"
    println "longCdr3=" + longCdr + " (>$MAX_CDR_LEN nt)"
}

//////////////////////////////////
// Stage III: Append read data //
////////////////////////////////
counter = 0
int goodReads = 0, goodEvents = 0, mappedReads = 0, mappedEvents = 0, totalReads = 0, totalEvents = 0
def writerCdr3Fastq = cdr3FastqFile ? Util.getWriter(cdr3FastqFile) : null,
    writerCdr3Umi = assembledInput && cdr3UmiTable ? Util.getWriter(cdr3UmiTable) : null

if (writerCdr3Umi)
    writerCdr3Umi.println("#umi\tmig_sz\tv\tcdr3nt")

def usedUmis = new HashSet<String>()

inputFileNames.each { inputFileName ->
    println "${timestamp()} Appending read data from $inputFileName"
    def reader = Util.getReader(inputFileName)
    new File(TMP_FOLDER + "/bad_reads.txt").withPrintWriter { pw ->
        while (((header = reader.readLine()) != null) && (nReads < 0 || counter < nReads)) {
            def seq = reader.readLine()
            reader.readLine()
            def qual = reader.readLine(), oldQual = qual
            int id = seqMap[seq]  // get unique read sequence identifier
            def clonotypeData = readId2ClonotypeData[id] // corresponding CDR3 extraction result

            int increment = 1
            String umi
            boolean duplicateUmi = false

            if (assembledInput) {
                def splitHeader = header.split("[@ ]")
                def umiFields = splitHeader.find { it.startsWith("UMI") }.split(":")
                umi = umiFields[1]
                increment = Integer.parseInt(umiFields[2])
            }

            // Has CDR3
            if (clonotypeData != null) {
                if (clonotypeData.rc) // found in RC
                    qual = qual.reverse()

                qual = qual.substring(clonotypeData.cdrFrom, clonotypeData.cdrTo)

                // Increment counters if quality is good
                if (Util.minQual(qual) >= qualThreshold) {
                    if (assembledInput && sameSample) {
                        def umiCdr = umi + "_" + clonotypeData.cdr3Seq
                        duplicateUmi = usedUmis.contains(umi)
                        if (!duplicateUmi)
                            usedUmis.add(umiCdr)
                    }

                    if (!duplicateUmi) {
                        // Protect from duplicate counting when UMI overlap
                        clonotypeData.counts[0]++
                    }
                    goodEvents++

                    clonotypeData.counts[2] += increment
                    goodReads += increment

                    if (cdr3FastqFile)
                        writerCdr3Fastq.writeLine(header + " CDR3:" + clonotypeData.cdr3Seq + ":" +
                                clonotypeData.cdrFrom +
                                " V:" + clonotypeData.vAllele.alleleId +
                                " J:" + clonotypeData.jAllele.alleleId + "\n" +
                                seq + "\n+\n" + oldQual)

                    if (writerCdr3Umi)
                        writerCdr3Umi.println([umi, increment, clonotypeData.vAllele.alleleId,
                                               clonotypeData.cdr3Seq].join("\t"))
                }

                // Total (good+bad) for cases in which CDR3 was extracted
                if (!duplicateUmi) {
                    clonotypeData.counts[1]++
                }
                mappedEvents++
                clonotypeData.counts[3] += increment
                mappedReads += increment
            } else if (DEBUG)
                pw.println(header + "\n" + seq + "\n+\n" + qual)

            // Total, including cases when CDR3 was not extracted
            totalEvents++
            totalReads += increment
            counter++

            if (counter % 1000000 == 0)
                println "[${new Date()} $SCRIPT_NAME] $counter reads processed"
        }
    }
}
if (cdr3FastqFile)
    writerCdr3Fastq.close()
if (writerCdr3Umi)
    writerCdr3Umi.close()

println "${timestamp()} Done"
println "EVENTS (good mapped total):\t$goodEvents\t$mappedEvents\t$totalEvents\t|\t" +
        "${(int) (100 * (double) goodEvents / (double) totalEvents)}%\t" +
        "${(int) (100 * (double) mappedEvents / (double) totalEvents)}%\t100%"
println "READS (good mapped total):\t$goodReads\t$mappedReads\t$totalReads\t|\t" +
        "${(int) (100 * (double) goodReads / (double) totalReads)}%\t" +
        "${(int) (100 * (double) mappedReads / (double) totalReads)}%\t100%"

// Combine clonotypes based on signature, not shortSignature
// shortSignature used earlier needed for backtracking reads and is redundant because of cdrFrom
def cloneMap = new HashMap<String, int[]>()
def nonCanonicalCounts = new HashMap<String, int[]>()

if (strictNcHandling)
    uniqueClonotypes.keySet().each {
        if (!it.isCanonical()) {
            def counts = nonCanonicalCounts.get(it.cdr3Seq)
            if (counts == null)
                nonCanonicalCounts.put(it.cdr3Seq, counts = new int[4])
            for (int i = 0; i < 4; i++)
                counts[i] += it.counts[i]
        }
    }

uniqueClonotypes.keySet().each {
    if (it.isCanonical()) {
        def signature = it.signature
        def counts = cloneMap.get(signature)
        if (counts == null)
            cloneMap.put(signature, counts = new int[4])
        for (int i = 0; i < 4; i++)
            counts[i] += it.counts[i]

        if (strictNcHandling) {
            def ncToRemove = new HashSet<String>()
            nonCanonicalCounts.entrySet().each { ncEntry ->
                if (ncEntry.key.contains(it.cdr3Seq)) {
                    for (int i = 0; i < 4; i++)
                        counts[i] += ncEntry.value[i]
                    ncToRemove.add(ncEntry.key)
                }
            }
            ncToRemove.each { nc2remove ->
                nonCanonicalCounts.remove(nc2remove)
            }
        }
    }
}

if (outputFile) {
    println "${timestamp()} Writing output"
    outputFile.withPrintWriter { pw ->
        pw.println("Count\tPercentage\t" +
                "CDR3 nucleotide sequence\tCDR3 amino acid sequence\t" +
                "V segments\tJ segments\tD segments\t" +
                "Last V nucleotide position\t" +
                "First D nucleotide position\tLast D nucleotide position\t" +
                "First J nucleotide position\t" +
                "Good events\tTotal events\tGood reads\tTotal reads")

        (doSort ? cloneMap.sort { -it.value[0] } : cloneMap).each {
            if (it.value[0] > 0)
                pw.println(it.value[0] + "\t" + (it.value[0] / goodEvents) + // percentage from good events only
                        "\t" + it.key +
                        "\t" + it.value.collect().join("\t"))
        }
    }
}

// Those were created by blast and have to be removed manually
if (!DEBUG)
    ResourceGroovyMethods.deleteDir(TMP_FOLDER_FILE)
//TMP_FOLDER_FILE.listFiles().each { it.deleteOnExit() }

// Append to log and report to batch runner
def logLine = [(assembledInput ? "asm" : "raw"), outputFile ? outputFile.absolutePath : "-",
               inputFileNames.join(","),
               goodEvents, mappedEvents, totalEvents,
               goodReads, mappedReads, totalReads].join("\t")

if (logFileName) {
    def logFile = new File(logFileName)

    if (logFile.exists()) {
        if (overwriteLog)
            logFile.delete()
    } else {
        logFile.absoluteFile.parentFile.mkdirs()
        logFile.withPrintWriter { pw ->
            pw.println(Util.CDRBLAST_LOG_HEADER)
        }
    }

    logFile.withWriterAppend { logWriter ->
        logWriter.println("$sampleName\t" + logLine)
    }
}

return logLine