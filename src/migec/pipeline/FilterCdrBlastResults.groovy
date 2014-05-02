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

def R_A_T = "1.0"
def cli = new CliBuilder(usage:
        'groovy FilterCdrBlastResults [options] inputAssembledResult inputRawResult outputResult')
cli.r(args: 1, argName: 'read accumulation threshold', "Only clonotypes that have a ratio of (reads after correction) / " +
        "(uncorrected reads) greater than that threshold are retained. Default: $R_A_T")
cli.s("Include clonotypes that are represented by single events (have only one associated MIG)")
cli.n("Include non-functional CDR3s")
cli.c("Include CDR3s that do not begin with a conserved C or end with a conserved W/F")

def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 3) {
    cli.usage()
    return
}

def scriptName = getClass().canonicalName
def readAccumulationThreshold = Double.parseDouble(opt.r ?: R_A_T)
def filterUnits = !opt.s, filterNonFunctional = !opt.n, includeNonCanonical = opt.c
def inputUmiFile = opt.arguments()[0], inputRawFile = opt.arguments()[1], outputFile = opt.arguments()[2]

new File(outputFile).mkdirs()

int NT_SEQ_COL = 2, AA_SEQ_COL = 3,
    READ_COUNT_COL = 13, READ_TOTAL_COL = 14, EVENT_COUNT_COL = 11, EVENT_TOTAL_COL = 12,
    DATA_FROM = 2, DATA_TO = 14
def rawReadCounts = new HashMap<String, Integer>()
def totalRawReads = 0, totalConsensusReads = 0
int n = 0
println "[${new Date()} $scriptName] Reading raw clonotypes from $inputRawFile.."
new File(inputRawFile).splitEachLine("\t") { line ->
    if (line[0].isInteger()) {
        if (!filterNonFunctional || line[AA_SEQ_COL].matches("[a-zA-Z]+"))
            rawReadCounts.put(line[NT_SEQ_COL], (rawReadCounts[line[NT_SEQ_COL]] ?: 0) +
                    Integer.parseInt(line[READ_COUNT_COL]))
        totalRawReads += Integer.parseInt(line[READ_TOTAL_COL])

        if (++n % 500000 == 0)
            println "[${new Date()} $scriptName] $n clonotypes read"
    }
}

// Set V and J segments for a given CDR3nt as the ones with top count, collapse by CDR3
def cdr2signature = new HashMap<String, String>()
def cdr2count = new HashMap<String, int[]>()
println "[${new Date()} $scriptName] Reading assembled clonotypes from $inputUmiFile and filtering"
n = 0
new File(inputUmiFile).splitEachLine("\t") { line ->
    if (line[0].isInteger()) {
        String cdrSeq = line[NT_SEQ_COL]
        def signature = cdr2signature[cdrSeq]
        int goodEvents = Integer.parseInt(line[EVENT_COUNT_COL]), eventsTotal = Integer.parseInt(line[EVENT_TOTAL_COL]),
            goodReads = Integer.parseInt(line[READ_COUNT_COL]), readsTotal = Integer.parseInt(line[READ_TOTAL_COL])

        if ((!filterNonFunctional || line[AA_SEQ_COL].matches(/[a-zA-Z]+/))
                && (includeNonCanonical || line[AA_SEQ_COL].matches(/^C(.+)[FW]$/))) {
            // first col in signature is counter
            if (signature == null || Integer.parseInt(signature.split("\t")[0]) < goodReads)
                cdr2signature.put(cdrSeq, [goodReads, line[DATA_FROM..DATA_TO]].flatten().join("\t"))

            int[] counters = cdr2count[cdrSeq]
            if (counters == null)
                cdr2count.put(cdrSeq, counters = new int[4])
            counters[0] += goodEvents
            counters[1] += eventsTotal
            counters[2] += goodReads
            counters[3] += readsTotal
        }

        totalConsensusReads += readsTotal

        if (++n % 500000 == 0)
            println "[${new Date()} $scriptName] $n clonotypes read"
    }
}

println "[${new Date()} $scriptName] Filtering.."
int totalUnits = 0
// Filter: at least 1 good event & read accumulation > 100% (by default)
def passFilter = { counters, rawReads ->
    (!filterUnits || counters[0] > 1) &&
            (rawReads == null ||
                    (rawReads != null &&
                            counters[2] > readAccumulationThreshold *
                            rawReads * totalConsensusReads / totalRawReads))
}
new File(outputFile).withPrintWriter { pw ->
    pw.println("Count\tPercentage\t" +
            "CDR3 nucleotide sequence\tCDR3 amino acid sequence\t" +
            "V segments\tJ segments\tD segments\t" +
            "Last V nucleotide position\t" +
            "First D nucleotide position\tLast D nucleotide position\t" +
            "First J nucleotide position\t" +
            "Good events\tTotal events\tGood reads\tTotal reads")

    // 1st pass - compute total
    cdr2count.each {
        def counters = it.value, rawReads = rawReadCounts[it.key]
        if (passFilter(counters, rawReads))
            totalUnits += counters[0]
    }

    cdr2count.sort { -it.value[0] }.each {
        def signature = cdr2signature[it.key].split("\t")[1..-1].join("\t") // omit counter
        def counters = it.value, rawReads = rawReadCounts[it.key]
        if (passFilter(counters, rawReads))
            pw.println(counters[0] + "\t" + (counters[0] / totalUnits) + "\t" + signature)
    }
}

println "[${new Date()} $scriptName] Finished"
