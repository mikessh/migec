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

def R_A_T = "1.0"
def scriptName = getClass().canonicalName
def cli = new CliBuilder(usage: "$scriptName [options] cdrblast_dir/ output_dir/")
cli.r(args: 1, argName: 'read accumulation threshold', "Only clonotypes that have a ratio of (reads after correction) / " +
        "(uncorrected reads) greater than that threshold are retained. Default: $R_A_T")
cli._(longOpt: 'collapse', "Collapse clonotypes by CDR3 and use top V and J chains")
cli.p(args: 1, "number of threads to use. Default: all available processors")
cli.s("Include clonotypes that are represented by single events (have only one associated MIG)")
cli.n("Include non-functional CDR3s")
cli.c("Include CDR3s that do not begin with a conserved C or end with a conserved W/F")

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(-1)
}

def baseArgs = [["-r", opt.r ?: R_A_T]]
if (opt.s)
    baseArgs = [baseArgs, ["-s"]]
if (opt.n)
    baseArgs = [baseArgs, ["-n"]]
if (opt.c)
    baseArgs = [baseArgs, ["-c"]]
if (opt.p)
    baseArgs = [baseArgs, ["-p", opt.p.toString()]]
if (opt.'collapse')
    baseArgs = [baseArgs, ["--collapse"]]

String inputDir = opt.arguments()[0], outputDir = opt.arguments()[1]

def cdrBlastBatchFile = new File("$inputDir/cdrblast.log.txt")
if (!cdrBlastBatchFile.exists()) {
    println "[ERROR] CdrBlast log file (cdrblast.log.txt) not found in input dir"
    System.exit(-1)
}

new File(outputDir).mkdirs()

// collect raw-assembled pairs
def rawSampleMap = new HashMap(), assembledSampleMap = new HashMap()

cdrBlastBatchFile.splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        if (splitLine[1].toLowerCase() == "asm") {
            assembledSampleMap.put(splitLine[0], splitLine[2])
        } else if (splitLine[1].toLowerCase() == "raw") {
            rawSampleMap.put(splitLine[0], splitLine[2])
        } else {
            println "[ERROR] Unknown sample type ${splitLine[1]}. Either \'raw\' or \'asm\' expected"
            System.exit(-1)
        }

        if (!new File(splitLine[2]).exists()) {
            println "[ERROR] File not found for sample ${splitLine[0..2]}."
            System.exit(-1)
        }
    }
}

def onlyInRaw = rawSampleMap.keySet().findAll { !assembledSampleMap[it] },
    onlyInAssembled = assembledSampleMap.keySet().findAll { !rawSampleMap[it] }

if (onlyInRaw.size() > 0) {
    println "[WARNING] The following samples have only raw data analyzed by CdrBlast and will be skipped: " +
            "${onlyInRaw.join(", ")}"
}
if (onlyInAssembled.size() > 0) {
    println "[WARNING] The following samples have only assembled data analyzed by CdrBlast and will be skipped: " +
            "${onlyInAssembled.join(", ")}"
}

if (onlyInRaw.size() == rawSampleMap.size()) {
    println "[ERROR] Np samples have both raw and assembled data analyzed by CdrBlast. Terminating"
    System.exit(-1)
}

println "[${new Date()} $scriptName] Running batch hot-spot error filtering for CdrBlast results.."
// START LOGGING
def logFile = new File(outputDir + "/cdrblastfilter.log.txt")
if (logFile.exists())
    logFile.delete()
else
    logFile.absoluteFile.parentFile.mkdirs()

logFile.withPrintWriter { pw ->

    pw.println(Util.CDRBLASTFILTER_LOG_HEADER)

    assembledSampleMap.each {
        def sampleId = it.key, asmFileName = it.value, rawFileName = rawSampleMap[sampleId]
        if (rawFileName) {
            def finalArgs = [baseArgs, [asmFileName, rawFileName, "$outputDir/${sampleId}.filtered.cdrblast.txt"]]
            def logLine = Util.run(new FilterCdrBlastResults(), finalArgs.flatten().join(" "))
            pw.println("$sampleId\t" + logLine)
        }
    }
}
