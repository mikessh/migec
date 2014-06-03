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

def scriptName = getClass().canonicalName

def cli = new CliBuilder(usage: "$scriptName [options] [checkout_dir/ or -] [assemble_dir/ or -] output_dir/")
cli.p(args: 1, 'number of threads to use')
cli._(longOpt: 'sample-info', args: 1, argName: 'file name',
        "A tab-delimited file indicating which samples to process and containing 6 columns:" +
                "sample_id (tab) species (tab) chain (tab) file_types (tab) mask (tab) quality_threshold_pair\n" +
                "Sample ids not included will be omitted by CdrBlast. If this argument is not provided, default values will be used.\n" +
                "Allowed species: ${Util.SPECIES.join(", ")}\nAllowed chains:${Util.CHAINS}\n" +
                "File types column is optional (it could be either omitted or set to \'-\'). " +
                "This columns should contain comma-separated list of file types to be processed for a given sample, e.g. \'paired, overlapped\'" +
                "Allowed file types are: ${Util.FILE_TYPES}. By default all file types are considered." +
                "Mask column is optional (it could be either omitted or set to \'-\') and will be ignored for unpaired and overlapped reads. " +
                "Mask which read(s) in pair will be searched for CDR3 region. " +
                "Allowed values are: 1:0 (R1 assembled), 0:1 (R2 assembled) and 1:1 (both reads assembled, [default])\n" +
                "Quality threshold column is optional (it could be either omitted or set to '-'). " +
                "Comma-separated pair of quality threshold values for Phred and CQS quality thresholds respectively" +
                "[default = 25,30]")
cli._(args: 1, longOpt: 'blast-path', 'Path to blast executable.')
cli._(longOpt: 'all-segments', 'Use full V/D/J segment library (including pseudogens, etc).')
cli._(longOpt: 'default-mask', args: 1, "Mask, default for all samples, see --sample-info")
cli._(longOpt: 'default-chain', args: 1, "Chain, default for all samples, see --sample-info")
cli._(longOpt: 'default-species', args: 1, "Species, default for all samples, see --sample-info")
cli._(longOpt: 'default-file-types', args: 1, "Accepted file types, default for all samples, see --sample-info")
cli._(longOpt: 'default-quality-threshold', args: 1, "Quality threshold pair, default for all samples, see --sample-info")

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 3) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(-1)
}

def blastPath = opt.'blast-path' ?: null

// INPUT AND OUTPUT FILES
String sampleInfoFileName = opt.'sample-info' ?: null, checkoutDir = opt.arguments()[0],
       assembleDir = opt.arguments()[1], outputPath = opt.arguments()[2]

if (!sampleInfoFileName && !opt.'default-chain') {
    println "[ERROR] Either --sample-info or --default-chain options should be provided"
    System.exit(-1)
}

boolean processRaw = checkoutDir != "-", processAssembled = assembleDir != "-"
if (!processRaw && !processAssembled) {
    println "[ERROR] At least one of assembled or raw (checkout) data directories should be specified"
    System.exit(-1)
}
boolean processBoth = processRaw && processAssembled

// LOAD SAMPLE ENTRIES
Map rawSamplesMap = null
if (processRaw) {
    def checkoutSampleListFile = new File(checkoutDir + "/checkout.filelist.txt")
    if (!checkoutSampleListFile.exists()) {
        println "[ERROR] Sample list (checkout.filelist.txt) absent in Checkout folder $checkoutDir"
        System.exit(-1)
    }

    rawSamplesMap = checkoutSampleListFile.readLines().findAll {
        !it.startsWith("#")   // no headers
    }.collectEntries {
        def splitLine = it.split('\t')
        def sampleName = splitLine[0], sampleType = splitLine[1]
        [("$sampleName\t$sampleType".toString()): splitLine[2..3]]
    }
}

Map assembledSampleMap = null
if (processAssembled) {
    def assembleSampleListFile = new File(assembleDir + "/assemble.log.txt")
    if (!assembleSampleListFile.exists()) {
        println "[ERROR] Assemble output (assemble.log.txt) absent in Assemble folder $assembleDir"
        System.exit(-1)
    }

    assembledSampleMap = assembleSampleListFile.readLines().findAll {
        !it.startsWith("#")   // no headers
    }.collectEntries {
        def splitLine = it.split('\t')
        def sampleName = splitLine[0], sampleType = splitLine[1], sampleKey = "$sampleName\t$sampleType".toString()

        [(sampleKey): [splitLine[2..3], splitLine[4..5]]]
    }
}

// GENERATE SAMPLE METADATA FOR CDRBLAST
def sampleInfoMap = new HashMap()
List<String> sampleInfoLines = []

if (sampleInfoFileName) {
    // Read sample info
    def chainsFile = new File(sampleInfoFileName)

    if (!chainsFile.exists()) {
        println "[ERROR] Sample info file not found"
        System.exit(-1)
    }

    sampleInfoLines = chainsFile.readLines()
} else {
    // Or generate it from default parameters
    def defaultMask = opt.'default-mask' ?: "1:1", defaultSpecies = opt.'default-species' ?: "human",
        defaultFileTypes = opt.'default-file-types' ?: Util.FILE_TYPES.join(","),
        defaultQualityThreshold = opt.'default-quality-threshold' ?: "25,30", defaultChain = opt.'default-chain'

    def defaultParameters = [defaultSpecies, defaultChain,
                             defaultFileTypes, defaultMask, defaultQualityThreshold].join("\t")

    (processRaw ? rawSamplesMap : assembledSampleMap).keySet().each {
        def sampleId = it.split("\t")[0]
        sampleInfoLines.add("$sampleId\t" + defaultParameters)
    }
}

// CHECK SAMPLE METADATA
sampleInfoLines.findAll { !it.startsWith("#") }.each { line ->
    def splitLine = line.split("\t")
    def sampleId = splitLine[0],
        species = splitLine[1].toLowerCase(), chain = splitLine[2].toUpperCase(),
        fileTypes = splitLine.length > 3 ? (splitLine[3] == '-' ? Util.FILE_TYPES.join(",") : splitLine[3]) : Util.FILE_TYPES.join(","),
        mask = (splitLine.length > 4 ? (splitLine[4] == '-' ? "1:1" : splitLine[4]) : "1:1"),
        qualityThreshold = splitLine.length > 5 ? (splitLine[5] == '-' ? null : splitLine[5]) : null

    qualityThreshold = qualityThreshold ? qualityThreshold.split(",") : null

    if (!Util.SPECIES.any { species == it }) {
        println "[ERROR] Bad species $species in line $line. Supported species are ${Util.SPECIES.join(", ")}"
        System.exit(-1)
    }
    if (!Util.CHAINS.any { chain == it }) {
        println "[ERROR] Bad chain $chain in line $line. Supported chains are ${Util.CHAINS.join(", ")}"
        System.exit(-1)
    }
    if (!Util.MASKS.any { mask == it }) {
        println "[ERROR] Bad mask $mask in line $line. Supported masks are ${Util.MASKS.join(", ")}"
        System.exit(-1)
    }
    mask.split(",").collect { it == "1" }
    fileTypes = fileTypes.split(",")
    fileTypes.each { fileType ->
        if (!Util.FILE_TYPES.any { fileType == it }) {
            println "[ERROR] Bad file type $fileType in line $line. Supported file types are ${Util.FILE_TYPES.join(", ")}"
            System.exit(-1)
        }
    }

    // check if sample exists
    boolean exists = false
    fileTypes.each {
        def sampleKey = "$sampleId\t$it".toString()
        if (processRaw && rawSamplesMap.containsKey(sampleKey))
            exists = true
        if (processAssembled && assembledSampleMap.containsKey(sampleKey))
            exists = true
    }
    if (!exists) {
        println "[ERROR] Sample $sampleId with file types ${fileTypes.join(",")} does not exist in checkout/assembled sample lists"
        System.exit(-1)
    }

    sampleInfoMap.put(sampleId, [species, chain, fileTypes, mask, qualityThreshold])
}

// START LOGGING
def logFile = new File(outputPath + "/cdrblast.log.txt")
if (logFile.exists())
    logFile.delete()
else
    logFile.absoluteFile.parentFile.mkdirs()

logFile.withPrintWriter { pw ->
    pw.println(Util.CDRBLAST_LOG_HEADER)

    // RUN CDRBLAST
    boolean allSegments = opt.'all-segments' ? true : false
    def baseArgs = [["-o"]]
    if (blastPath)
        baseArgs = [baseArgs, ["--blast-path", blastPath]]
    if (allSegments)
        baseArgs = [baseArgs, ["--all-segments"]]
    String stats
    if (processAssembled) {
        if (processBoth)
            println "[${new Date()} $scriptName] Running CdrBlast for raw and assembled data.."
        else
            println "[${new Date()} $scriptName] Running CdrBlast for assembled data.."

        sampleInfoMap.each { chainData ->
            def sampleId = chainData.key, species = chainData.value[0], chain = chainData.value[1],
                fileTypes = chainData.value[2], mask = chainData.value[3], qualityThreshold = chainData.value[4]

            def assemblyFiles = [], rawFiles = []

            fileTypes.each { fileType ->
                def sampleKey = "$sampleId\t$fileType".toString()

                if (assembledSampleMap.containsKey(sampleKey)) {
                    def correspondingFiles = assembledSampleMap[sampleKey]

                    // apply mask
                    if (fileType == 'paired') {
                        boolean zeroMask = true
                        if (mask[0]) {
                            rawFiles.add(correspondingFiles[0][0])
                            assemblyFiles.add(correspondingFiles[1][0])
                            if (correspondingFiles[1][0] != '-')
                                zeroMask = false
                        }
                        if (mask[1]) {
                            rawFiles.add(correspondingFiles[0][1])
                            assemblyFiles.add(correspondingFiles[1][1])
                            if (correspondingFiles[1][1] != '-')
                                zeroMask = false
                        }
                        if (zeroMask) {
                            println "[ERROR] Sample $sampleId with mask ${mask.collect { it ? 1 : 0 }} " +
                                    "does not have any corresponding assembled files, " +
                                    "of ${correspondingFiles[1]}"
                            System.exit(-1)
                        }
                        assemblyFiles.removeAll { it == '-' } // cleanup from files masked on assembly stage
                    } else {
                        rawFiles.add(correspondingFiles[0][0])
                        assemblyFiles.add(correspondingFiles[1][0])
                    }
                }
            }

            def baseArgs1 = [baseArgs, ["-R", chain], ["-S", species]]

            stats = Util.run(new CdrBlast(), [baseArgs1, ["-a"], ["-q", qualityThreshold[1]],
                                              assemblyFiles, "$outputPath/${sampleId}.asm.cdrblast.txt"].flatten().join(" "))
            pw.println(sampleId + "\tasm\t" + stats)

            if (processBoth) {
                stats = Util.run(new CdrBlast(), [baseArgs1, ["-q", qualityThreshold[0]],
                                                  rawFiles, "$outputPath/${sampleId}.raw.cdrblast.txt"].flatten().join(" "))
                pw.println(sampleId + "\traw\t" + stats)
            }
        }
    } else {
        println "[${new Date()} $scriptName] Running CdrBlast for raw data.."
        sampleInfoMap.each { chainData ->
            def sampleId = chainData.key, species = chainData.value[0], chain = chainData.value[1],
                fileTypes = chainData.value[2], mask = chainData.value[3], qualityThreshold = chainData.value[4]

            def rawFiles = []

            fileTypes.each { fileType ->
                def sampleKey = "$sampleId\t$fileType".toString()

                if (rawSamplesMap.containsKey(sampleKey)) {
                    def correspondingFiles = rawSamplesMap[sampleKey]

                    // apply mask
                    if (fileType == 'paired') {
                        if (mask[0])
                            rawFiles.add(correspondingFiles[0])

                        if (mask[1])
                            rawFiles.add(correspondingFiles[1])
                    } else
                        rawFiles.add(correspondingFiles[0])
                }
            }

            stats = Util.run(new CdrBlast(), [baseArgs, ["-R", chain], ["-S", species], ["-q", qualityThreshold[0]],
                                              rawFiles, "$outputPath/${sampleId}.raw.cdrblast.txt"].flatten().join(" "))
            pw.println(sampleId + "\traw\t" + stats)
        }
    }
}