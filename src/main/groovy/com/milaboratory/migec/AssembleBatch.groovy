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

def DEFAULT_MODE = "0:1"
def cli = new CliBuilder(usage: 'AssembleBatch [options] checkout.filelist.txt histogram.estimates.txt output_directory')
cli.p(args: 1, 'number of threads to use')
cli._(longOpt: 'assembly-mode', args: 1, argName: '0:1, 1:0 or 1:1',
        "Mask for read(s) in pair that should be assembled be assembled. Will be used for pair-end samples only." +
                " Default: \"$DEFAULT_MODE\".")
cli.c("compressed output")

def scriptName = getClass().canonicalName
def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(-1)
}

def filelistFileName = opt.arguments()[0], estimatesFileName = opt.arguments()[1],
    outputPath = opt.arguments()[2], logFile = outputPath + '/assembly.log.txt'

// Some messy argument passing to Assemble
def assemblyModePaired = opt.'assembly-mode' ?: DEFAULT_MODE
def assemblyIndices = assemblyModePaired.split(':').collect { Integer.parseInt(it) > 0 }

if (!assemblyIndices.any()) {
    println "ERROR Bad assembly mode ${opt.'assembly-mode'}. At least one of reads should be indicated for assembly"
    System.exit(-1)
}

def baseArgs = []
if (opt.c)
    baseArgs.add(['-c'])
if (opt.p)
    baseArgs.add(['-p', opt.p])

// File names from Checkout output
def sampleFileNamesMap = new File(filelistFileName).readLines().findAll {
    !it.startsWith("#")   // no headers
}.collectEntries {
    def splitLine = it.split('\t')
    [(splitLine[0..1].join('\t')): splitLine[2..3]]
}

// Remove existing assemble log
if (new File(logFile).exists())
    new File(logFile).delete()

double collisionFactorThreshold = 0.05

println "[${new Date()} $scriptName] Starting batch assembly.."
new File(estimatesFileName).readLines().findAll { !it.startsWith("#") }.each {   // skip header
    def splitLine = it.split('\t')

    // Parse threshold estimates, check if it is safe to filter collisions
    def totalUmis = splitLine[3].toInteger(), overseqThreshold = splitLine[4].toInteger(),
        collThreshold = splitLine[5].toInteger(), umiQualThreshold = Byte.parseByte(splitLine[6]),
        umiLen = splitLine[7].toInteger(), filterCollisions = false

    if (collThreshold >= overseqThreshold &&
            totalUmis < collisionFactorThreshold * Math.pow(4, umiLen - 1)) // # collisions << # starting molecules
        filterCollisions = true

    // More messy argument passing
    def assembleArgs = [baseArgs, ['-m', collThreshold], ['-q', umiQualThreshold]]

    if (filterCollisions)
        assembleArgs.add(['--filter-collisions'])

    def sampleName = splitLine[0], sampleType = splitLine[1]
    if (sampleType == 'overlapped')
        assembleArgs.add(['--assembly-mode', '0:0'])
    else if (sampleType == 'paired')
        assembleArgs.add(['--assembly-mode', assemblyModePaired])

    // Pass filenames for I/O
    def sampleFileNames = sampleFileNamesMap["$sampleName\t$sampleType".toString()]
    assembleArgs.add(sampleFileNames)
    assembleArgs.add([outputPath + "/" + sampleName])
    assembleArgs.add([logFile])

    Util.run(new Assemble(), assembleArgs.flatten().join(" "))
}