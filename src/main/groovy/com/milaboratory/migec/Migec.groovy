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

import java.util.jar.JarFile

def version = (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
        /Implementation-Version: (.+)/)[0][1]

def printHelp = {
    println " MiGEC Pipeline V$version "
    println ""
    println "[Main pipeline]"
    println "Run as \$java -jar migec-${version}.jar SCRIPT_NAME arguments"
    println ""
    println "where SCRIPT_NAME is one of the following:"
    println "Checkout"
    println "CheckoutBatch"
    println "Histogram"
    println "Assemble"
    println "AssembleBatch"
    println "CdrBlast"
    println "CdrBlastBatch"
    println "FilterCdrBlastResults"
    println "FilterCdrBlastResultsBatch"
    println "CreateCdrHypermGraph"
    println ""
    println ""
    println "[Miscellaneous]"
    println "Executed via classpath"
    println "Run as \$java -cp migec-${version}.jar SCRIPT_NAME arguments"
    println ""
    println "where SCRIPT_NAME is one of the following:"
    println "com.milaboratory.migec.GroupByCdr"
    println "com.milaboratory.migec.SplitSpikeClonotypes"
    println "com.milaboratory.migec.BacktrackSequence"
    println "com.milaboratory.migec.HistQ"
    println "com.milaboratory.migec.UmiFrequencyTable"
}

def getScript = { String scriptName ->
    switch (scriptName.toUpperCase()) {
        case "CHECKOUT":
            return new Checkout()
        case "CHECKOUTBATCH":
            return new CheckoutBatch()
        case "HISTOGRAM":
            return new Histogram()
        case "ASSEMBLE":
            return new Assemble()
        case "ASSEMBLEBATCH":
            return new AssembleBatch()
        case "CDRBLAST":
            return new CdrBlast()
        case "CDRBLASTBATCH":
            return new CdrBlastBatch()
        case "FILTERCDRBLASTRESULTS":
            return new FilterCdrBlastResults()
        case "FILTERCDRBLASTRESULTSBATCH":
            return new FilterCdrBlastResultsBatch()
        case "CREATECDRHYPERMGRAPH":
            return new CreateCdrHypermGraph()
        case "-H":
        case "H":
        case "-HELP":
        case "HELP":
            printHelp()
            println ""
            System.exit(-1)
            break

        default:
            printHelp()
            println ""
            println "Unknown MAIN PIPELINE script $scriptName"
            System.exit(-1)
    }
}

if (args.length == 0)
    printHelp()
else {
    def script = getScript(args[0])
    try {
        Util.CMD_LINE = "migec-${version}.jar ${args.join(" ")}"
        Util.run(script, args.length > 1 ? args[1..-1].join(" ") : "")
    } catch (Exception e) {
        println "[ERROR] $e.message, see _migec_error.log for details"
        new File("_migec_error.log").withWriterAppend { writer ->
            writer.println("[${new Date()}]")
            writer.println("[Script]")
            writer.println(args[0])
            writer.println("[CommandLine]")
            writer.println("executing $Util.CMD_LINE")
            writer.println("[Message]")
            writer.println(e.message)
            writer.println("[StackTrace-Short]")
            writer.println(e.stackTrace.findAll { it.toString().contains("com.milaboratory.migec") }.join("\n"))
            writer.println("[StackTrace-Full]")
            e.printStackTrace(new PrintWriter(writer))
        }
        System.exit(-1)
    }
}