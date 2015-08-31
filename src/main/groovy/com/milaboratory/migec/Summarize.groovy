/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.milaboratory.migec

import java.nio.file.FileSystems
import java.nio.file.Files

def homeDir = "."

new File("$homeDir/summary/").deleteDir()
new File("$homeDir/summary/").mkdir()

boolean firstSample = true
new File(homeDir).eachDir { sampleDir ->
    def sampleName = sampleDir.name

    println "Processing $sampleName"

    sampleDir.eachDir { outputDir ->
        def outputName = outputDir.name

        // copy logs

        def logFile = outputDir.listFiles().find { it.name.endsWith(".log.txt") }

        if (logFile != null) {
            new File("$homeDir/summary/${outputName}.log.txt").withWriterAppend { w ->
                def logLines = logFile.readLines()

                if (firstSample)
                    w.println(logLines[0])

                logLines[1..-1].each { w.println(it) }
            }
        }

        // copy final clonotype output

        if (outputName == "cdrfinal")
            Files.copy(
                    FileSystems.default.getPath(outputDir.listFiles().find {
                        it.name.endsWith(".filtered.cdrblast.txt")
                    }.absolutePath),
                    FileSystems.default.getPath("$homeDir/summary/${sampleName}.filtered.cdrblast.txt")
            )

        // histograms - just stack all of them
        if (outputName == "histogram") {
            outputDir.listFiles().each { histFile ->
                new File("$homeDir/summary/$histFile.name").withWriterAppend { w ->
                    def histLines = histFile.readLines()

                    if (firstSample)
                        w.println(histLines[0])

                    histLines[1..-1].each { w.println(it) }
                }
            }
        }
    }

    firstSample = false
}