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

def cli = new CliBuilder(usage: 'SplitSpikeClonotypes file_with_spikes file_with_cdrs output_prefix')
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 3) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}
def spikeFileName = opt.arguments()[0], inputFileName = opt.arguments()[1], outputFilePrefix = opt.arguments()[2]

def spikeList = new File(spikeFileName).readLines()
def spikeHash = new HashSet<String>(spikeList)
// add 1 mm variants also
spikeList.each {
    def chars = it.toCharArray()
    def oldChar
    for (int i = 0; i < chars.length; i++) {
        oldChar = chars[i]
        [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
            if (nt != oldChar) {
                chars[i] = nt
                spikeHash.add(new String(chars))
            }
        }
        chars[i] = oldChar
    }
}

// Iterate through file
int NT_SEQ_COL = 2
new File(outputFilePrefix + ".sample.txt").withPrintWriter { pw1 ->
    new File(outputFilePrefix + ".spikes.txt").withPrintWriter { pw2 ->
        new File(inputFileName).splitEachLine("\t") { line ->
            if (line[0].isInteger()) {
                String seq = line[NT_SEQ_COL]
                def chars = seq.toCharArray()
                def oldChar
                boolean isSpike = false
                for (int i = 0; i < chars.length; i++) {
                    oldChar = chars[i]
                    // two mismatches total
                    [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
                        chars[i] = nt
                        if (spikeHash.contains(new String(chars)))
                            isSpike = true
                    }
                    chars[i] = oldChar
                    if (isSpike)
                        break
                }
                if (isSpike)
                    pw2.println(line.join("\t"))
                else
                    pw1.println(line.join("\t"))
            } else {
                pw1.println(line.join("\t"))
                pw2.println(line.join("\t"))
            }
        }
    }
}