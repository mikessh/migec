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

//
// Internal routine checks how many different CDR3 sequences are there per a single MIG
// To be applied to the output of CdrBlast processing of raw data with --cdr3-fastq-file option
// CDR3s with a single mismatch are collapsed, as here we're interested only in spotting events when
// two different clonotypes got marked by the same UMI, not errors
//

package com.milaboratory.migec

def (inputFileName, outputFileName) = args[0..1]

def cdrCountsByUmi = new HashMap<String, HashMap<String, Integer>>()

println "[${new Date()}] Started reading"
new File(inputFileName).withReader { reader ->
    def header
    def counter = 0

    while ((header = reader.readLine()) != null) {
        def splitHeader = header.split(" "),
            umi = splitHeader.find { it.startsWith("UMI:") }.split(":")[1],
            cdr = splitHeader.find { it.startsWith("CDR3:") }.split(":")[1]

        def cdrCounts = cdrCountsByUmi[umi]
        if (!cdrCounts)
            cdrCountsByUmi.put(umi, cdrCounts = new HashMap<String, Integer>())
        cdrCounts.put(cdr, (cdrCounts[cdr] ?: 0) + 1)

        reader.readLine()
        reader.readLine()
        reader.readLine()

        if (++counter % 100000 == 0)
            println "[${new Date()}] $counter reads processed"
    }
}

println "[${new Date()}] Processing and writing output"
new File(outputFileName).withPrintWriter { pw ->
    pw.println("#umi\ttotal\tcdr3\tcount")

    int counter = 0

    cdrCountsByUmi.sort { -it.value.values().sum() }.each {
        def umi = it.key, cdrCounts = it.value
        int total = cdrCounts.values().sum()

        def cdrCountsCollapsed = new HashMap<String, Integer>()

        // collapse by 1mm
        cdrCounts.each {
            def cdr = it.key, cdrChars = cdr.toCharArray(),
                count = it.value
            boolean bad = false

            for (int i = 0; i < cdrChars.length; i++) {
                char oldChar = cdrChars[i]
                [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
                    if (nt != oldChar && !bad) {
                        cdrChars[i] = nt
                        def otherCdr = new String(cdrChars)
                        def otherCount = cdrCounts[otherCdr]

                        if (count < otherCount) {
                            cdrCountsCollapsed.put(otherCdr, (cdrCountsCollapsed[otherCdr] ?: 0) + count)
                            bad = true
                        }
                    }
                }
                cdrChars[i] = oldChar
                if (bad)
                    break
            }

            if (!bad)
                cdrCountsCollapsed.put(cdr, (cdrCountsCollapsed[cdr] ?: 0) + count)
        }

        cdrCountsCollapsed.sort { it.key }.sort { -it.value }.each {
            pw.println([umi, total, it.key, it.value].join("\t"))
        }

        if (++counter % 1000 == 0)
            println "[${new Date()}] $counter MIGs processed"
    }
}
