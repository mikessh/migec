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

package migec.control

import java.util.zip.GZIPInputStream

if (args.length < 5) {
    println "Usage: HistQ inputFastqFile output from(e.g. 0) to(e.g. read length, 100)"
    System.exit(0)
}

def infiles = args[0].split(","), ofile = args[1]
int from = Integer.parseInt(args[2]), to = Integer.parseInt(args[3])

double[][] qhist = new double[42][to - from]
int n = 0

infiles.each { infile ->
    def reader = new BufferedReader(new InputStreamReader(infile.endsWith(".gz") ?
            new GZIPInputStream(new FileInputStream(infile)) : new FileInputStream(infile)))

    while (reader.readLine() != null) {
        reader.readLine()
        reader.readLine()
        def q = reader.readLine().toCharArray()
        for (int i = from; i < to; i++) {
            qhist[((int) q[i] - 33)][i - from]++
        }
        if (n++ % 1000000 == 0)
            println "[${new Date()}] Processed $n reads"
    }
}

if (new File(ofile).parentFile)
    new File(ofile).parentFile.mkdirs()

new File(ofile).withPrintWriter { pw ->
    pw.println(qhist.collect { it.collect().join("\t") }.join("\n"))
}
