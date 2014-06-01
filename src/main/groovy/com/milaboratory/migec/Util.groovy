package com.milaboratory.migec

import java.util.concurrent.atomic.AtomicIntegerArray
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

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
class Util {
    static final byte DEFAULT_UMI_QUAL_THRESHOLD = (byte) 15
    static final char[] NTS = ['A', 'T', 'G', 'C']

    static byte qualFromSymbol(char symbol) {
        (int) symbol - 33
    }

    static byte[] qualFromString(String qual) {
        byte[] qualArr = new byte[qual.length()]
        for (int i = 0; i < qual.length(); i++)
            qualArr[i] = qualFromSymbol(qual.charAt(i))
        qualArr
    }

    static char symbolFromQual(int qual) {
        qual = qual < 2 ? 2 : qual
        qual = qual > 40 ? 40 : qual
        (char) (qual + 33)
    }

    static char code2nt(int code) {
        switch (code) {
            case 0:
                return 'A'
            case 1:
                return 'T'
            case 2:
                return 'G'
            case 3:
                return 'C'
        }
    }

    static byte nt2code(char symbol) {
        switch (symbol) {
            case 'A':
                return 0
            case 'T':
                return 1
            case 'G':
                return 2
            case 'C':
                return 3
        }
    }

    static BufferedReader getReader(String fname) {
        new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fname)) :
                new FileInputStream(fname)))
    }

    static BufferedWriter getWriter(String outfile, boolean compressed) {
        if (compressed)
            outfile += ".gz"
        new BufferedWriter(new OutputStreamWriter(compressed ?
                new GZIPOutputStream(new FileOutputStream(outfile)) : new FileOutputStream(outfile)))
    }

    static BufferedWriter getWriter(String outfile) {
        boolean compressed = outfile.endsWith(".gz")
        new BufferedWriter(new OutputStreamWriter(compressed ?
                new GZIPOutputStream(new FileOutputStream(outfile)) : new FileOutputStream(outfile)))
    }


    static String getUmi(String header, byte umiQualThreshold) {
        def splitHeader = header.split(" ")
        def umiEntry = splitHeader.find { it.startsWith("UMI:") }
        if (umiEntry == null) {
            println "[ERROR] no UMI header in input. Terminating"
            System.exit(-1)
        }
        String umi = umiEntry.split(":")[1] // quality can contain :
        for (int i = umiEntry.length() - umi.length(); i < umiEntry.length(); i++)
            if (qualFromSymbol(umiEntry.charAt(i)) < umiQualThreshold)
                return null
        umi
    }

    static String toString(AtomicIntegerArray arr) {
        (0..<arr.length()).collect { arr.get(it) }.join("\t")
    }

    static void run(Script script, String args) {
        script.binding.setVariable("args", args.split(" "))
        script.run()
    }
}
