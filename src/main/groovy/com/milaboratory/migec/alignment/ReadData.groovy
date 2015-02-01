package com.milaboratory.migec.alignment

import com.milaboratory.migec.Util

/**
 * Created by mikesh on 1/30/15.
 */
class ReadData {
    int count = 0
    final long[] qualArr

    ReadData(String seq) {
        this.qualArr = new long[seq.length()]
    }

    void append(String qual, boolean group) {
        count++

        if (group)
            for (int i = 0; i < qual.length(); i++)
                qualArr[i] += Util.qualFromSymbol(qual.charAt(i))
    }

    String finalizeQual() {
        char[] qual = new char[qualArr.length]

        for (int i = 0; i < qualArr.length; i++)
            qual[i] = Util.symbolFromQual((byte) (qualArr[i] / count))

        new String(qual)
    }
}
