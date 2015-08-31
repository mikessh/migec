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

package com.milaboratory.migec.alignment

import com.milaboratory.migec.Util

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
