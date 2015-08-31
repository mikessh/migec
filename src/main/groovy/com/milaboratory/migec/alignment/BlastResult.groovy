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

import com.milaboratory.migec.segment.Allele

class BlastResult {
    int qFrom, qTo, aFrom, aTo
    boolean[] gaps
    Allele allele
    int nIdent

    BlastResult(Allele allele,
                int qFrom, int qTo,
                int aFrom, int aTo,
                String qseq,
                int nIdent) {
        this.allele = allele
        this.qFrom = qFrom
        this.qTo = qTo
        this.aFrom = aFrom
        this.aTo = aTo
        gaps = new boolean[qseq.length()]
        // we dont need any sequences further, only gaps in query
        for (int i = 0; i < qseq.length(); i++)
            gaps[i] = qseq.charAt(i) == '-'
        this.nIdent = nIdent
    }
}
