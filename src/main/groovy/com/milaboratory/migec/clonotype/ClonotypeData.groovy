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

package com.milaboratory.migec.clonotype

import com.milaboratory.migec.Util
import com.milaboratory.migec.dalign.SegmentSearcher
import com.milaboratory.migec.segment.Allele

import static com.milaboratory.migec.Util.BLANK_FIELD

class ClonotypeData {
    Allele vAllele, jAllele
    Allele[] dAlleles
    String cdr3Seq
    int cdrFrom, cdrTo,
        vEndRel, jStartRel,
        dStartRel = -1, dEndRel = -1
    boolean rc
    int[] counts = new int[4] // counters: good events/total events/good reads/total reads with this V-J-CDR3
    String aaSeq

    ClonotypeData(Allele vAllele, Allele jAllele,
                  String cdr3Seq,
                  int cdrFrom, int cdrTo, int vEnd, int jStart,
                  boolean rc) {
        this.vAllele = vAllele
        this.jAllele = jAllele
        this.cdr3Seq = cdr3Seq
        this.cdrFrom = cdrFrom
        this.cdrTo = cdrTo
        this.vEndRel = vEnd - cdrFrom
        this.jStartRel = jStart - cdrFrom
        this.rc = rc
    }

    boolean isCanonical() {
        aaSeq = aaSeq ?: Util.translate(cdr3Seq)
        aaSeq.contains("?") || aaSeq.contains("*") || aaSeq.endsWith("F") || aaSeq.endsWith("W")
    }

    String getdAllele() {
        (dAlleles && dAlleles.length > 0) ? dAlleles.collect { it.alleleId }.join(",") : BLANK_FIELD
    }

    String getSignature() {
        [cdr3Seq, aaSeq = aaSeq ?: Util.translate(cdr3Seq),
         vAllele.alleleId, jAllele.alleleId, dAllele,
         vEndRel, dStartRel, dEndRel, jStartRel].join("\t") // for tabular output
    }

    String getShortSignature() {
        vAllele.alleleId + "|" + jAllele.alleleId + "|" + cdrFrom + "|" + cdr3Seq // for hash
    }

    void appendDInfo(SegmentSearcher dSegmentSearcher) {
        def result = dSegmentSearcher.scan(cdr3Seq, vEndRel, jStartRel)
        if (result) {
            dStartRel = result.segmentStart
            dEndRel = result.segmentEnd
            dAlleles = result.alleles as Allele[]
        }
    }

    @Override
    boolean equals(Object obj) {
        this.shortSignature == (obj as ClonotypeData).shortSignature
    }

    @Override
    int hashCode() {
        this.shortSignature.hashCode()
    }
}
