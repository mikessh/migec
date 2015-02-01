package com.milaboratory.migec.alignment

import com.milaboratory.migec.segment.Allele

/**
 * Created by mikesh on 1/30/15.
 */
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
