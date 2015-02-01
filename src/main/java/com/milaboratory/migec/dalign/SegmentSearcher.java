package com.milaboratory.migec.dalign;

import com.milaboratory.migec.segment.Allele;

import java.util.LinkedList;
import java.util.List;

/**
 * Created by mikesh on 1/30/15.
 */
public class SegmentSearcher {
    private final int minHitSize;
    private final List<HitTracker> hitTrackers;

    public SegmentSearcher(List<Allele> alleles) {
        this(alleles, 2);
    }

    public SegmentSearcher(List<Allele> alleles, int minHitSize) {
        this.hitTrackers = new LinkedList<>();
        for (Allele allele : alleles) {
            hitTrackers.add(new HitTracker(allele, minHitSize));
        }
        this.minHitSize = minHitSize;
    }

    public SearchResult scan(String cdr3seq, int vEndRel, int jStartRel) {
        if (jStartRel - vEndRel - 1 < minHitSize)
            return null;

        final String subSeq = cdr3seq.substring(vEndRel + 1, jStartRel);
        final List<Allele> hitAlleles = new LinkedList<>();

        // iterate from largest window to smallest one
        for (int i = subSeq.length(); i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < subSeq.length() - i; j++) {
                String kmer = subSeq.substring(j, j + i);

                boolean hasHit = false;
                for (HitTracker hitTracker : hitTrackers) {
                    if (hitTracker.hasHit(kmer)) {
                        hasHit = true;
                        hitAlleles.add(hitTracker.getAllele());
                    }
                }

                if (hasHit) {
                    int from = vEndRel + 1 + j;
                    return new SearchResult(from,
                            from + i - 1,
                            i,
                            hitAlleles
                    );
                }
            }
        }

        return null;
    }
}
