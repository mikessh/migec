package com.milaboratory.migec.dalign;

import com.milaboratory.migec.segment.Allele;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by mikesh on 1/30/15.
 */
public class HitTracker {
    private final Allele allele;
    private final Set<String> kmers = new HashSet<>();

    public HitTracker(Allele allele, int minHitSize) {
        this.allele = allele;
        for (int i = minHitSize; i < allele.getSeq().length(); i++) {
            for (int j = 0; j < allele.getSeq().length() - i; j++) {
                String kmer = allele.getSeq().substring(j, j + i);
                kmers.add(kmer);
            }
        }
    }

    public boolean hasHit(String queryKmer) {
        return kmers.contains(queryKmer);
    }

    public Allele getAllele() {
        return allele;
    }
}
