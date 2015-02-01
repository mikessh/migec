package com.milaboratory.migec.dalign;

import com.milaboratory.migec.segment.Allele;

import java.util.List;

public class SearchResult {
    private final int segmentStart, segmentEnd;
    private final double score;
    private final List<Allele> alleles;

    public SearchResult(int segmentStart, int segmentEnd, double score, List<Allele> alleles) {
        this.segmentStart = segmentStart;
        this.segmentEnd = segmentEnd;
        this.score = score;
        this.alleles = alleles;
    }

    public int getSegmentStart() {
        return segmentStart;
    }

    public int getSegmentEnd() {
        return segmentEnd;
    }

    public double getScore() {
        return score;
    }

    public List<Allele> getAlleles() {
        return alleles;
    }
}
