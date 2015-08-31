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

package com.milaboratory.migec.dalign;

import com.milaboratory.migec.segment.Allele;

import java.util.LinkedList;
import java.util.List;

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
