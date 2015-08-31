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
