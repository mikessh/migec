package com.milaboratory.migec

import java.util.regex.Pattern

class ReadOverlapper {
    final int maxOverlapOffset, overlapSeedSize, overlapFuzzySize, maxConsMms
    final boolean allowPartialOverlap
    final double maxOverlapMismatchRatio

    ReadOverlapper(int maxOverlapOffset, int overlapSeedSize, int overlapFuzzySize, 
                   boolean allowPartialOverlap = false, int maxConsMms = 2, double maxOverlapMismatchRatio = 0.1) {
        this.maxOverlapOffset = maxOverlapOffset
        this.overlapSeedSize = overlapSeedSize
        this.overlapFuzzySize = overlapFuzzySize
        this.maxConsMms = maxConsMms
        this.allowPartialOverlap = allowPartialOverlap
        this.maxOverlapMismatchRatio = maxOverlapMismatchRatio
    }

    String[] overlap(String r1, String r2, String q1, String q2) {
        String[] result = null

        for (int i = 0; i < maxOverlapOffset; i++) {
            if (i + overlapSeedSize > r2.length())
                return null

            def kmer = r2.substring(i, i + overlapSeedSize)
            def pattern = Pattern.compile(kmer)
            def matcher = pattern.matcher(r1)
            // Find last match
            int position
            while (matcher.find()) {
                position = matcher.start()
                if (position >= 0) {
                    // Start fuzzy align
                    boolean alignedAll = true
                    int nConsMms = 0, nMms = 0, actualFuzzyOverlapSize = overlapFuzzySize

                    for (int j = 0; j < overlapFuzzySize; j++) {
                        def posInR1 = position + overlapSeedSize + j, posInR2 = i + overlapSeedSize + j
                        if (posInR1 + 1 > r1.length() || posInR2 + 1 > r2.length()) {
                            actualFuzzyOverlapSize = j + 1
                            alignedAll = false
                            break     // went to end of r1
                        }
                        if (r1.charAt(posInR1) != r2.charAt(posInR2)) {
                            nMms++
                            if (++nConsMms >= maxConsMms)
                                break  // several consequent mismatches
                        } else {
                            nConsMms = 0 // zero counter
                        }
                    }

                    if (nConsMms < maxConsMms &&
                            (allowPartialOverlap || alignedAll) &&
                            (nMms / (double) actualFuzzyOverlapSize) <= maxOverlapMismatchRatio) {
                        // take best qual nts
                        def seq = new StringBuilder(r1.substring(0, position)),
                            qual = new StringBuilder(q1.substring(0, position))

                        int pos2 = i - 1
                        for (int j = position; j < r1.length(); j++) {
                            pos2++
                            if (pos2 == r2.length())
                                break // should not happen

                            seq.append(Util.qualFromSymbol(q1.charAt(j)) > Util.qualFromSymbol(q2.charAt(pos2)) ?
                                    r1.charAt(j) : r2.charAt(pos2))

                            qual.append(Util.qualFromSymbol(q1.charAt(j)) > Util.qualFromSymbol(q2.charAt(pos2)) ?
                                    q1.charAt(j) : q2.charAt(pos2))
                        }
                        for (int j = pos2 + 1; j < r2.length(); j++) {
                            // fill the remainder
                            seq.append(r2.charAt(j))

                            qual.append(q2.charAt(j))
                        }

                        // report overlap
                        result = new String[2]
                        result[0] = seq.toString()
                        result[1] = qual.toString()

                        return result
                    }
                }
            }
        }

        result // failed
    }
}
