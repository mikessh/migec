.. _logs:

MIGEC log structure
-------------------

Below is the description of log files produced by various MIGEC routines.

**Checkout**

De-multiplexing, barcode extraction and overlapping:

* ``INPUT_FILE_1`` first input file containing R1 reads
* ``INPUT_FILE_2`` second input file containing R2 reads
* ``SAMPLE`` sample name
* ``MASTER`` number of reads where primary (master) barcode was detected
* ``SLAVE`` number of reads where secondary (slave) barcode was detected
* ``MASTER+SLAVE`` number of reads where both barcodes were
* ``OVERLAPPED`` number of succesfully overlapped reads

**Histogram**

The routine produces a number of histograms for UMI coverage, i.e. statistics of the number of reads tagged with a given UMI:

* ``overseq.txt`` contains sample id and sample type (single/paired/overlapped) in the header, followed by UMI coverage (MIG size). Each row has total read counts for UMIs corresponding to a given UMI coverage
* ``overseq-units.txt`` same as ``overseq.txt``, but lists numbers of unique UMIs, not total read counts
* ``estimates.txt`` contains sample id, sample type, total number of reads (``TOTAL_READS``) and UMIs (``TOTAL_MIGS``) in the sample and selected thresholds: ``OVERSEQ_THRESHOLD`` - UMI coverage threshold, ``COLLISION_THRESHOLD`` - if greater or equal to ``OVERSEQ_THRESHOLD`` will search for UMIs that differ by a single mismatch and have a huge count difference and treat them as being the same UMI, ``UMI_QUAL_THRESHOLD`` - threshold for min UMI sequence quality, ``UMI_LEN`` - UMI length
* ``collision1.txt`` - same as ``overseq.txt``, but lists only UMIs that are likely to be erroneous (i.e. have a 1-mismatch UMI neighbour with a substantially higher count)
* ``collision1-units.txt`` - same as ``collision1.txt``, but lists numbers of unique UMIs, not total read counts
* ``pwm.txt`` and ``pwm-units.txt`` - a position weight matrix (PWM) representation of all UMI sequences

**Assemble**

Statistics of MIG (group of reads tagged with the same UMI) consensus sequence assembly. Note that it also contains summary of pre-filtering steps, e.g. UMIs with low coverage are filtered at this stage:

* ``SAMPLE_ID`` sample name
* ``SAMPLE_TYPE`` sample type (single/paired/overlapped)
* ``INPUT_FASTQ1`` first input file containing R1 reads
* ``INPUT_FASTQ2`` second input file containing R2 reads
* ``OUTPUT_ASSEMBLY1`` first output file containing R1 consensuses
* ``OUTPUT_ASSEMBLY2`` second output file containing R2 consensuses
* ``MIG_COUNT_THRESHOLD`` UMI coverage threshold used in assemble procedure
* ``MIGS_GOOD_FASTQ1`` number of succesfully assembled consensuses from R1
* ``MIGS_GOOD_FASTQ2`` same for R2
* ``MIGS_GOOD_TOTAL`` number of succesfully assembled consensuses that have both R1 and R2 parts
* ``MIGS_TOTAL`` total number of input UMIs prior to coverage filtering
* ``READS_GOOD_FASTQ1`` number of reads in succesfully assembled consensuses from R1
* ``READS_GOOD_FASTQ2`` same for R2
* ``READS_GOOD_TOTAL`` number of paired reads in succesfully assembled consensuses that have both R1 and R2 parts. If a given assembled consensus contains inequal number of reads in R1 and R2, an average number is added to this statistic
* ``READS_TOTAL`` total number of input reads prior to coverage filtering
* ``READS_DROPPED_WITHIN_MIG_1`` number of reads dropped during consensus assembly as they had high number of mismatches to the consensus in R1
* ``READS_DROPPED_WITHIN_MIG_2`` same for R2
* ``MIGS_DROPPED_OVERSEQ_1`` number of UMIs dropped due to insufficient coverage in R1
* ``MIGS_DROPPED_OVERSEQ_2`` same for R2
* ``READS_DROPPED_OVERSEQ_1`` number of reads in UMIs dropped due to insufficient coverage in R1
* ``READS_DROPPED_OVERSEQ_2`` same for R2
* ``MIGS_DROPPED_COLLISION_1`` number of UMIs dropped due to being an erroneous (1-mismatch) variant of some UMI with higher count in R1
* ``MIGS_DROPPED_COLLISION_2`` same for R2
* ``READS_DROPPED_COLLISION_1`` number of reads in UMIs dropped due to being an erroneous (1-mismatch) variant of some UMI with higher count in R1
* ``READS_DROPPED_COLLISION_2`` same for R2

**CdrBlast**

Statistics of V(D)J mapping with BLAST algorithm:

* ``SAMPLE_ID`` sample name
* ``DATA_TYPE`` raw reads (raw) or assembled consensuses (asm)
* ``OUTPUT_FILE`` output file name
* ``INPUT_FILES`` list of input files
* ``EVENTS_GOOD`` number of MIGs (group of reads tagged with the same UMI, equals to number of reads for raw data) that were V(D)J mapped and passed the quality threshold
* ``EVENTS_MAPPED`` number of MIGs that were V(D)J mapped
* ``EVENTS_TOTAL`` number of input MIGs
* ``READS_GOOD`` number of reads that were V(D)J mapped and passed the quality threshold
* ``READS_MAPPED`` number of reads that were V(D)J mapped
* ``READS_TOTAL`` number of input reads


**FilterCdrBlastResults**

Statistics of the second round of TCR/Ig clonotype filtering that considers the number of supporting reads before and after consensus assembly:

* ``SAMPLE_ID`` sample name
* ``OUTPUT_FILE`` output file name
* ``INPUT_RAW`` input file containing CdrBlast results for raw reads
* ``INPUT_ASM`` input file containing CdrBlast results for assembled consensuses
* ``CLONOTYPES_FILTERED`` number of clonotypes (unique TCR/Ig V+CDR3 nucleotide+J combinations) that were filtered
* ``CLONOTYPES_TOTAL`` number of input clonotypes
* ``EVENTS_FILTERED`` number of MIGs in filtered clonotypes
* ``EVENTS_TOTAL`` number of input MIGs
* ``READS_FILTERED`` number of reads in filtered clonotypes
* ``READS_TOTAL`` number of input reads
* ``NON_FUNCTIONAL_CLONOTYPES`` number of non-functional clonotypes that contain stop codon/frameshift in CDR3
* ``NON_FUNCTIONAL_EVENTS`` number of MIGs in non-functional clonotypes
* ``NON_FUNCTIONAL_READS`` number of reads in non-functional clonotypes
