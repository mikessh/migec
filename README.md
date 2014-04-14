=====================================================================
  MiGEC: Molecular Identifier Group-based Error Correction pipeline  
=====================================================================

This pipeline provides several useful tools for analysis of immune repertoire sequencing data. The pipeline utilizes unique nucleotide tags (UMIs) in order to filter experimental errors from resulting sequences. Those tags are attached to molecules before sequencing library preparation and allow to backtrack the original sequence of molecule. This pipeline is applicable for Illumina MiSeq and HiSeq 2500 reads. Sequencing libraries targeting CDR3 locus of immune receptor genes with high over-sequencing, i.e. ones that have at least 10 reads (optimally 30+ reads) per each starting molecule, should be used.

The data from 454 platform should be used with caution, as it contains homopolymer errors which (in present framework) result in reads dropped during consensus assembly. The 454 platform has a relatively low read yield, so additional read dropping could result in over-sequencing level below required threshold. If you still wish to give it a try, we would recommend filtering off all short reads and repairing indels with Coral (http://www.cs.helsinki.fi/u/lmsalmel/coral/), the latter should be run with options ```-mr 2 -mm 1000 -g 3```.

Features:
- Flexible de-multiplexing of NGS data and extraction of UMI sequence
- Assembly of consensuses of original molecules
- Extraction of CDR3 regions and determination of V/J genes for human and mouse immune receptors (TRA/TRB/TRG/TRD and IGH/IGL/IGK)
- Additional filtering of hot-spot errors

To run the steps of the pipeline either download the scripts and install Java + Groovy and run as (e.g. for the first step of pipeline, "Checkout") 

>$groovy Checkout.groovy

or simply download a standalone jar and execute

>$java -cp migec-v1.0.2.jar Checkout


STANDARD PIPELINE
=================

1. Checkout
==============================
Description: A script to perform de-multiplexing and UMI tag extraction

Standard usage: 
>$java -cp migec-v1.0.2.jar Checkout -cu barcodes.txt R1.fastq.gz R2.fastq.gz ./checkout/

For unpaired library:
>$java -cp migec-v1.0.2.jar Checkout -cu barcodes.txt R.fastq.gz - ./checkout/

barcodes.txt format is the following, 
>SAMPLE-ID (tab) MASTER-ADAPTER-SEQUENCE (tab) SLAVE-ADAPTER-SEQUENCE

A sequencing read is scanned for master adapter and then, if found, its mate is scanned for slave adapter. R2.fastq.gz and slave adapter sequence could be omitted.
Adaptor sequnce is accepted with any IUPAC DNA letters. Upper and lower case letters mark seed and fuzzy-search region parts respectively. 'N' characters mark UMI region to be extracted.
E.g. 
>S1 (tab) acgtacgtAAGGTTcadkgagNNNNNN

will search for AAGGTT seed exact match, then for the remaining adapter sequence with two mismatches allowed and output the NNNNNN region to header.

Additional parameters:

```-o``` could speed up if reads are oriented (i.e. master adapter should be in R1).

```-r``` will apply a custom RC mask. By default it assumes Illumina reads with mates on different strands, so it reverse-complements read with slave adapter so that output reads will be on master strand.

```--rc-barcodes``` also searches for both adapter sequences in reverse complement. Use it if unsure of your library structure.




2. Histogram
==============================
Description: A script to generate over-sequencing statistics

Standard usage:
>$java -cp migec-v1.0.2.jar Histogram ./checkout/checkout.filelist.txt ./checkout/histogram

Will generate several files, the one important for basic data processing is ./checkout/histogram.overseq.txt. The header contains MIG sizes, while each row for a sample is the number of reads in MIGs of a given size. Plot of MIG size (log coordinates) vs number of reads in MIG should display a clear peak. This plot should be used to set a MIG size cutoff in Assemble. 




3. Assemble
==============================
Description: A script to perform UMI-guided assembly

Standard usage:

>$java -cp migec-v1.0.2.jar Assemble -c ./checkout/S1_R1.fastq.gz ./checkout/S1_R2.fastq.gz ./assembly/S1 ./assembly/assembly.log

For unpaired library:

>$java -cp migec-v1.0.2.jar Assemble -c ./checkout/S1_R1.fastq - ./assembly/S1 ./assembly/assembly.log


All reads are grouped by their UMI and then read groups (aka molecular identifier groups, MIGs) with >10 reads (default value, see Histogram.groovy for details on setting it) are assembled. Multiple alignment is performed and consensus sequence is generated.

By default both only ./assembly/S1_R2.fastq.gz is assembled. To assemble both reads, use ```-m``` option, execute 
>$groovy Assemble -c -m 1:1 ./checkout/S1_R2.fastq.gz ./checkout/S1_R2.fastq.gz ./assembly/S1

In case of library with overlapping reads, the script can try to overlap them prior to assembly to generate a single output: 
>$groovy Assemble -c -m 0:0 ./checkout/S1_R2.fastq.gz ./checkout/S1_R2.fastq.gz ./assembly/S1

which will generate ./assembly/S1_RO.fastq.gz, containing assembly results _only_ for overlapping reads.

The ```--min-count``` option sets minimum number of reads in MIG.




4. CdrBlast
===============================
Description: A script to extract CDR3 sequences

Standard usage (assuming library contains T-cell Receptor Alpha Chain sequences)

For assembled data:

>$java -cp migec-v1.0.2.jar CdrBlast -a -C TRA ./assembly/S1_R2.fastq.gz ./cdr3blast/S1_asm.cdr3blast.txt 

For raw data:

>$java -cp migec-v1.0.2.jar CdrBlast -C TRA ./checkout/S1_R2.fastq.gz ./cdr3blast/S1_raw.cdr3blast.txt


NOTE:

1) NCBI-BLAST+ package required

2) Both raw and assembled data should be processed to apply the last step of filtration.




5. FilterCdrBlastResults
============================================
Description: A script to filter erroneous CDR3 sequences produced due to hot-spot PCR and NGS errors

Standard usage: 

>$java -cp migec-v1.0.2.jar FilterCdrBlastResults -s ./cdr3blast/S1_asm.cdr3blast.txt ./cdr3blast/S1_raw.cdr3blast.txt ./final/S1.cdr3blast.txt

The ```-s``` option tells to filter CDR3s represented by single MIGs, as for deep profiling (with our protocol) they could be associated with reverse transcription errors and experimental artifacts. Now the file S1.cdr3blast.txt contains a filtered CDR3/V/J clonotype table. Note that translated CDR3 sequences are obtained by simultaneously translating codons in two directions: from V and J segments to the middle of CDR3. If a frameshift is detected, the incomplete codon is added in lower case, with missing nucleotides marked as "?"; stop codons are marked by "*". CDR3 that contain either frameshift or stop codon are non-functional and could be filtered using ```-n``` option. You could additionally build a graph of hypermutations for the sample using

>$java -cp migec-v1.0.2.jar CreateCdrHypermGraph ./final/S1.cdr3blast.txt ./net

which will generate files that allow fast network construction using Cytoscape's network from table and import table routines.
