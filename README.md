# MiGEC: Molecular Identifier Group-based Error Correction pipeline  

This pipeline provides several useful tools for analysis of immune repertoire sequencing data. The pipeline utilizes unique nucleotide tags (UMIs) in order to filter experimental errors from resulting sequences. Those tags are attached to molecules before sequencing library preparation and allow to backtrack the original sequence of molecule. This pipeline is applicable for Illumina MiSeq and HiSeq 2500 reads. Sequencing libraries targeting CDR3 locus of immune receptor genes with high over-sequencing, i.e. ones that have at least 10 reads (optimally 30+ reads) per each starting molecule, should be used.

### FEATURES

- Flexible de-multiplexing of NGS data and extraction of UMI sequence

- Assembly of consensuses of original molecules

- Extraction of CDR3 regions and determination of V/J genes for human and mouse immune receptors (TRA/TRB/TRG/TRD and IGH/IGL/IGK)

- Additional filtering of hot-spot errors

### INSTALLATION AND RUNNING

The pipeline is written in Groovy (a Java scripting language) and distributed as an executable JAR. To install it get the latest [JRE](http://www.oracle.com/technetwork/java/javase/downloads/index.html) and download the executable from [releases section](https://github.com/mikessh/migec/releases). Then to ran a specific script from the pipeline, say **Checkout**, execute

```
$java -cp migec.jar com.milaboratory.migec.Checkout [arguments]
```

alternatively you can download the repository and compile it from source using [Maven](http://maven.apache.org/) (requires Maven version 3.0)

```
$git clone https://github.com/mikessh/migec.git
$cd migec/
$mvn clean install
```

### NOTE

The data from 454 platform should be used with caution, as it contains homopolymer errors which (in present framework) result in reads dropped during consensus assembly. The 454 platform has a relatively low read yield, so additional read dropping could result in over-sequencing level below required threshold. If you still wish to give it a try, we would recommend filtering off all short reads and repairing indels with [Coral](http://www.cs.helsinki.fi/u/lmsalmel/coral/), the latter should be run with options ```-mr 2 -mm 1000 -g 3```.

### IMPORTANT

NCBI-BLAST+ package is required. Could be directly installed on Linux using a command like $sudo apt-get ncbi-blast+ or downloaded and installed directly from here: <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>


## STANDARD PIPELINE

### 1. Checkout

**Description**

A script to perform de-multiplexing and UMI tag extraction

**Usage**

For paired-end data:

```
$java -cp migec.jar com.milaboratory.migec.Checkout -cute barcodes.txt R1.fastq.gz R2.fastq.gz ./checkout/
```

For unpaired library:

```
$java -cp migec.jar com.milaboratory.migec.Checkout -cute barcodes.txt R.fastq.gz - ./checkout/
```

For overlapping paired reads:

```
$java -cp migec.jar com.milaboratory.migec.Checkout -cute --overlap barcodes.txt R1.fastq.gz R2.fastq.gz - ./checkout/
```

accepted *barcodes.txt* format is a tab-delimited table with the following structure: 

Sample ID | Master barcode sequence     | Slave barcode sequence
----------|-----------------------------|-----------------------
S0        | acgtacgtAGGTTAcadkgag       |
S1        | acgtacgtGGTTAAcadkgag       | ctgkGTTCaat
S2        | acgtacgtAAGGTTcadkgagNNNNNN |
S3        | acgtacgtTAAGGTcadkgagNNNNNN | NNNNNNctgkGTTCaat

A sequencing read is scanned for master adapter and then, if found, its mate is reverse-complemented to get on the same strand as master read and scanned for slave adapter.

* Slave adapter sequence could be omitted.

* Adaptor sequence could contain any IUPAC DNA letters.

* Upper and lower case letters mark seed and fuzzy-search region parts respectively.

* *N* characters mark UMI region to be extracted.

For example, in case *S2* **Checkout** will search for *AAGGTT* seed exact match, then for the remaining adapter sequence with two mismatches allowed and output the *NNNNNN* region to header. In case *S3* in addition the slave read is scanned for *GTTC* seed, fuzzy match to the rest of barcode is performed and *NNNNNN* region is extracted and concatenated with UMI region of master read.

**Parameters**

General:

```-c``` compressed output (gzip compression).

```-t``` trim adapter sequence from output.

```-e``` also remove trails of template-switching (poly-G) for the case when UMI-containing adapter is added using reverse-transcription (cDNA libraries).

```--overlap``` will try to overlap reads (paired-end data only), non-overlapping and overlapping reads will be placed to *_R1/_R2* and *_R12* FASTQ files respectively.

Barcode search:

```-o``` speed up by assuming that reads are oriented, i.e. master adapter should be in R1

```-r``` will apply a custom RC mask. By default it assumes Illumina reads with mates on different strands, so it reverse-complements read with slave adapter so that output reads will be on master strand.

```--rc-barcodes``` also searches for both adapter sequences in reverse complement. Use it if unsure of your library structure.


### 2. Histogram

**Description**

A script to generate over-sequencing statistics

**Usage**

```
$java -cp migec.jar com.milaboratory.migec.Histogram ./checkout/checkout.filelist.txt ./histogram/run
```

Will generate several files, the one important for basic data processing is *./checkout/histogram.overseq.txt*. The header contains MIG sizes (in log2 scale), while each row for a sample contains the number of reads in MIGs of a given size (cumulative abundance). 

For a decent dataset the plot of cumulative abundance display a small peak at MIG size of 1 that could be attributed to erroneous MIGs and has an exponential decline, and a clear peak at MIG size of 10+ containing amplified MIGs. Those erroneous MIGs could arise as experimental artifacts, however the most common reason for their presence is an error event in UMI sequence itself. Note that the latter is only valid when number of distinct UMIs is far lower than theoretically possible UMI diversity (e.g. 4^12 for 12-letter UMI regions)!
 
MIG size cutoff in **Assemble** should be set to dissect erroneous MIGs while retaining amplified ones. If peaks overlap collision filtering should be considered.


### 3. Assemble

**Description**

A script to perform UMI-guided assembly

**Usage**

```
$java -cp migec.jar com.milaboratory.migec.Assemble -c ./checkout/S1_R1.fastq.gz ./assembly/S1 ./assembly/assembly.log
```

All reads are grouped by their UMI and then read groups (aka molecular identifier groups, MIGs) with >10 reads (default value, see Histogram.groovy for details on setting it) are assembled. Multiple alignment is performed and consensus sequence is generated.

**Settings**

The ```-m``` option sets minimum number of reads in MIG. This should be set according to Histogram script output to separate two peaks: over-sequenced MIGs and erroneous MIGs that cluster around MIG size of 1.

To inspect the effect of such single-mismatch erroneous UMI sub-variants see "collisions" output of Histogram script. Such collision events could interfere with real MIGs when over-sequencing is relatively low. In this case collisions could be filtered during MIG consensus assembly using ```--filter-collisions``` option.


### 4. CdrBlast

**Description** 

A script to extract CDR3 sequences

**Usage**

Standard, assuming library contains T-cell Receptor Alpha Chain sequences

in case of MIG-assembled data:

```
$java -cp migec.jar com.milaboratory.migec.CdrBlast -a -C TRA ./assembly/S1_R2.fastq.gz ./cdrblast/S1_asm.cdrblast.txt 
```
      
for raw data:

```
$java -cp migec.jar com.milaboratory.migec.CdrBlast -C TRA ./checkout/S1_R2.fastq.gz ./cdrblast/S1_raw.cdrblast.txt
```

To get a sorted output use ```-o``` option, otherwise sorting will be performed at **FilterCdrBlastResults** step. Note that both raw and assembled data should be processed to apply the last step of filtration.


### 5. FilterCdrBlastResults

**Description**

A script to filter erroneous CDR3 sequences produced due to hot-spot PCR and NGS errors

**Usage** 

```
$java -cp migec.jar com.milaboratory.migec.FilterCdrBlastResults -s ./cdrblast/S1_asm.cdrblast.txt ./cdrblast/S1_raw.cdrblast.txt ./final/S1.cdrblast.txt
```

The ```-s``` option tells to include CDR3s represented by single MIGs. Those are filtered by default as for deep profiling (with our protocol) they could be associated with reverse transcription errors and experimental artifacts.

Now the file *S1.cdrblast.txt* contains a filtered and sorted CDR3/V/J clonotype table.

You could additionally build a graph of hypermutations for the sample using

```
$java -cp migec.jar com.milaboratory.migec.CreateCdrHypermGraph ./final/S1.cdr3blast.txt ./net
```

which will generate files that allow fast network construction using Cytoscape's network from table and import table routines for further data exploration.

Note that translated CDR3 sequences are obtained by simultaneously translating codons in two directions: from V and J segments to the middle of CDR3. If a frameshift is detected, the incomplete codon is added in lower case, with missing nucleotides marked as ```?```; stop codons are marked by ```*```. CDR3 that contain either frameshift or stop codon are non-functional and are filtered by default. To include them into your output use ```-n``` option.