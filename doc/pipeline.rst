The pipeline
------------

All routines in the pipeline are available in **manual** and **batch**
variants. Batch variants are designed to automatically handle several
input samples with minimal shell scripting glue between analysis steps, 
this is also the recommended way to use with MIGEC. If the "barcodes" file 
is set properly, the entire pipeline can be written as following:

.. code:: bash

    MIGEC="java -Xmx8G -jar MIGEC-$VERSION.jar"
    $MIGEC CheckoutBatch -cu barcodes.txt checkout/
    $MIGEC Histogram checkout/ histogram/
    $MIGEC AssembleBatch -c checkout/ histogram/ assemble/
    $MIGEC CdrBlastBatch -R TRB checkout/ assemble/ cdrblast/
    $MIGEC FilterCdrBlastResultsBatch cdrblast/ cdrfinal/
    
Manual usage example
~~~~~~~~~~~~~~~~~~~~

An example for a 300bp paired-end MiSeq run of IGH library on a 16Gb RAM
Unix server. Such sequencing read length allows complete IGH sequencing,
thus mate pairs overlap. First *barcodes.txt* should be created
containing adapter sequences, see the section below for guidelines.
Then, assuming that the corresponding FASTQ files are
*IGH\_SAMPLE\_R1.fastq.gz* and *IGH\_SAMPLE\_R2.fastq.gz*, UMI- and
multiplex index-containing adapter is near 5'UTR of V segment (so the
CDR3 is in mate#2 after reads are oriented) and NCBI-BLAST+ is
installed, run all 5 stages of the pipeline using the following command:

.. code:: bash

    $MIGEC Checkout -cute --overlap barcodes.txt IGH_S1-10_R1.fastq.gz IGH_S1-10_R2.fastq.gz checkout/
    $MIGEC Histogram checkout/ histogram/
    $MIGEC Assemble -c --mask 0:0 checkout/S1_R12.fastq.gz . assembly/
    $MIGEC CdrBlast -R IGH checkout/S1_R12.fastq.gz cdrblast/S1_raw.txt
    $MIGEC CdrBlast -a -R IGH assembly/S1_R12.fastq.gz cdrblast/S1_asm.txt
    $MIGEC FilterCdrBlastResults cdrblast/S1_asm.txt cdrblast/S1_raw.txt cdrfinal/S1.txt
    
De-multiplexing
~~~~~~~~~~~~~~~

Checkout-batch
^^^^^^^^^^^^^^

**Description**

A script to perform de-multiplexing and UMI tag extraction for a set of
FASTQ files that were previously split using Illumina sample indices.

**Usage**

General:

.. code:: bash

    java -jar migec.jar CheckoutBatch [options] barcodes_file output_dir

The barcodes file specifies sample multiplexing and UMI (NNN.. region)
extraction rules. It has the same structure as for "manual" Checkout
(see section below), with additional two columns that specify input
FASTQ file names.

+-------------+-------------------------------+--------------------------+---------------------------+---------------------------+
| Sample ID   | Master barcode sequence       | Slave barcode sequence   | Read#1 FASTQ              | Read#2 FASTQ              |
+=============+===============================+==========================+===========================+===========================+
| S0          | acgtacgtAGGTTAcadkgag         |                          |                           |                           |
+-------------+-------------------------------+--------------------------+---------------------------+---------------------------+
| S1          | acgtacgtGGTTAAcadkgag         | ctgkGTTCaat              | ILM1\_R1\_L001.fastq.gz   | ILM1\_R2\_L001.fastq.gz   |
+-------------+-------------------------------+--------------------------+---------------------------+---------------------------+
| S1          | acgtacgtAAGGTTcadkgagNNNNNN   |                          | ILM2\_R1\_L001.fastq.gz   | ILM2\_R2\_L001.fastq.gz   |
+-------------+-------------------------------+--------------------------+---------------------------+---------------------------+
| S3          | acgtacgtTAAGGTcadkgagNNNNNN   | NNNNNNctgkGTTCaat        | ILM1\_R1\_L001.fastq.gz   | ILM1\_R2\_L001.fastq.gz   |
+-------------+-------------------------------+--------------------------+---------------------------+---------------------------+

The following rules apply:

-  All specified FASTQ files are sequentially processed using Checkout
-  If no FASTQ file is specified for a given barcode, it will be
   searched in all FASTQ files
-  CheckoutBatch will properly aggregate reads from multiple FASTQ files
   that have the same sample id
-  Still there should not be the case when a FASTQ file has the same
   barcode specified more than once

**Parameters**

Same as in manual version of Checkout, see below.

Checkout-manual
^^^^^^^^^^^^^^^

**Description**

A script to perform de-multiplexing and UMI tag extraction

**Usage**

General:

.. code:: bash

    java -jar migec.jar Checkout [options] barcodes_file R1.fastq[.gz] [. or R2.fastq[.gz]] output_dir

For paired-end data:

.. code:: bash

    java -jar migec.jar Checkout -cute barcodes.txt R1.fastq.gz R2.fastq.gz ./checkout/

For unpaired library:

.. code:: bash

    java -jar migec.jar Checkout -cute barcodes.txt R.fastq.gz . ./checkout/

For overlapping paired reads:

.. code:: bash

    java -jar migec.jar Checkout -cute --overlap barcodes.txt R1.fastq.gz R2.fastq.gz . checkout/

accepted *barcodes.txt* format is a tab-delimited table with the
following structure:

+-------------+-------------------------------+--------------------------+
| Sample ID   | Master barcode sequence       | Slave barcode sequence   |
+=============+===============================+==========================+
| S0          | acgtacgtAGGTTAcadkgag         |                          |
+-------------+-------------------------------+--------------------------+
| S1          | acgtacgtGGTTAAcadkgag         | ctgkGTTCaat              |
+-------------+-------------------------------+--------------------------+
| S2          | acgtacgtAAGGTTcadkgagNNNNNN   |                          |
+-------------+-------------------------------+--------------------------+
| S3          | acgtacgtTAAGGTcadkgagNNNNNN   | NNNNNNctgkGTTCaat        |
+-------------+-------------------------------+--------------------------+

A sequencing read is scanned for master adapter and then, if found, its
mate is reverse-complemented to get on the same strand as master read
and scanned for slave adapter.

-  Slave adapter sequence could be omitted.

-  Adaptor sequence could contain any IUPAC DNA letters.

-  Upper and lower case letters mark seed and fuzzy-search region parts
   respectively.

-  *N* characters mark UMI region to be extracted.

-  Multiple rows could correspond to the same sample

-  In order to be able to run batch pipeline operations, all samples
   should contain UMI region of the same size

For example, in case *S2* **Checkout** will search for *AAGGTT* seed
exact match, then for the remaining adapter sequence with two mismatches
allowed and output the *NNNNNN* region to header. In case *S3* in
addition the slave read is scanned for *GTTC* seed, fuzzy match to the
rest of barcode is performed and *NNNNNN* region is extracted and
concatenated with UMI region of master read.

**Parameters**

General:

``-c`` compressed output (gzip compression).

``-u`` perform UMI region extraction and output it to the header of
de-multiplexed FASTQ files

``-t`` trim adapter sequence from output.

``-e`` also remove trails of template-switching (poly-G) for the case
when UMI-containing adapter is added using reverse-transcription (cDNA
libraries).

``--overlap`` will try to overlap reads (paired-end data only),
non-overlapping and overlapping reads will be placed to \*\_R1/\_R2\*
and \*\_R12\* FASTQ files respectively. While overlapping the nucleotide
with higher quality will be taken thus improving overall data quality.

``--overlap-max-offset X`` controls to which extent overlapping region
is searched. **IMPORTANT** If the read-through extent is high (reads are
embedded) should be set to ~40.

Barcode search:

``-o`` speed up by assuming that reads are oriented, i.e. master adapter
should be in R1

``-r`` will apply a custom RC mask. By default it assumes Illumina reads
with mates on different strands, so it reverse-complements read with
slave adapter so that output reads will be on master strand.

``--rc-barcodes`` also searches for both adapter sequences in reverse
complement. Use it if unsure of your library structure.

``--skip-undef`` will not store reads that miss adapter sequence to save
drive space. **NOTE** When there is a huge number of unassigned/unused
reads this option greatly speeds up de-multiplexing. However, take care
to carefully investigate the reasons behind low barcode extraction rate
if it is a case.

MIG statistics
~~~~~~~~~~~~~~

**Description**

A script to generate consensus coverage statistics, i.e. molecular 
identifier group (MIG) size distribution.

**Usage**

General:

.. code:: bash

    java -jar migec.jar Histogram checkout/ histogram/

Running this script will generate several files in *histogram* folder,
the one important for basic data processing is *overseq.txt*. The header
of table contains MIG sizes (in log2 scale), while each row corresponds
to a de-multiplexed sample contains the number of reads in MIGs of a
given size (cumulative abundance).

For a decent dataset the plot of cumulative abundance display a small
peak at MIG size of 1 that could be attributed to erroneous MIGs and has
an exponential decline, and a clear peak at MIG size of 10+ containing
amplified MIGs. Those erroneous MIGs could arise as experimental
artifacts, however the most common reason for their presence is an error
event in UMI sequence itself. Note that the latter is only valid when
number of distinct UMIs is far lower than theoretically possible UMI
diversity (e.g. 4^12 for 12-letter UMI regions)!

MIG size cutoff in **Assemble** should be set to dissect erroneous MIGs
while retaining amplified ones. If peaks overlap collision filtering
should be considered.

A simple plotting routine written in R can facilitate visualization of
MIG size distributions, available
`here <https://github.com/mikessh/migec/tree/master/util>`__.

Consensus assembly
~~~~~~~~~~~~~~~~~~

Assemble-match
^^^^^^^^^^^^^^

**Description**

A script to perform UMI-guided assembly

**Usage**

General:

.. code:: bash

    java -jar migec.jar AssembleBatch [options] checkout_output_folder/ histogram_output_folder/ output_folder/

Performs a batch assembly for all FASTQ files produced by checkout, all
assembly parameters are set according to **Histogram** output.

One can specify a default mask telling for paired-end reads which
mate(s) to assemble. The mask is provided by
``--default-mask <R1=[0,1]:R2=[0,1]>`` argument, i.e. to assemble only
second mate use ``--default-mask 0:1``. This speeds-up the assembly.
Also, by default the mask is ``1:1``, so for each MIG an output
consensus pair is created only if both consensuses are successfully
assembled. In case of ``0:0`` mask will process only overlapped reads.
Remember that during **Checkout** reads get re-oriented so they are on
the same strand, corresponding to the strand of *Master* barcode and the
read with *Master* barcode is assigned with \*\_R1\* index.

A sample metadata file could also be provided with
``--sample-metadata <file_name>`` argument to guide the batch assembly.
This file should have the following tab-separated table structure:

+-------------+--------------+--------+
| Sample ID   | File type    | Mask   |
+=============+==============+========+
| S0          | paired       | 1:0    |
+-------------+--------------+--------+
| S0          | overlapped   |        |
+-------------+--------------+--------+
| S1          | unpaired     |        |
+-------------+--------------+--------+
| S2          | paired       | 0:1    |
+-------------+--------------+--------+

Note that *S0* is present with two file types, as when performing read
overlap **Checkout** stores non-overlapped reads in \*\_R1/\_R2\* files,
which could be then incorporated into data processing.

The ``--force-overseq X`` and ``--force-collision-filter`` will force a
MIG size threshold of ``X`` and filtering of 1-mm UMI collisions for all
samples being processed.

**IMPORTANT** In most cases, the automatic MIG size threshold selected
by Histogram routine is ok. However we strongly recommend manual
inspection of Histogram output files and considering to manually specify
an appropriate MIG size threshold for input samples.

Assemble-manual
^^^^^^^^^^^^^^^

**Description**

A script to perform UMI-guided assembly

**Usage**

General:

.. code:: bash

    java -jar migec.jar Assemble [options] R1.fastq[.gz] [. or R2.fastq[.gz]] output_folder

Unpaired and overlapped FASTQ:

.. code:: bash

    java -jar migec.jar Assemble -c checkout/S1_R0.fastq.gz . assembly/

Paired FASTQ:

.. code:: bash

    java -jar migec.jar Assemble -c checkout/S1_R1.fastq.gz checkout/S1_R2.fastq.gz ./assembly/

Paired FASTQ with only second read to be assembled:

.. code:: bash

    java -jar migec.jar Assemble -c --mask 0:1 checkout/S1_R1.fastq.gz checkout/S1_R2.fastq.gz assembly/

All reads are grouped by their UMI and then read groups (aka molecular
identifier groups, MIGs) with >10 reads (default value, see
**Histogram** section for details on setting it) are assembled. Multiple
alignment is performed and consensus sequence is generated. Note that
for paired reads both consensuses should be successfully assembled,
otherwise the pair is dropped.

Automatic output file naming convention is used for compatibility with
batch operations. Output file name will be appended with \_R0 for
unpaired FASTQ file, with either \_R1 and \_R2 for the corresponding
paired FASTQ file and with \_R12 for overlapped FASTQ file. Output file
name will also include MIG size threshold used.

**Settings**

The ``--mask <R1=[0,1]:R2=[0,1]>`` parameter indicates FASTQ
files to be assembled in paired-end data. By default both reads are
assembled. In case of ``0:0`` mask will process only overlapped reads.

The ``-c`` option indicates compressed output.

The ``-m`` option sets minimum number of reads in MIG. This should be
set according to Histogram script output to separate two peaks:
over-sequenced MIGs and erroneous MIGs that cluster around MIG size of
1.

To inspect the effect of such single-mismatch erroneous UMI sub-variants
see "collisions" output of Histogram script. Such collision events could
interfere with real MIGs when over-sequencing is relatively low. In this
case collisions could be filtered during MIG consensus assembly using
``--filter-collisions`` option.

V(D)J junction mapping
~~~~~~~~~~~~~~~~~~~~~~

CdrBlast-batch
^^^^^^^^^^^^^^

**Description**

A script to perform UMI-guided assembly

**Usage**

General:

.. code:: bash

    java -jar migec.jar CdrBlastBatch [options] -R gene [checkout_output_folder/ or .] [assemble_output_folder/ or .] output_folder

Performs CDR3 extraction and V/J segment determination for both raw
(**Checkout** output) and assembled-data. Gene parameter ``-R`` is
required unless metadata (``--sample-metadata``) is provided that
specifies gene for each sample; supported genes are *TRA*, *TRB*, *TRG*,
*TRD*, *IGH*, *IGK* and *IGL*. If either of *assembly\_output\_folder*
or *checkout\_output\_folder* is not specified, the processing will be
done solely for the remaining input, this is useful e.g. if one wants
quickly process the assembled data. Otherwise only samples and file
types (paired, overlapped or single) that are present in both outputs
will be used. Processing both raw and assembled data is required for
second stage error correction (removal of hot-spot errors).

Several default **CdrBlast** parameters could be set,

``--default-mask <R1=[0,1]:R2=[0,1]>`` - mask which specifies for which
read(s) in paired-end data to perform CDR3 extraction. In case of
``0:0`` mask will process only overlapped reads ``--default-species`` -
default species to be used for all samples, *human* (used by default) or
*mouse* ``--default-file-types`` - default file types (paired,
overlapped or single) to be processed for each sample. If several file
types are specified, the corresponding raw and assembled files will be
combined and used as an input to CdrBlast
``--default-quality-threshold <Phred=[2..40],CQS=[2..40]>`` - quality
threshold pair, default for all samples. First threshold in pair is used
for raw sequence quality (sequencing quality phred) and the second one
is used for assembled sequence quality (CQS score, the fraction of reads
in MIG that contain dominant letter at a given position) ``--no-sort`` -
no sorting is performed for output files which speeds up processing.
Could be safely used in full pipeline as FilterCdrBlastResults will
provide final clonotype table in sorted format

A sample metadata file could also be provided with
``--sample-metadata <file_name>`` argument to guide the batch CDR3
extraction. This file should have the following tab-separated table
structure:

+-------------+-----------+--------+----------------------+--------+--------------------------+
| Sample ID   | Species   | Gene   | File types           | Mask   | Quality threshold pair   |
+=============+===========+========+======================+========+==========================+
| S0          | human     | TRA    | paired, overlapped   | 1:0    | 25,30                    |
+-------------+-----------+--------+----------------------+--------+--------------------------+
| S1          | human     | TRB    | unpaired             | -      | 25,30                    |
+-------------+-----------+--------+----------------------+--------+--------------------------+
| S2          | mouse     | TRB    | paired               | 0:1    | 20,25                    |
+-------------+-----------+--------+----------------------+--------+--------------------------+

CdrBlast-manual
^^^^^^^^^^^^^^^

**Description**

A script to extract CDR3 sequences

**Usage**

General:

.. code:: bash

    java -jar migec.jar CdrBlast [options] -R gene file1.fastq[.gz] [file2.fastq[.gz] ...] output_file 

Standard, assuming an example of a library containing T-cell Receptor
Alpha Chain sequences

in case of MIG-assembled data:

.. code:: bash

    java -jar migec.jar CdrBlast -a -R TRA assembly/S1_R2.fastq.gz cdrblast/S1_asm.cdrblast.txt 

for raw data:

.. code:: bash

    java -jar migec.jar CdrBlast -R TRA checkout/S1_R2.fastq.gz cdrblast/S1_raw.cdrblast.txt

to concatenate and process two or more FASTQ files at once:

.. code:: bash

    java -jar migec.jar CdrBlast -R TRA checkout/S1_R2.fastq.gz checkout/S2_R2.fastq.gz cdrblast/S12_raw.cdrblast.txt

Gene parameter ``-R`` is required, supported genes are *TRA*, *TRB*,
*TRG*, *TRD*, *IGH*, *IGK* and *IGL*. Species could be provided with
``-S`` parameter, by default uses *HomoSapiens*, supported species are
*HomoSapiens*, *MusMusculus* and others. Assembled data should be passed
to the script with ``-a`` option. ``--same-sample`` option should be
used if several assembled files are provided from the same sample, so
duplicate UMIs will be discarded and not counted twice.

To get a sorted output use ``-o`` option, otherwise sorting will be
performed at **FilterCdrBlastResults** step. Note that both raw and
assembled data should be processed to apply the last step of filtration.

Result filtering
~~~~~~~~~~~~~~~~

FilterCdrBlastResults-batch
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description**

A script to filter erroneous CDR3 sequences produced due to hot-spot PCR
and NGS errors. It can also use a hybrid error correction method that
includes frequency-based filtering of singleton clonotypes (i.e.
clonotypes represeted by a single MIG).

**Usage**

General:

.. code:: bash

    java -jar migec.jar FilterCdrBlastResultsBatch [options] cdrblast_batch_folder/ output_folder/

Perform hot-spot error filtration for data process with
**CdrBlastBatch**. Options are the same as for manual version below.

FilterCdrBlastResults-manual
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Usage**

General:

.. code:: bash

    java -jar migec.jar FilterCdrBlastResults [options] cdrblast_result_assembled_data cdrblast_result_raw_data output_file

Example:

.. code:: bash

    java -jar migec.jar FilterCdrBlastResults cdrblast/S1_asm.cdrblast.txt cdrblast/S1_raw.cdrblast.txt final/S1.cdrblast.txt

The ``-s`` option tells to filter CDR3s represented by single MIGs. The
rationale for this is that the deep repertoire profiling (at least with
our protocol) can generate spurious singletons that are associated with
reverse transcription errors and experimental artifacts. Filtering is a
non-greedy procedure and filters single-MIG clonotypes only if a 1- or
2-mismatch parent clonotype exists at ratio 1:20 and 1:400 respectively.
This is done to preserve diversity for samples with shallow sequencing,
e.g. ran on MiSeq.

Other options:

-  ``-n`` - output non-coding clonotypes that contain either a stop
   codon or a frameshift within CDR3 sequence.

-  ``-c`` - include non canonical clonotypes that have a CDR3 region
   that does not start with conserved C residue, or end with a conserved
   F/W residue.

-  ``-r`` - sets the read accumulation threshold (default is ``1.0``)
   used for hot-spot error correction, see MiGEC paper for details.

Now the file *S1.cdrblast.txt* contains a filtered and sorted CDR3/V/J
clonotype table.
