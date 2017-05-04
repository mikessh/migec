De-multiplexing
---------------

.. _checkoutbatch:

Checkout-batch
~~~~~~~~~~~~~~

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

**Output format**

The Checkout routine produces files in the FASTQ format that have a specific
``UMI`` field added to the header. Each read successfully matched by Checkout
will be output as follows:

.. code::

    @ILLUMINA_HEADER UMI:NNNN:QQQQ
    ATAGATTATGAGTATG
    +
    ##II#IIIIIIIIIII

The original read header (``ILLUMINA_HEADER`` here) is preserved, the
appended ``UMI:NNNN:QQQQ`` contains the sequence of the UMI tag (``NNNN`` bases)
and its quality string (``QQQQ``).

.. _checkoutmanual:

Checkout-manual
~~~~~~~~~~~~~~~

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
drive space.

.. note::

    When there is a huge number of unassigned/unused reads ``--skip-undef`` option
    greatly speeds up de-multiplexing. However, take care to carefully investigate
    the reasons behind low barcode extraction rate if it is a case.

.. important::

    The ``--overlap`` option may not perform well for poor quality reads, which is
    a typical situation for 300+300bp MiSEQ sequencing. In this case, merging reads
    using external software after Assemble stage is recommended.
