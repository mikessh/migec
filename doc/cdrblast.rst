V(D)J junction mapping
----------------------

.. _cdrblastbatch:

CdrBlast-batch
~~~~~~~~~~~~~~

**Description**

A script to extract CDR3 sequences. Will properly combine reads coming from 
paired and overlapped data and perform analysis for both raw and assembled data.

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

**Usage**

General:

.. code:: bash

    java -jar migec.jar CdrBlastBatch [options] -R gene [checkout_output_folder/ or .] [assemble_output_folder/ or .] output_folder

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

+-------------+-----------+---------+----------------------+--------+--------------------------+
| Sample ID   | Species   | Gene    | File types           | Mask   | Quality threshold pair   |
+=============+===========+=========+======================+========+==========================+
| S0          | human     | TRA     | paired, overlapped   | 1:0    | 25,30                    |
+-------------+-----------+---------+----------------------+--------+--------------------------+
| S1          | human     | TRB     | unpaired             | -      | 25,30                    |
+-------------+-----------+---------+----------------------+--------+--------------------------+
| S2          | mouse     | TRB,TRA | paired               | 0:1    | 20,25                    |
+-------------+-----------+--------=+----------------------+--------+--------------------------+

See section below for more details.

.. _cdrblastmanual:

CdrBlast-manual
~~~~~~~~~~~~~~~

**Description**

A script to map V-(D)-J junctions, extract CDR3 sequences and assemble clonotypes.

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
*TRG*, *TRD*, *IGH*, *IGK* and *IGL*. Several chains can be specified, 
for example ``-R TRA,TRB`` or ``-R IGH,IGL``. Species could be provided with
``-S`` parameter, by default uses *HomoSapiens*, supported species are
*HomoSapiens*, *MusMusculus* and others. Assembled data should be passed
to the script with ``-a`` option. ``--same-sample`` option should be
used if several assembled files are provided from the same sample, so
duplicate UMIs will be discarded and not counted twice.

To get a sorted output use ``-o`` option, otherwise sorting will be
performed at **FilterCdrBlastResults** step. Note that both raw and
assembled data should be processed to apply the last step of filtration.

.. note::

    In order to use all alleles, not just the major (*01 ones), use the
    ``--all-alleles`` option. To include non-coding segments (V segment 
    pseudogenes) use the ``--all-segments`` option.