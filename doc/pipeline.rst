The pipeline
------------

All routines in the pipeline are available in **manual** and **batch**
variants. Batch variants are designed to automatically handle several
input samples with minimal shell scripting glue between analysis steps, 
this is also the recommended way to use with MIGEC.

List of MIGEC batch routines:

- :ref:`checkoutbatch`
- :ref:`histogram`
- :ref:`assemblebatch`
- :ref:`cdrblastbatch`
- :ref:`cdrfinalbatch`

If the "barcodes" file is set properly, 
the entire pipeline can be written as following:

.. code-block:: bash

    MIGEC="java -Xmx8G -jar MIGEC-$VERSION.jar"
    $MIGEC CheckoutBatch -cu barcodes.txt checkout/
    $MIGEC Histogram checkout/ histogram/
    $MIGEC AssembleBatch -c checkout/ histogram/ assemble/
    $MIGEC CdrBlastBatch -R TRB checkout/ assemble/ cdrblast/
    $MIGEC FilterCdrBlastResultsBatch cdrblast/ cdrfinal/

Manual usage
~~~~~~~~~~~~

List of MIGEC manual routines:

- :ref:`checkoutmanual`
- :ref:`histogram`
- :ref:`assemblemanual`
- :ref:`cdrblastmanual`
- :ref:`cdrfinalmanual`

An example for a 300bp paired-end MiSeq run of IGH library on a 16Gb RAM
Unix server. Such sequencing read length allows complete IGH sequencing,
thus mate pairs overlap. First *barcodes.txt* should be created
containing adapter sequences, see the section below for guidelines.
Then, assuming that the corresponding FASTQ files are
*IGH\_SAMPLE\_R1.fastq.gz* and *IGH\_SAMPLE\_R2.fastq.gz*, UMI- and
multiplex index-containing adapter is near 5'UTR of V segment (so the
CDR3 is in mate#2 after reads are oriented) and NCBI-BLAST+ is
installed, run all 5 stages of the pipeline using the following command:

.. code-block:: bash

    $MIGEC Checkout -cute --overlap barcodes.txt IGH_S1-10_R1.fastq.gz IGH_S1-10_R2.fastq.gz checkout/
    $MIGEC Histogram checkout/ histogram/
    $MIGEC Assemble -c --mask 0:0 checkout/S1_R12.fastq.gz . assembly/
    $MIGEC CdrBlast -R IGH checkout/S1_R12.fastq.gz cdrblast/S1_raw.txt
    $MIGEC CdrBlast -a -R IGH assembly/S1_R12.fastq.gz cdrblast/S1_asm.txt
    $MIGEC FilterCdrBlastResults cdrblast/S1_asm.txt cdrblast/S1_raw.txt cdrfinal/S1.txt
    
.. note::
    
    As dot (`.`) is reserved by MIGEC to specify empty metadata fields and 
    unused arguments, use `./` in case you want to point to current directory.
    

Full-length immunoglobulin data analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MIGEC can be used as pre-processing tool to assemble full-length consensuses 
for further post-analysis with `HIgBlast <https://github.com/mikessh/higblast>`__ 
tool. Note that due to typically poor quality of MiSEQ 300+300bp reads, the 
``--overlap`` mode of Checkout routine is not guaranteed to perform well. Instead, we recommend 
to assemble consensuses first and then perform overlapping using 
external tools. For example, `MiTools <https://github.com/milaboratory/mitools>`__ 
``merge`` action can be used with the ``--same-strand`` option specified, the latter is 
critical as assembled consensuses are on the same strand in output in contrast to 
normal orientation of Illumina reads.

Consensus quality and overlap efficiency can be greatly improved using 
the ``--only-first-read`` option of Histogram and Assemble routines. If set, 
this option instructs routines to use only the first read that typically has higher quality 
than the second one. This applies to non-oriented reads and works better for 
asymmetric sequencing design, e.g. 400+200bp reads.
