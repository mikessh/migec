Consensus assembly
------------------

.. _assemblebatch:

Assemble-batch
~~~~~~~~~~~~~~

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

.. warning::

    In most cases, the automatic MIG size threshold selected by Histogram routine is ok. 
    However we strongly recommend manual inspection of Histogram output files and considering 
    to manually specify an appropriate MIG size threshold for input samples.

.. _assemblemanual:

Assemble-manual
~~~~~~~~~~~~~~~

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

.. note::
    
    To inspect the effect of such single-mismatch erroneous UMI sub-variants
    see "collisions" output of Histogram script. Such collision events could
    interfere with real MIGs when over-sequencing is relatively low. In this
    case collisions could be filtered during MIG consensus assembly using
    ``--filter-collisions`` option.