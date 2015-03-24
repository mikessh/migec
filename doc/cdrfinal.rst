CdrFinal
--------

Batch result filtering
~~~~~~~~~~~~~~~~~~~~~~

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

Manual result filtering
~~~~~~~~~~~~~~~~~~~~~~~

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
