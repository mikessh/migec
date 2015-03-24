MIGEC: Molecular Identifier Group-based Error Correction pipeline
=================================================================

FEATURES
~~~~~~~~

-  De-multiplexing, adapter trimming and read overlapping for NGS data.
   Extraction of UMI sequences

-  Assembly of consensuses for original molecules that entered library
   preparation by grouping reads with identical molecular identifiers

-  Mapping of V, D and J segments, extraction of CDR3 regions and
   clonotype assembly for all human and mouse immune receptors 
   (TRA/TRB/TRG/TRD and IGH/IGL/IGK)

-  Additional filtering of hot-spot errors in CDR3 sequence

-  Flexible and straightforward batch processing

-  Currently all species-gene pairs that have germline sequences (based
   on `IMGT <http://www.imgt.org/>`__) that allow CDR3 identification
   are supported

+------------------------+-------------------------------------+
| Species                | Gene                                |
+========================+=====================================+
| HomoSapiens            | TRA, TRB, TRG, TRD, IGL, IGK, IGH   |
+------------------------+-------------------------------------+
| MusMusculus            | TRB, TRG, TRD, IGL, IGK, IGH        |
+------------------------+-------------------------------------+
| MacacaMulatta          | TRB, IGK, IGH                       |
+------------------------+-------------------------------------+
| OryctolagusCuniculus   | IGL, IGK, IGH                       |
+------------------------+-------------------------------------+
| RattusNorvegicus       | IGL, IGH                            |
+------------------------+-------------------------------------+
| CanisLupusFamiliaris   | TRB, TRG                            |
+------------------------+-------------------------------------+
| SusScrofa              | IGL, IGK                            |
+------------------------+-------------------------------------+
| BosTaurus              | TRD                                 |
+------------------------+-------------------------------------+
| MusSpretus             | IGL                                 |
+------------------------+-------------------------------------+
| GallusGallus           | TRB                                 |
+------------------------+-------------------------------------+
| AnasPlatyrhynchos      | TRB                                 |
+------------------------+-------------------------------------+

INSTALLATION AND RUNNING
~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline is written in Groovy (a Java scripting language) and
distributed as an executable JAR. To install it get the latest
`JRE <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`__
and download the executable from `releases
section <https://github.com/mikessh/MIGEC/releases>`__.

To ran a specific script from the pipeline, say **Checkout**, execute

.. code:: bash

    java -jar MIGEC-$VERSION.jar Checkout [arguments]

Where ``$VERSION`` stands for pipeline version (e.g. 1.2.1), this notation is 
omitted in MIGEC routine documentation. 
To view the list of available scripts execute:

.. code:: bash

    java -jar MIGEC-$VERSION.jar -h

alternatively you can download the repository and compile it from source
using `Maven <http://maven.apache.org/>`__ (requires Maven version 3.0)

.. code:: bash

    git clone https://github.com/mikessh/MIGEC.git
    cd MIGEC/
    mvn clean install
    java -jar target/MIGEC-$VERSION.jar

NOTE
~~~~

The data from 454 platform should be used with caution, as it contains
homopolymer errors which (in present framework) result in reads dropped
during consensus assembly. The 454 platform has a relatively low read
yield, so additional read dropping could result in over-sequencing level
below required threshold. If you still wish to give it a try, we would
recommend filtering off all short reads and repairing indels with
`Coral <http://www.cs.helsinki.fi/u/lmsalmel/coral/>`__, the latter
should be run with options ``-mr 2 -mm 1000 -g 3``.

IMPORTANT
~~~~~~~~~

NCBI-BLAST+ package is required. Could be directly installed on Linux
using a command like $sudo apt-get ncbi-blast+ or downloaded and
installed directly from here:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Consider providing sufficient memory for the pipeline, i.e. 8Gb for
MiSeq or 36Gb for HiSeq sample, depending on sample sequence diversity
and current script (CdrBlast requires has the highest memory
requirements). To do so, execute the script with ``-Xmx`` argument:

::

    java -Xmx8G -jar MIGEC-$VERSION.jar CdrBlast [arguments]

If insufficient amount memory is allocated, the Java Virtual Machine
could drop with a *Java Heap Space Out of Memory* error.

THE PIPELINE
------------

All routines in the pipeline are available in "manual" and "batch"
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
    
MANUAL USAGE EXAMPLE
--------------------

An example for a 300bp paired-end MiSeq run of IGH library on a 16Gb RAM
Unix server. Such sequencing read length allows complete IGH sequencing,
thus mate pairs overlap. First *barcodes.txt* should be created
containing adapter sequences, see **Checkout** section for guidelines.
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

Table of contents:
==================

.. toctree::
   :maxdepth: 1
   checkout
   histogram
   assemble
   cdrblast
   cdrfinal
   post

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

