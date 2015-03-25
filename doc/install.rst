Installing and running
----------------------

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

.. note:: 

    The data from 454 platform should be used with caution, as it contains 
    homopolymer errors which (in present framework) result in reads dropped
    during consensus assembly. The 454 platform has a relatively low read
    yield, so additional read dropping could result in over-sequencing level
    below required threshold. If you still wish to give it a try, we would
    recommend filtering off all short reads and repairing indels with
    `Coral <http://www.cs.helsinki.fi/u/lmsalmel/coral/>`__, the latter
    should be run with options ``-mr 2 -mm 1000 -g 3``.

.. warning:: 

    NCBI-BLAST+ package is required. Could be directly installed on Linux
    using a command like $sudo apt-get ncbi-blast+ or downloaded and
    installed directly from here: 
    ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

.. warning:: 

    Consider providing sufficient memory for the pipeline, i.e. 8Gb for
    MiSeq or 36Gb for HiSeq sample, depending on sample sequence diversity
    and current script (CdrBlast requires has the highest memory
    requirements). To do so, execute the script with ``-Xmx`` argument: 
    ``java -Xmx8G -jar MIGEC-$VERSION.jar CdrBlast [arguments]``. 
    If insufficient amount memory is allocated, the Java Virtual Machine
    could drop with a *Java Heap Space Out of Memory* error.
