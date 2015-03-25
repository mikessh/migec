.. _histogram:

MIG statistics
--------------

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