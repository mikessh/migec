## Various utilities

This folder contains a number of scripts for processing MIGEC output:

* ``GroomConsensusFastq.groovy input.fastq[.gz] output.fastq[.gz]`` will ensure compatibility of assembled consensus FASTQ file by grooming headers and correcting read orientation (use ``-r`` flag for second read FASTQ for this purpose).

* ``histogram.R`` and ``pwm.R`` are scripts to produce some fancy graphics for ``Histogram`` routine output.

* ``migec_summary.Rmd`` is the MIGEC summary template (same as the one used by ``Report`` routine).