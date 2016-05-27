[![Build Status](https://travis-ci.org/mikessh/migec.svg?branch=master)](https://travis-ci.org/mikessh/migec)

# MiGEC: Molecular Identifier Guided Error Correction pipeline  

This pipeline provides several useful tools for analysis of immune repertoire sequencing data. Its main feature is the ability to use information from unique nucleotide tags (UMIs, see this [paper](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html) for details), which are attached to molecules before sequencing library preparation and allow to backtrack the original sequence of molecule. UMIs make it possible to computationally filter nearly all experimental errors from resulting immune receptor sequences. 

This pipeline was designed for libraries sequenced using Illumina MiSeq and HiSeq and the main requirement for sequencing reads is that they should contain the entire CDR3 region of immune receptor gene. Sequencing libraries with high over-sequencing, i.e. ones that have 5+ reads per starting molecule (unique UMI tag), should be used for optimal error elimination.

Several modules of the pipeline, such as de-multiplexing and CDR3 extraction could be utilized for a wider range of datasets.

Compiled binaries are available from [here](https://github.com/mikessh/migec/releases/latest).
The software is cross-platform and requires Java v1.8 to run.
Easy installation on MacOS/Linux via [Homebrew](http://brew.sh/) or [Linuxbrew](http://linuxbrew.sh/):
```bash
brew tap mikessh/repseq
brew install migec
migec Checkout ...
```

For more details please see the [paper](http://www.nature.com/nmeth/journal/v11/n6/abs/nmeth.2960.html) describing MiGEC.

Full documentation is provided via [ReadTheDocs](http://migec.readthedocs.org/en/latest/index.html). You might be also interested in taking the following [tutorial](http://repseq-tutorial.readthedocs.org/en/latest/).

Please cite the tool as:

> Shugay M *et al.* Towards error-free profiling of immune repertoires. Nature Methods **11**, 653–655 (2014)