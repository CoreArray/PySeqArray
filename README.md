PySeqArray: Python Interface to SeqArray Files for Data Management of Whole-Genome Sequence Variant Calls
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)

pre-release version


## Features

Big data management of whole-genome sequence variant calls with thousands of individuals: genotypic data (e.g., SNVs, indels and structural variation calls) and annotations in SeqArray files are stored in an array-oriented and compressed manner, with efficient data access using the Python programming language.

The SeqArray format is built on top of Genomic Data Structure (GDS) data format, and defines required data structure. GDS is a flexible and portable data container with hierarchical structure to store multiple scalable array-oriented data sets. It is suited for large-scale datasets, especially for data which are much larger than the available random-access memory. It also offers the efficient operations specifically designed for integers of less than 8 bits, since a diploid genotype usually occupies fewer bits than a byte. Data compression and decompression are available with relatively efficient random access.


## Installation

```sh
pip install git+git://github.com/CoreArray/pygds.git
pip install git+git://github.com/CoreArray/PySeqArray.git
```


## Citation

#### Original paper (implemented in an [R/Bioconductor](http://bioconductor.org/packages/SeqArray) package):

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## SeqArray File Download

* [1000 Genomes Project](http://bochet.gcc.biostat.washington.edu/seqarray/1000genomes)


## Examples

More examples will be posted here.

