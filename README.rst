========
Overview
========

SCbrowse is an interactive browser for single-cell ATAC-seq data.
It allows to simultaneously explore the accessibility patterns
of selected cells in the embedding space and for selected genomic regions.
The tools freely available under a GNU Lesser General Public License v3 or later (LGPLv3+).

Installation
============

The following lines show how to install the requirements and scbrowse into
a new conda environment.

::

    conda create -n scbrowse python=3.7
    conda activate scbrowse
    conda install -c bioconda bedtools

    pip install janggu[tf]
    pip install https://github.com/BIMSBbioinfo/scregseg/archive/master.zip


    git clone https://github.com/BIMSBbioinfo/scbrowse
    pip install -e scbrowse[gunicorn]

::

Data preparation
================

SCbrowse needs two required ingredients:

1. A genome-wide count matrix
2. 2D embedding of the cells
3. Gene annotation in bed format 

The count matrix can be created from a BAM-file.
For example for a 1000 kp resolution countmatrix use:

::

    scregseg make_tile --bamfile <bam> --resolution 1000 --bedfile <outputbed>
    scregseg bam_to_counts --bamfile <bam> --counts <countmatrix> --regions <outputbed>

In order to obtain a 2D embedding of the cells,
a number of tools have been proposed already.
For instance, one might want to follow cisTopic, snapATAC, Seurat, etc.
The tutorials of the respective tools illustrate how UMAP embeddings
of the cells can be obtains.
Eventually, scbrowse expects a tsv table with the format

::

   barcode      dim1    dim2    annot.annot1    annot.annot2   ...

The first three columns (barcode, dim1 and dim2) are required.
An arbitrary number of additional columns with prefix `annot.`
can be added with categorical labels. They can be selected for
coloring the cells.

Using scbrowse (development mode)
=================================

To launch scbrowse in development mode, use first set the environment variables
to point to the input files

::

    export SCBROWSE_EMBEDDING=<embedding.tsv>
    export SCBROWSE_MATRIX=<matrix.npz>
    export SCBROWSE_REGIONS=<regions.bed>
    export SCBROWSE_GENES=<genes.bed>
    export SCBROWSE_LOGS=<scbrowse.log>
    scbrowse


Afterwards you can browse the data locally in a web-browser by opening
https://localhost:8051

Deploy scbrowse in production mode
==================================

The simplest way to deploy scbrowse using gunicorn is
to adjust the :code:`startup.sh`.
It requires to specify the input file locations using a set of
environment variables.
Subsequently, it launches the application using gunicorn.

::

    sh startup.sh

Afterwards you can browse the data locally in a web-browser by opening
https://localhost:8000


