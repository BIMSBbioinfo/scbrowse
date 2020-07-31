========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - |
        |
    * - package
      - | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/python-scbrowse/badge/?style=flat
    :target: https://readthedocs.org/projects/python-scbrowse
    :alt: Documentation Status

.. |commits-since| image:: https://img.shields.io/github/commits-since/wkopp/python-scbrowse/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/wkopp/python-scbrowse/compare/v0.0.0...master



.. end-badges

SCbrowse is an interactive browser for single-cell ATAC-seq data.
It allows to simultaneously explore the accessibility patterns
of selected cells in the embedding space and for selected genomic regions.
The tools freely available under a GNU Lesser General Public License v3 or later (LGPLv3+).

Installation
============

Before installing scbrowser, please install bedtools and scregseg.

Scregseg can be install using

::

    pip install janggu[tf]
    pip install https://github.com/BIMSBbioinfo/scregseg/archive/master.zip


::

#    pip install scbrowse
#
#You can also install the in-development version with::

    pip install https://github.com/wkopp/python-scbrowse/archive/master.zip


Data preparation
================

SCbrowse needs two required ingredients:

1. A genome-wide count matrix
2. 2D embedding of the cells

The count matrix can be created from a BAM-file.
For example for a 1000 kp resolution countmatrix use:

::

    scregseg make_tile --bamfile <bam> --resolution 500 --bedfile <outputbed>
    scregseg bam_to_counts --bamfile <bam> --counts <countmatrix> --regions <outputbed>

In order to obtain a 2D embedding of the cells,
a number of tools have been proposed already.
For instance, one might want to follow cisTopic, snapATAC, Seurat, etc.
The tutorials of the respective tools illustrate how UMAP embeddings
of the cells can be obtains.
Eventually, scbrowse expects a tsv table with the format

::

   barcode      dim1    dim2    annot.annot1    annot.annot2

The first three columns (barcode, dim1 and dim2) are required.
An arbitrary number of additional columns with prefix `annot.`
can be added with categorical labels. They can be selected for
coloring the cells.

Using scbrowse
=================

scbrowse can be launched according to

::

    scbrowse -embedding <embedding.tsv> -matrix <countmatrix.npz> -bed <tile.bed> -genes <genes.bed>


Afterwards you can browse the data in a web-browser by opening
https://localhost:8051


Documentation
=============


https://python-scbrowse.readthedocs.io/


