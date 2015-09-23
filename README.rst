genedataset
===========

**genedataset** is a package to store and access gene expression datasets
and gene definitions. It consists of two main classes, geneset and
dataset.

geneset stores gene information combined from both Ensembl and
NCBI/Entrez (mouse and human only), so that you can query it:

::

    $ gs = geneset.Geneset().subset(queryStrings='ccr3')
    $ print gs.geneIds()
     ['ENSG00000183625', 'ENSMUSG00000035448']

dataset can store gene expression data so that it can be queried. The
stored data consists of expression values (microarray and rna-seq) and
sample data packaged into HDF5 format.

::

    $ ds = dataset.Dataset("mydataset.h5")
    $ print ds.platform_type
     'rna-seq'

More details to come.

Version
-------

-  Version 0.1.4

Contact
-------

Jarny Choi, Walter + Eliza Hall Institute
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  e-mail: jchoi@wehi.edu.au
