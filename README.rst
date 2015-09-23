genedataset
===========

**genedataset** is a package to store and access gene expression
datasets and gene definitions. It consists of two main classes, geneset
and dataset.

geneset
-------

geneset stores gene information combined from both Ensembl and
NCBI/Entrez (mouse and human only), so that you can query it:

::

    $ gs = geneset.Geneset().subset(queryStrings='ccr3')
    $ print gs.geneIds()
     ['ENSG00000183625', 'ENSMUSG00000035448']
    $ gs.dataframe()
     | EnsemblId          | Species     | EntrezId | GeneSymbol | Synonyms                     | Description                      | MedianTranscriptLength | Orthologue              |
     |--------------------|-------------|----------|------------|------------------------------|----------------------------------|------------------------|-------------------------|
     | ENSG00000183625    | HomoSapiens | 1232     | CCR3       | CC-CKR-3|CD193|CKR3|CMKBR3   | chemokine (C-C motif),receptor 3 | 1242.5                 | ENSMUSG00000035448:Ccr3 |
     | ENSMUSG00000035448 | MusMusculus | 12771    | Ccr3       | CC-CKR3|CKR3|Cmkbr1l2|Cmkbr3 | chemokine (C-C motif),receptor 3 | 3273                   | ENSG00000183625:CCR3    |

dataset
-------

dataset can store gene expression data so that it can be queried. The
stored data consists of expression values (microarray and rna-seq) and
sample data packaged into HDF5 format.

::

    $ ds = dataset.Dataset("genedataset/data/testdataset.h5")
    $ ds
     <Dataset name:testdata species:MusMusculus, platform_type:microarray>
    $ ds.expressionMatrix()
     | probeId | s01  | s02  | s03  | s04  |
     |---------|------|------|------|------|
     | probe1  | 3.45 | 4.65 | 2.65 | 8.23 |
     | probe2  | 5.54 | 0.00 | 1.43 | 6.43 |
     | probe3  | 0.00 | 0.00 | 4.34 | 5.44 |
    $ ds.sampleTable()
     | sampleId | celltype | tissue |
     |----------|----------|--------|
     | s01      | B1       | BM     |
     | s02      | B1       | BM     |
     | s03      | B2       | BM     |
     | s04      | B2       | BM     |

Contact
-------

Jarny Choi, Walter + Eliza Hall Institute
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  e-mail: jchoi@wehi.edu.au

Changes
-------

-  v0.1.x - Initial release with minor adjustments to test pypi and github upload/download.

License
-------

`MIT License`_

.. _MIT License: LICENSE.txt
