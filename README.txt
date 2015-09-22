genedataset
===========

genedataset is a package to store and access gene expression datasets and gene definitions. It consists of two main classes, geneset and dataset.

geneset can be used to handle gene data:
>> import geneset
>> gs = geneset.Geneset().subset(queryStrings='ccr3')
>> print gs.geneIds()
>> ['ENSG00000183625', 'ENSMUSG00000035448']

dataset can be used to handle expression datasets:
>> import dataset
>> ds = dataset.Dataset("mydataset.h5")
>> print ds.platform_type
>> 'rna-seq'

These classes interact with stored data inside the data/ sub directory.
