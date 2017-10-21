genedataset
======
**genedataset** is a package to store and access gene expression datasets and gene definitions. It consists of two main classes, geneset and dataset.

## Changes in version 1.0
Some significant changes have been made in this version:
1. "MedianTranscriptLength" property of Geneset has been replaced with "TranscriptLengths", which holds a list of transcript lengths, one per transcript (eg. [2000,1530]). Note that the value contained is a string, so convert into a list of integers before using the value.
2. The Gene annotation has been upgraded to Ensembl version 88.
3. Dataset no longer supports microarrays, so it's been simplified.
4. The package supports both python 3 and 2 - tested on 2.7.14 and 3.6.3.

## geneset
geneset stores gene information combined from both Ensembl and NCBI/Entrez (mouse and human only), so that you can query it:
```python
from genedataset import geneset
gs = geneset.Geneset().subset(queryStrings='ccr3')
print gs.geneIds()
 ['ENSG00000183625', 'ENSMUSG00000035448']
gs.dataframe()
 | EnsemblId          | Species     | EntrezId | GeneSymbol | Synonyms                     | Description                      | TranscriptLengths                             | Orthologue              |
 |--------------------|:-----------:|---------:|-----------:|-----------------------------:|---------------------------------:|----------------------------------------------:|------------------------:|
 | ENSG00000183625    | HomoSapiens | 1232     | CCR3       | CC-CKR-3|CD193|CKR3|CMKBR3   | C-C motif chemokine receptor 3   | [2000, 1581, 400, 436, 212, 1284, 1201, 1786] | ENSMUSG00000035448:Ccr3 |
 | ENSMUSG00000035448 | MusMusculus | 12771    | Ccr3       | CC-CKR3|CKR3|Cmkbr1l2|Cmkbr3 | chemokine (C-C motif) receptor 3 | [3272]                                        | ENSG00000183625:CCR3    |
```

## dataset
dataset can store gene expression data so that it can be queried. The stored data consists of expression values (usually rna-seq) and sample data packaged into HDF5 format.
```python
from genedataset import dataset
ds = dataset.Dataset("genedataset/data/testdataset.h5")
ds
 <Dataset name:testdata 4 samples>
ds.expressionMatrix()
 | featureId | s01  | s02  | s03  | s04  |
 |---------|------|------|------|------|
 | gene1  | 3.45 | 4.65 | 2.65 | 8.23 |
 | gene2  | 5.54 | 0.00 | 1.43 | 6.43 |
 | gene3  | 0.00 | 0.00 | 4.34 | 5.44 |
ds.sampleTable()
 | sampleId | celltype | tissue |
 |----------|----------|--------|
 | s01      | B1       | BM     |
 | s02      | B1       | BM     |
 | s03      | B2       | BM     |
 | s04      | B2       | BM     |
```

## Dataset creation example
Here is an example to create a Dataset file from text files. Once the file has been created, it can be accessed through the Dataset instance. The advantage of this is to store all related information for a dataset in one file, and gives you a python object that can be used for analyses and for application development.
```python
import pandas
from genedataset import geneset, dataset

# Define dataset attributes dictionary 
attributes = {"name":"mydataset",
			  "fullname": "MyDataset",
			  "version": "1.0",
			  "description": "RNA-seq expression of murine haematopoietic cells",
			  "expression_type": "cpm",
			  "pubmed_id": None,
			  "species": "MusMusculus"}
					  
# Read expression matrix and sample table
expression = pandas.read_csv("cpm.txt", sep="\t", index_col=0)   # data frame must have index
samples = pandas.read_csv("samples.txt", sep="\t", index_col=0)  # sample ids form index, should match columns of expression (different ordering OK)

ds = dataset.createDatasetFile("/datasets", attributes=attributes, samples=samples, expression=expression)
        
```

## Contact
#### Jarny Choi, University of Melbourne
* e-mail: jarnyc@unimelb.edu.au

## Changes 
* v1.0 - Major upgrade.
* v0.1.x - Initial release with minor adjustments to test pypi and github upload/download.

## License
[MIT License](LICENSE.txt)

