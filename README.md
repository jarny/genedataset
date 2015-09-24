genedataset
======
**genedataset** is a package to store and access gene expression datasets and gene definitions. It consists of two main classes, geneset and dataset.

## geneset
geneset stores gene information combined from both Ensembl and NCBI/Entrez (mouse and human only), so that you can query it:
```python
$ gs = geneset.Geneset().subset(queryStrings='ccr3')
$ print gs.geneIds()
 ['ENSG00000183625', 'ENSMUSG00000035448']
$ gs.dataframe()
 | EnsemblId          | Species     | EntrezId | GeneSymbol | Synonyms                     | Description                      | MedianTranscriptLength | Orthologue              |
 |--------------------|-------------|----------|------------|------------------------------|----------------------------------|------------------------|-------------------------|
 | ENSG00000183625    | HomoSapiens | 1232     | CCR3       | CC-CKR-3|CD193|CKR3|CMKBR3   | chemokine (C-C motif),receptor 3 | 1242.5                 | ENSMUSG00000035448:Ccr3 |
 | ENSMUSG00000035448 | MusMusculus | 12771    | Ccr3       | CC-CKR3|CKR3|Cmkbr1l2|Cmkbr3 | chemokine (C-C motif),receptor 3 | 3273                   | ENSG00000183625:CCR3    |
```

## dataset
dataset can store gene expression data so that it can be queried. The stored data consists of expression values (microarray and rna-seq) and sample data packaged into HDF5 format.
```python
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
```

## Full example with details
Here is a full example to download a microarray dataset from GEO using GEOparse package, and create a Dataset file.
Once the file has been created, it can be accessed through the Dataset instance. The advantage of this is to store
all related information for a dataset in one file, and gives you a python object that can be used for analyses and 
for application development.
```python
# Import the dataset from GEO (http://www.ncbi.nlm.nih.gov/geo) using GEOparse.
# You can also download the files from the website and parse it yourself.
$ import GEOparse
$ gse = GEOparse.get_GEO(geo='GSE68009')

# Gather various information needed to create the HDF5 file needed for Dataset instance
$ platformId = gse.metadata['platform_id'][0]  # only one platform in this dataset: 'GPL6246'

# Define dataset attributes dictionary 
$ attributes = {"fullname": "Sefik",
                "version": "1.0",   # version number for Dataset instance (customise for individual preference)
                "description": gse.metadata['summary'][0],  # "The colonic lamina propria contains..."
                "platform_type": 'microarray',   # 'microarray' or 'rna-seq'
                "platform_details": gse.gpls[platformId].metadata['title'][0],
                "pubmed_id": '26272906',
                "species": 'MusMusculus'}   # 'MusMusculus' or 'HomoSapiens'
    
# Define sample table (pandas DataFrame) with sampleIds as index and various sample groups as columns.
# Usually, some manual curation is needed here to select useful fields from metadata.
$ sampleData = []
$ for gsm_name, gsm in gse.gsms.iteritems():
	  surfaceMarkers = [item.split(':')[1].strip() for item in gsm.metadata['characteristics_ch1'] if item.startswith('phenotype markers')][0]
	  sampleData.append([gsm_name, gsm.metadata['title'][0], surfaceMarkers])
$ samples = pandas.DataFrame(sampleData, columns=['sampleId','name','surfaceMarkers']).set_index('sampleId')
    
# Define expression matrix (pandas DataFrame) with probeId as index and sampleIds as columns
$ expressionData = {}
$ for gsm_name, gsm in gse.gsms.iteritems():
	  expressionData[gsm_name] = gsm.table["VALUE"].tolist()
	  expression = pandas.DataFrame(expressionData, index=gse.gsms[gsm_name].table["ID_REF"].astype(str))  # ID_REF column looks like integers
		  # - assumes that index is identical for all samples
$ expression.index.name = 'probeId'   # expression.shape is (24922, 13)
    
# Define mapping between probes and genes (dict).
# Note that sometimes it may be easier to get this info from a different source, such as vendor website.
# First extract pandas DataFrame of probe ids and genes from gse.gpls
$ df = gse.gpls[platformId].table[["ID","gene_assignment"]].set_index("ID")
$ df.index = df.index.astype(str)
# Note that gene_assignment column is a bit of a mixed field which looks like this:
# "NM_008866 // Lypla1 // lysophospholipase 1 // 1 A1|1 // 18777 /// ENSMUST00000027036 // Lypla1 // lysophospholipase 1 // 1 A1|1 // 18777 /// BC013536 // Lypla1 // lysophospholipase 1 // 1 A1|1 // 18777"
# We need ensembl gene id - use genedataset.geneset to help with this
$ geneSymbolsFromProbeId = dict([(probeId,value.split("//")[1].strip()) for probeId,value in df["gene_assignment"].iteritems() if probeId in expression.index])
$ gs = geneset.Geneset().subset(queryStrings=set(geneSymbolsFromProbeId.values()), searchColumns=["GeneSymbol"], caseSensitive=True, matchSubstring=False)
$ gsdf = gs.dataframe()
    
$ genesFromProbes = {}   # will look like {'10345442': ['ENSMUSG00000045216'], ...}
$ for probeId,geneSymbol in geneSymbolsFromProbeId.iteritems():
	  index = gsdf[gsdf['GeneSymbol']==geneSymbol].index
	  if len(index)>0: genesFromProbes[probeId] = index.tolist()

$ probesFromGenes = dict([(geneId,[]) for geneId in set(sum(genesFromProbes.values(),[]))])
	  # - will look like {'ENSMUSG00000045216':['10345442'], ...}
$ for probeId,geneIds in genesFromProbes.iteritems():
	  for geneId in geneIds:
		  probesFromGenes[geneId].append(probeId)

$ ds = dataset.createDatasetFile(
	  "../DataOutput",
	  name="sefik", 
	  attributes=attributes, 
	  samples=samples,
	  expression=expression, 
	  probeIdsFromGeneId=probesFromGenes,
	  geneIdsFromProbeId=genesFromProbes)
```

## Contact
#### Jarny Choi, Walter + Eliza Hall Institute
* e-mail: jchoi@wehi.edu.au

## Changes 
* v0.1.x - Initial release with minor adjustments to test pypi and github upload/download.

## License
[MIT License](LICENSE.txt)

