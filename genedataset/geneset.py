"""
Geneset class is a list of information about genes. It is not currently a collection of individual
"Gene" objects, because I find most controller code is concerned about fetching a collection of genes
while view code is often concerned about individual gene properties (to show in a table for example),
but by this stage we can just pass a json object anyway.

object with self.dataframe which contains information about each gene per row. This data frame has columns:
['EnsemblId','Species','EntrezId','GeneSymbol','Synonyms','Description','MedianTranscriptLength','Orthologue']
where 'EnsemblId' is the index.
"""
import pandas, os, sys, re, copy
from . import dataDirectory

# ------------------------------------------------------------
# Geneset class
# ------------------------------------------------------------
class Geneset(object):
	"""
	Defines a Geneset object, which has name, description and a collection of gene properties.
	There is also metadata associated with the Geneset which specifies origins of data and 
	"""
	def __init__(self, name='unnamed', description=''):
		self._datafile = '%s/Genes.h5' % dataDirectory()
		self._dataframe = pandas.read_hdf(self._datafile, 'data')
		self._metadata = pandas.read_hdf(self._datafile, 'metadata')

		self.name = name
		self.description = description
	
	def __repr__(self):
		return "<Geneset name:%s size:%s>" % (self.name, self.size())
		
	def size(self):
		# Better to have this as a method rather than property, as self._dataframe may change
		"""
		Return the size of this Geneset, which is equivalent to the number of unique gene ids.
		"""
		return self._dataframe.shape[0]
		
	def geneIds(self):
		"""
		Return a list of gene ids found in this Geneset.
		"""
		return self._dataframe.index.tolist()
		
	def geneSymbols(self, returnType="list"):
		"""
		Return a list of gene symbols in this geneset if returnType=="list".
		If returnType=="dict", it returns a dictionary of gene symbols keyed on gene id
		"""
		return self._dataframe['GeneSymbol'].to_dict() if returnType=="dict" else self._dataframe['GeneSymbol'].tolist()
		
	def species(self, ignoreMixed=True):
		"""
		Return species for this Geneset.

		Parameters:
			ignoreMixed: boolean. If True, the method will only look at the first gene, so this is more efficient
				if you already know that this Geneset comprises of one species. If False, the method will look at
				all species and return a list.
			 
		Returns:
			a string, either of {"MusMusculus", "HomoSapiens"} if ignoreMixed==True
			a list, such as ["MusMusculus"] or ["MusMusculus","HomoSapiens"] if ignoreMixed==False.
		"""
		if ignoreMixed:
			return self._dataframe['Species'].iloc[0] if self.size()>0 else None
		else:
			return list(set(self._dataframe['Species']))
		
	def subset(self, *args, **kwargs):
		"""
		Return a copy of this instance, but with a subset of genes.
		
		Parameters:
			list of query strings: such as Geneset().subset(['myb','ccr3']) or keyworded arguments
			
			queryStrings: list of query strings such as ['myb','ENSMUSG00000039601']
			species: one of ['MusMusculus','HomoSapiens'] to restrict search space; may be left out;
			caseSensitive: boolean, default False;
			searchColumns: a list whose members may be any of ['EnsemblId','EntrezId','GeneSymbol','Synonyms','Description'].
				Leaving it out will use this full list. ('GeneId' can also be used instead of 'EnsemblId')
			matchSubstring: boolean, default True; If True, 'Ccr' will match 'Ccr2' and 'Ccr3' for example;
				Use caseSensitive=True, matchSubstring=False to create exact match conditions
			
		Returns:
			A new Geneset instance with the subset of genes.
			
		Examples:
			>>> print Geneset().subset(['ENSMUSG00000039601','ccr3']).geneSymbols()
			EnsemblId
			ENSG00000183625       CCR3
			ENSMUSG00000035448    Ccr3
		"""
		gs = copy.copy(self)
		queryStrings = kwargs['queryStrings'] if 'queryStrings' in kwargs else args[0] if args else []
		if isinstance(queryStrings, str) or isinstance(queryStrings, unicode):	# assume queryString was specified as a string
			queryStrings = [queryStrings]
	
		if len(queryStrings)==0:
			return gs
			
		df = gs._dataframe
		
		if 'species' in kwargs and kwargs['species'] in ['MusMusculus','HomoSapiens']:
			df = df[df['Species']==kwargs['species']]
			
		caseSensitive = kwargs.get('caseSensitive', False)
		if not caseSensitive: queryStrings = [item.lower() for item in queryStrings]
		
		matchSubstring = kwargs.get('matchSubstring', True)
		
		# determine which columns to search
		searchColumns = kwargs.get('searchColumns')
		allColumns = ['EnsemblId','EntrezId','GeneSymbol','Synonyms','Description']
		if searchColumns:
			searchColumns = ['EnsemblId' if item=='GeneId' or item=='geneId' else item for item in searchColumns]
		if searchColumns and len(set(allColumns).intersection(set(searchColumns)))>0:
			searchColumns = list(set(allColumns).intersection(set(searchColumns)))
		else:
			searchColumns = allColumns
		
		rowsToKeep = set()
		df = df.reset_index()
		for column in searchColumns:
			for rowIndex,value in df[column].iteritems():
				if rowIndex in rowsToKeep or not value: continue
				if not caseSensitive: value = value.lower()
				for queryString in queryStrings:
					if matchSubstring and queryString in value or (not matchSubstring and queryString==value):
						rowsToKeep.add(rowIndex)
						break

		gs._dataframe = df.loc[list(rowsToKeep),:].set_index('EnsemblId')
		return gs
		
	def orthologueGeneIds(self):
		"""
		Return a list of gene ids comprising of orthologues of the current set.
		"""
		geneIds = []
		for geneId,row in self._dataframe.iterrows():
			for item in row['Orthologue'].split(','):	# looks like 'ENSG00003435:Gene1,ENSG00002525:Gene2' (multiple orthologues possible)
				if item.split(':')[0]: geneIds.append(item.split(':')[0])
		return list(set(geneIds))

	def to_json(self):
		"""
		Return an array of dictionaries for all genes in this geneset. 
		Example:
			>>> gs.to_json()
			[{"EnsemblId":"ENSMUSG00000047591",
			  "Species":"MusMusculus",
			  "EntrezId":"378435",
			  "GeneSymbol":"Mafa",
			  "Synonyms":"RIPE3b1",
			  "Description":"v-maf musculoaponeurotic fibrosarcoma oncogene family, protein A (avian) ",
			  "MedianTranscriptLength":1080.0,
			  "Orthologue":"ENSG00000182759:MAFA"},...]
		"""
		return self._dataframe.reset_index().to_json(orient="records")					

	def dataframe(self):
		"""
		Return a pandas DataFrame object corresponding to gene information in this Geneset.
		Index of the DataFrame returned is Ensembl gene id. Columns are:
		['Species','EntrezId','GeneSymbol','Synonyms','Description','MedianTranscriptLength','Orthologue']
		"""
		return self._dataframe
		
	'''
	def medianTranscriptLengths(self):
		"""Return a dictionary of median transcript length keyed on gene id
		"""
		return dict(zip(self._dataframe.index, self._dataframe['MedianTranscriptLength']))
	'''
	
'''
Create the Genes.h5 file from source files:
	Ensembl files downloaded from http://ensembl.org/biomart
	Entrez files downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/
	Mouse-human orthologue file downloaded from ftp://ftp.informatics.jax.org/pub/reports/index.html#homology

Format of the Ensembl files: (Be careful when selecting attributes to download from biomart - columns are in the same order as selection.)
	Ensembl Gene ID	EntrezGene ID	Associated Gene Name	Description
	ENSMUSG00000064372		mt-Tp	mitochondrially encoded tRNA proline [Source:MGI Symbol;Acc:MGI:102478]
	ENSMUSG00000064371		mt-Tt	mitochondrially encoded tRNA threonine [Source:MGI Symbol;Acc:MGI:102473]
	ENSMUSG00000064370	17711	mt-Cytb	mitochondrially encoded cytochrome b [Source:MGI Symbol;Acc:MGI:102501]

Format of the Entrez files:
	#Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description type_of_gene Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority Nomenclature_status Other_designations Modification_date (tab is used as a separator, pound sign - start of a comment)
	10090	11287	Pzp	-	A1m|A2m|AI893533|MAM	MGI:MGI:87854|Ensembl:ENSMUSG00000030359|Vega:OTTMUSG00000022212	6	6 F1-G3|6 63.02 cM	pregnancy zone protein	protein-coding	Pzp	pregnancy zone protein	O	alpha 1 macroglobulin|alpha-2-M|alpha-2-macroglobulin	20140927
Note that column headers are NOT separated by tabs in this file, while columns themselves are! 
This is annoying for parsing the file and I simply used a normal list rather than a pandas data frame for reading this file.

Format of transcript length files:
	Ensembl Gene ID	Ensembl Transcript ID	Transcript length
	ENSG00000210049	ENST00000387314	71
	ENSG00000211459	ENST00000389680	954

Format of the mouse-human orthologue file:
	Human Marker Symbol	Human Entrez Gene ID	HomoloGene ID	Mouse Marker Symbol	MGI Marker Accession ID	High-level Mammalian Phenotype ID (space-delimited)
	A1BG	  1	11167	A1bg	  MGI:2152878		
	A1CF	  29974	16363	A1cf	  MGI:1917115	  MP:0010768 MP:0005384	

'''
def create_full_geneset(outfile, ensemblFiles, entrezFiles, transcriptLengthFiles, orthologueFile, 
	ensemblPattern = re.compile(".*Ensembl:(\w+)|"), ensemblVersion = 'Ensembl Genes 77',
	entrezDataDate = '', entrezDataUrl = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/'):
	
	import numpy
	
	# Read transcript length file and store transcript lengths on ensembl gene id
	transcriptLengths = {}
	for file in transcriptLengthFiles:
		filelines = open(file).read().split('\n')
		for line in filelines[1:]:			
			cols = line.split('\t')
			if len(cols)<3 or not cols[0].startswith('ENS'): continue   # there are some LRG entries in human file
			if cols[0] not in transcriptLengths: transcriptLengths[cols[0]] = []
			transcriptLengths[cols[0]].append(int(cols[2]))
	
	# Read orthologue file. Use mouse gene symbol and human Entrez gene ids, because orthologue file does not have mouse Entrez ids.
	# Multiple matches are allowed so values of these dictionaries are lists.
	mouseSymbolsFromHumanId, humanIdsFromMouseSymbol = {}, {}
	filelines = open(orthologueFile).read().split('\n')
	for line in filelines:
		cols = line.split('\t')
		if len(cols)<4: continue
		humanId = cols[1].strip()
		mouseSymbol = cols[3].strip()
		if humanId and mouseSymbol:
			# multiple matches are possible 226 out of 16815 human ids have more than 1 mouse symbols,
			# 104 out of 17145 mouse symbols have more than 1 human ids
			if humanId not in mouseSymbolsFromHumanId: mouseSymbolsFromHumanId[humanId] = []
			if mouseSymbol not in humanIdsFromMouseSymbol: humanIdsFromMouseSymbol[mouseSymbol] = []
			mouseSymbolsFromHumanId[humanId].append(mouseSymbol)
			humanIdsFromMouseSymbol[mouseSymbol].append(humanId)

	geneInfoFromEntrezId = {}	# {'11287':{'ensemblId':'ENSMUSG00000030359', 'synonyms':'A1m|A2m|AI893533|MAM'},...}
	geneInfoFromEnsemblId = {}	# {'ENSMUSG00000030359':{'entrezIds':['11287'], 'synonyms':['A1m|A2m|AI893533|MAM'], 'symbol':'Pzp'},...}
	ensemblData = {}
	frameData = []	# will be used as the final data frame with these columns
	frameColumns = ['EnsemblId','Species','EntrezId','GeneSymbol','Synonyms','Description','MedianTranscriptLength','Orthologue']

	# Loop through each species
	for fileIndex in range(len(entrezFiles)):
		species = 'MusMusculus' if fileIndex==0 else 'HomoSapiens'

		# Read entrezFiles and store Ensembl ids and synonyms keyed on entrez id
		filelines = open(entrezFiles[fileIndex]).read().split('\n')		
				
		for line in filelines[1:]:
			cols = line.split('\t')
			if len(cols)<5: continue
			# extract ensembl id from column 6 which looks like MGI:MGI:87854|Ensembl:ENSMUSG00000030359|Vega:OTTMUSG00000022212
			m = ensemblPattern.match(cols[5])
			ensemblId = m.group(1) if m and m.group(1) else '' # actual ensembl id portion of the pattern match
			geneInfoFromEntrezId[cols[1]] = {'ensemblId':ensemblId, 'synonyms':cols[4], 'symbol':cols[2]}
		
		# "invert" geneInfoFromEntrezId so that ensemblId is the key - may produce multiple matches of entrez ids
		# For cases where there is no ensembl id found for the entrez gene, match on symbol
		geneInfoFromSymbol = {}	# {'Rnu11':{'entrezIds':['353373'], 'synonyms':['MBII-407|Mirn801|mir-801|mmu-mir-801']},...}
		for entrezId,value in geneInfoFromEntrezId.iteritems():
			if value['ensemblId']:
				if value['ensemblId'] not in geneInfoFromEnsemblId: geneInfoFromEnsemblId[value['ensemblId']] = {'entrezIds':[], 'synonyms':[]}
				geneInfoFromEnsemblId[value['ensemblId']]['entrezIds'].append(entrezId)
				geneInfoFromEnsemblId[value['ensemblId']]['synonyms'].append(value['synonyms'])
			else:
				if value['symbol'] not in geneInfoFromSymbol: geneInfoFromSymbol[value['symbol']] = {'entrezIds':[], 'synonyms':[]}
				geneInfoFromSymbol[value['symbol']]['entrezIds'].append(entrezId)
				geneInfoFromSymbol[value['symbol']]['synonyms'].append(value['synonyms'])

		# Read ensemblFiles and construct the gene info. Note that duplicate ensembl ids are possible in this data frame
		ensemblInfo = pandas.read_csv(ensemblFiles[fileIndex], sep='\t', index_col=0, dtype={'EntrezGene ID':str, 'Description':str})
		print 'species: %s, ensembl data rows: %s' % (species, ensemblInfo.shape[0])
		
		for ensemblId,row in ensemblInfo.iterrows():
			if not ensemblId.startswith('ENS'): continue  # there are some LRG_ genes in there
			
			entrezId = row['EntrezGene ID']
			geneSymbol = row['Associated Gene Name']
			
			if pandas.isnull(geneSymbol): continue
			
			if pandas.isnull(entrezId):	# substitute entrez id coming from Entrez, but only if there is a unique match
				if ensemblId in geneInfoFromEnsemblId and len(geneInfoFromEnsemblId[ensemblId]['entrezIds'])==1:
					entrezId = geneInfoFromEnsemblId[ensemblId]['entrezIds'][0]
				elif geneSymbol in geneInfoFromSymbol and len(geneInfoFromSymbol[geneSymbol]['entrezIds'])==1:
					entrezId = geneInfoFromSymbol[geneSymbol]['entrezIds'][0]
				else:
					continue
			else:	# there is entrez id here but check to see if it matches what Entrez says
				if ensemblId in geneInfoFromEnsemblId and (entrezId not in geneInfoFromEnsemblId[ensemblId]['entrezIds']):
					continue

			# synonym comes from matching entrez id; first look in geneInfoFromEnsemblId since it has already resolved multiple entrez id matches
			if ensemblId in geneInfoFromEnsemblId:
				synonyms = geneInfoFromEnsemblId[ensemblId]['synonyms'][geneInfoFromEnsemblId[ensemblId]['entrezIds'].index(entrezId)]
			elif entrezId in geneInfoFromEntrezId:
				synonyms = geneInfoFromEntrezId[entrezId]['synonyms']
			elif geneSymbol in geneInfoFromSymbol and entrezId in geneInfoFromSymbol[geneSymbol]['entrezIds']:
				synonyms = geneInfoFromSymbol[geneSymbol]['synonyms'][geneInfoFromSymbol[geneSymbol]['entrezIds'].index(entrezId)]
			else:	# no synonyms
				synonyms = ''
				
			# remove source comment from description
			description = row['Description'].split('[')[0] if pandas.notnull(row['Description']) else ''
			
			# add median transcript length
			medianTranscriptLength = numpy.median(transcriptLengths[ensemblId]) if ensemblId in transcriptLengths else None
			
			if ensemblId not in ensemblData: ensemblData[ensemblId] = []
			ensemblData[ensemblId].append([ensemblId, species, entrezId, geneSymbol, synonyms, description, medianTranscriptLength])
	
	# Merge rows with duplicate ensembl ids, create dictionaries useful for orthologue
	filteredEnsemblData = {}
	idSymbolFromId, idSymbolFromSymbol = {}, {}
	for ensemblId,_list in ensemblData.iteritems():
		listToUse = []
		if len(_list)==1:	# only one row
			listToUse = _list[0]
		else:	# multiple rows
			if len(set([item[3] for item in _list]))!=len(_list):	# different gene symbols for single ensembl id - drop
				continue
			elif len(set([item[2] for item in _list]))==1:	# all point to one entrez id - just take first match
				listToUse = _list[0]
	
		if listToUse:
			filteredEnsemblData[ensemblId] = listToUse
			idSymbolFromId[ensemblId] = '%s:%s' % (ensemblId, listToUse[3])
			idSymbolFromSymbol[listToUse[3]] = '%s:%s' % (ensemblId, listToUse[3])
	
	frameData = []
	for ensemblId,_list in filteredEnsemblData.iteritems():
		# Insert orthologue info
		orth = ''
		if _list[1]=='MusMusculus':
			humanEntrezIds = humanIdsFromMouseSymbol.get(_list[3],[])	# list of orthologous human entrez ids
			if humanEntrezIds:	# convert these into human ensembl ids, and fetch matching human symbols
				ensemblIds = [geneInfoFromEntrezId[entrezId]['ensemblId'] for entrezId in humanEntrezIds if entrezId in geneInfoFromEntrezId]
				orth = ','.join([idSymbolFromId[id] for id in ensemblIds if id in idSymbolFromId])
		elif _list[1]=='HomoSapiens':
			mouseSymbols = mouseSymbolsFromHumanId.get(_list[2],[])
			if mouseSymbols:	# fetch matching mouse ids
				orth = ','.join([idSymbolFromSymbol[symbol] for symbol in mouseSymbols if symbol in idSymbolFromSymbol])
		frameData.append(_list+[orth])

	df = pandas.DataFrame(frameData, columns=frameColumns)
	df = df.set_index('EnsemblId')
	df.fillna('')	# hdf file writing doesn't like nan
	print 'total data rows: %s (unique rows:%s)' % (df.shape[0], len(set(df.index)))

	# Write to file
	df.to_hdf(outfile, '/data')
	
	# Write description to file also
	desc = {'Ensembl Version':ensemblVersion, 'Entrez Data Source Date':entrezDataDate, 'Entrez Data Source URL':entrezDataUrl}
	pandas.Series(desc).to_hdf(outfile,'/metadata')


def create_full_geneset_20160202():
	create_full_geneset(outfile='/Users/jchoi/projects/dev/genedataset/genedataset/data/Genes.h5',
						entrezFiles = ['/Users/jchoi/projects/received/Entrez/Mus_musculus.gene_info', '/Users/jchoi/projects/received/Entrez/Homo_sapiens.gene_info'],
						ensemblFiles = ['/Users/jchoi/projects/received/Ensembl/EnsemblGenes_MusMusculus.v83.txt', '/Users/jchoi/projects/received/Ensembl/EnsemblGenes_HomoSapiens.v83.txt'],
						transcriptLengthFiles = ['/Users/jchoi/projects/received/Ensembl/TranscriptLength_MusMusculus.v83.txt', '/Users/jchoi/projects/received/Ensembl/TranscriptLength_HomoSapiens.v83.txt'],
						orthologueFile = '/Users/jchoi/projects/received/MGI/HMD_HumanPhenotype.rpt',
						ensemblVersion = 'Ensembl Genes 83',
						entrezDataDate = '2016-02-02')


# ------------------------------------------------------------
# Tests - eg. nosetests dataset.py
# ------------------------------------------------------------
def test_subset():
	gs = Geneset()
	assert gs.size>60000
	
	# subset methods
	gs = gs.subset(queryStrings=['ccr3'])
	assert set(gs.geneIds())==set(['ENSG00000183625', 'ENSMUSG00000035448'])
	assert '"EnsemblId":"ENSG00000183625"' in gs.to_json()
	
	gs = Geneset().subset(queryStrings='ccr3')
	assert set(gs.geneIds())==set(['ENSG00000183625', 'ENSMUSG00000035448'])

	gs = Geneset().subset(queryStrings='Ccr3', caseSensitive=True)
	assert gs.geneIds()==['ENSMUSG00000035448']

	gs = Geneset().subset(queryStrings='ccr1', searchColumns=["GeneSymbol"])
	assert set(gs.geneSymbols())==set(['CCR10', 'CCR12P', 'Ccr1', 'CCR1', 'Ccr10', 'Ccr1l1'])
	
	gs = Geneset().subset(queryStrings='ccr1', searchColumns=["GeneSymbol"], matchSubstring=False)
	assert set(gs.geneSymbols())==set(['Ccr1', 'CCR1'])

	gs = Geneset().subset(queryStrings='CCR1', searchColumns=["GeneSymbol"], matchSubstring=False, caseSensitive=True)
	assert set(gs.geneSymbols())==set(['CCR1'])

	gs = Geneset(name='Search: %s' % 'myb', description='') \
				.subset(queryStrings=['myb'], searchColumns=None, species=None)
	assert gs.size()==44
	
	gs = Geneset().subset(queryStrings=['ENSMUSG00000019982', 'ENSMUSG00000047591'], searchColumns=["EnsemblId"])	
	assert gs.size()==2

	gs = Geneset().subset(queryStrings=['ENSMUSG00000019982', 'ENSMUSG00000047591'], searchColumns=["GeneId"], species="MusMusculus")	
	assert gs.size()==2	
	
def test_geneSymbols():
	gs = Geneset().subset(queryStrings=['ENSMUSG00000019982', 'ENSMUSG00000047591'], searchColumns=["GeneId"])	
	assert set(gs.geneSymbols()) == set(['Mafa', 'Myb'])
	assert gs.geneSymbols(returnType="dict")['ENSMUSG00000047591'] == 'Mafa'
	
def test_species():
	gs = Geneset().subset(queryStrings=['ENSMUSG00000019982', 'ENSMUSG00000047591'], searchColumns=["GeneId"])	
	assert gs.species()=='MusMusculus'
	assert gs.species(ignoreMixed=False)==['MusMusculus']
	gs = Geneset().subset(queryStrings='ccr3')
	assert gs.species(ignoreMixed=False)==['MusMusculus','HomoSapiens']
