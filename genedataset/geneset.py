"""
Geneset class is a list of information about genes. It comes with annotated gene information sourced from 
Ensembl and Entrez (NCBI) for mouse and human genes, to make it easier to work with python scripts in bioinformatics.

The main object which contains the information is a pandas data frame, with each gene per row. This data frame has columns:
['EnsemblId','Species','EntrezId','GeneSymbol','Synonyms','Description','TranscriptLengths','Orthologue']
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
		self._datafile = '{}/Genes.h5'.format(dataDirectory())
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
		
	def subsetFromGeneIds(self, geneIds):
		"""
		Return a copy of this instance, but with a subset of genes based on geneIds. Gene ids not found in this Geneset
		will be ignored.
		
		Parameters:			
			geneIds: list of gene ids such as ['ENSMUSG00000039601', ...]
			
		Returns:
			A new Geneset instance with the subset of genes.			
		"""
		gs = copy.copy(self)
		gs._dataframe = gs._dataframe.loc[[item for item in geneIds if item in gs._dataframe.index]]
		return gs
		
	def subsetFromGeneSymbols(self, geneSymbols):
		"""
		Return a copy of this instance, but with a subset of genes based on gene symbols. This method is for exact matching -
		use subset() for more loose matching.
		
		Parameters:			
			geneSymbols: list of gene symbols such as ['Gata1', ...]
			
		Returns:
			A new Geneset instance with the subset of genes.			
		"""
		gs = copy.copy(self)
		gs._dataframe = gs._dataframe[gs._dataframe['GeneSymbol'].isin(geneSymbols)]
		return gs
		
	def subset(self, *args, **kwargs):
		"""
		Return a copy of this instance, but with a subset of genes based on query strings. Slower than subsetFromGeneIds
		or subsetFromGeneSymbols but more generic.
		
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
		if isinstance(queryStrings, str) or isinstance(queryStrings, bytes):	# assume queryString was specified as a string
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
Create the Genes.h5 file from source files (Updated 2017-05-26). This requires downloading various data files in text format
and running this function, which will create the .h5 file containing a pandas dataframe. An example:
	create_full_geneset(outfile='/Users/jchoi/projects/dev/genedataset/python3/genedataset/data/Genes.h5',
						entrezFiles = ['/Users/jchoi/projects/received/Entrez/Mus_musculus.gene_info', 
									   '/Users/jchoi/projects/received/Entrez/Homo_sapiens.gene_info'],
						ensemblFiles = ['/Users/jchoi/projects/received/Ensembl/EnsemblGenes_MusMusculus.v88.GRCh38.txt',
										'/Users/jchoi/projects/received/Ensembl/EnsemblGenes_HomoSapiens.v88.GRCh38.txt'],
						transcriptLengthFiles = ['/Users/jchoi/projects/received/Ensembl/TranscriptLength_MusMusculus.v88.GRCh38.txt',
												 '/Users/jchoi/projects/received/Ensembl/TranscriptLength_HomoSapiens.v88.GRCh38.txt'],
						orthologueFile = '/Users/jchoi/projects/received/Ensembl/MouseHumanOrthologues.v88.GRCh38.txt',
						ensemblVersion = 'Ensembl Genes 88, GRCh38',
						entrezDataDate = '2017-05-26')
Note that mouse file should always come before the human file within each pair.
Beware of Entrez gene id being parsed as integer by default - this can cause problems in various places.

The source data files come from:
	Ensembl files downloaded from http://ensembl.org/biomart
	Entrez files downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/

Format of the Ensembl gene files: (Be careful when selecting attributes to download from biomart - columns are in the same order as selection.)
	Gene stable ID	NCBI gene ID	Gene name	Gene description
	ENSMUSG00000064336		mt-Tf	mitochondrially encoded tRNA phenylalanine [Source:MGI Symbol;Acc:MGI:102487]
	ENSMUSG00000064337		mt-Rnr1	mitochondrially encoded 12S rRNA [Source:MGI Symbol;Acc:MGI:102493]

Format of Ensembl transcript length files:
	Gene stable ID	Transcript stable ID	Transcript length (including UTRs and CDS)
	ENSMUSG00000064372	ENSMUST00000082423	67
	ENSMUSG00000064371	ENSMUST00000082422	67

Format of the mouse-human orthologue file:
	Gene stable ID	Gene name	Human gene stable ID	Human gene name
	ENSMUSG00000064372	mt-Tp		
	ENSMUSG00000064370	mt-Cytb	ENSG00000198727	MT-CYB

Format of the Entrez files:
	#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date
	9606	1	A1BG	-	A1B|ABG|GAB|HYST2477	MIM:138670|HGNC:HGNC:5|Ensembl:ENSG00000121410|Vega:OTTHUMG00000183507	19	19q13.43	alpha-1-B glycoprotein	protein-coding	A1BG	alpha-1-B glycoprotein	O	alpha-1B-glycoprotein|HEL-S-163pA|epididymis secretory sperm binding protein Li 163pA	20170508
	9606	2	A2M	-	A2MD|CPAMD5|FWP007|S863-7	MIM:103950|HGNC:HGNC:7|Ensembl:ENSG00000175899|Vega:OTTHUMG00000150267	12	12p13.31	alpha-2-macroglobulin	protein-coding	A2M	alpha-2-macroglobulin	O	alpha-2-macroglobulin|C3 and PZP-like alpha-2-macroglobulin domain-containing protein 5|alpha-2-M	20170513

'''
def create_full_geneset(outfile, ensemblFiles, entrezFiles, transcriptLengthFiles, orthologueFile, 
	ensemblPattern = re.compile(".*Ensembl:(\w+)|"), ensemblVersion = 'Ensembl Genes 88',
	entrezDataDate = '', entrezDataUrl = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/'):
	
	import pandas, numpy
			
	# Read orthologue file. Multiple matches are allowed so values of these dictionaries are lists.
	print("Parsing orthologue file:", orthologueFile)
	orth = pandas.read_csv(orthologueFile, sep='\t').fillna('')
	
	mouseIdsFromHumanId, humanIdsFromMouseId = {}, {}
	for index,row in orth.iterrows():
		humanId = row['Human gene stable ID'].strip()
		mouseId = row['Gene stable ID'].strip()
		if humanId and mouseId:	# multiple matches are possible
			if humanId not in mouseIdsFromHumanId: mouseIdsFromHumanId[humanId] = []
			if mouseId not in humanIdsFromMouseId: humanIdsFromMouseId[mouseId] = []
			mouseIdsFromHumanId[humanId].append(mouseId)
			humanIdsFromMouseId[mouseId].append(humanId)
	
	# These will store the info at different stages
	geneInfoFromEntrezId = {}	# {'11287':{'ensemblId':'ENSMUSG00000030359', 'synonyms':'A1m|A2m|AI893533|MAM'},...}
	geneInfoFromEnsemblId = {}	# {'ENSMUSG00000030359':{'entrezIds':['11287'], 'synonyms':['A1m|A2m|AI893533|MAM'], 'symbol':'Pzp'},...}
	ensemblData = {}
	frameData = []	# will be used as the final data frame with these columns
	frameColumns = ['EnsemblId','Species','EntrezId','GeneSymbol','Synonyms','Description','TranscriptLengths','Orthologue']

	# Loop through each species
	for fileIndex,species in enumerate(['MusMusculus','HomoSapiens']):

		# Read transcript length file and store transcript lengths on ensembl gene id
		print("Parsing transcript length file:", transcriptLengthFiles[fileIndex])	
		df = pandas.concat([pandas.read_csv(filepath, sep='\t') for filepath in transcriptLengthFiles], ignore_index=True)
		transcriptLengths = {}
		for index,row in df.iterrows():
			if row['Gene stable ID'] not in transcriptLengths:
				transcriptLengths[row['Gene stable ID']] = []
			transcriptLengths[row['Gene stable ID']].append(row[2])
		print("Transcript lengths of ENSMUSG00000000031 is", transcriptLengths['ENSMUSG00000000031'])
	
		# Read entrez file and store Ensembl ids and synonyms keyed on entrez id
		entrez = pandas.read_csv(entrezFiles[fileIndex], sep='\t', index_col=1, dtype={'GeneID':str})
		entrez.index = entrez.index.astype(str)  # note that dtype={'GeneID':str} within read_csv doesn't work!
		print("\nParsing Entrez file:", entrezFiles[fileIndex])
		print("Number of rows:", len(entrez), "Are there duplicate gene ids?", len(entrez.index)!=len(set(entrez.index)))
		
		for index,row in entrez.iterrows():
			# extract ensembl id from dbXrefs columns which looks like MGI:MGI:87854|Ensembl:ENSMUSG00000030359|Vega:OTTMUSG00000022212
			m = ensemblPattern.match(row['dbXrefs'])
			ensemblId = m.group(1) if m and m.group(1) else '' # actual ensembl id portion of the pattern match
			geneInfoFromEntrezId[index] = {'ensemblId':ensemblId, 'synonyms':row['Synonyms'], 'symbol':row['Symbol']}
			
		# "invert" geneInfoFromEntrezId so that ensemblId is the key - may produce multiple matches of entrez ids
		# For cases where there is no ensembl id found for the entrez gene, match on symbol
		geneInfoFromSymbol = {}	# {'Rnu11':{'entrezIds':['353373'], 'synonyms':['MBII-407|Mirn801|mir-801|mmu-mir-801']},...}
		for entrezId,value in geneInfoFromEntrezId.items():
			if value['ensemblId']:
				if value['ensemblId'] not in geneInfoFromEnsemblId: geneInfoFromEnsemblId[value['ensemblId']] = {'entrezIds':[], 'synonyms':[]}
				geneInfoFromEnsemblId[value['ensemblId']]['entrezIds'].append(entrezId)
				geneInfoFromEnsemblId[value['ensemblId']]['synonyms'].append(value['synonyms'])
			else:
				if value['symbol'] not in geneInfoFromSymbol: geneInfoFromSymbol[value['symbol']] = {'entrezIds':[], 'synonyms':[]}
				geneInfoFromSymbol[value['symbol']]['entrezIds'].append(entrezId)
				geneInfoFromSymbol[value['symbol']]['synonyms'].append(value['synonyms'])
		
		# Read ensembl file and construct the gene info. Note that duplicate ensembl ids are possible in this data frame
		print("\nParsing Ensembl gene file:", ensemblFiles[fileIndex])
		ensembl = pandas.read_csv(ensemblFiles[fileIndex], sep='\t', index_col=0, dtype={'NCBI gene ID':str, 'Gene description':str})
		print("Number of rows:", len(ensembl), "Are there duplicate gene ids?", len(ensembl.index)!=len(set(ensembl.index)))
		
		missingGeneSymbol, missingEntrezId, nonMatchingEntrezId = 0, 0, 0
		for ensemblId,row in ensembl.iterrows():
			if not ensemblId.startswith('ENS'): continue  # there may be some LRG_ genes in there

			entrezId = row['NCBI gene ID']
			geneSymbol = row['Gene name']
			
			if pandas.isnull(geneSymbol):
				missingGeneSymbol
				continue
			
			if pandas.isnull(entrezId):	
				# substitute entrez id coming from Entrez, but what if there are multiple matches?
				# I used to insist on keeping the gene only if there was a unique match here, but now just take the first one
				# - note this is quite arbitrary!
				if ensemblId in geneInfoFromEnsemblId:
					entrezId = geneInfoFromEnsemblId[ensemblId]['entrezIds'][0]
				elif geneSymbol in geneInfoFromSymbol:
					entrezId = geneInfoFromSymbol[geneSymbol]['entrezIds'][0]
				else:
					# These genes with missing entrez ids are mostly pseudogenes (eg 'Gm26540'), which we will leave out here
					missingEntrezId += 1
					continue
					
			else:	# there is entrez id here but check to see if it matches what Entrez says
				if ensemblId in geneInfoFromEnsemblId and (entrezId not in geneInfoFromEnsemblId[ensemblId]['entrezIds']):
					nonMatchingEntrezId += 1
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
			synonyms = synonyms.strip()
			if synonyms=='-': synonyms = ''
				
			# remove source comment from description
			description = row['Gene description'].split('[')[0] if pandas.notnull(row['Gene description']) else ''
			
			# add transcript lengths as a numpy array
			transLengths = numpy.asarray(transcriptLengths.get(ensemblId,[]))
			
			if ensemblId not in ensemblData: ensemblData[ensemblId] = []
			ensemblData[ensemblId].append([ensemblId, species, entrezId, geneSymbol.strip(), synonyms, description.strip(), transLengths])
			
		print("Length of ensemblData:", len(ensemblData))
		print("Rows dropped due to missing gene symbol:", missingGeneSymbol)
		print("Rows dropped due to missing Entrez id:", missingEntrezId)
		print("Rows dropped due to inconsistent Entrez id between Ensembl and Entrez:", nonMatchingEntrezId)
		
	# Merge rows with duplicate ensembl ids; also create a dictionary which will be used for orthologue info later
	print("\nMerging rows with duplicate Ensembl ids")
	filteredEnsemblData = {}
	idAndSymbolFromId = {}
	for ensemblId,_list in ensemblData.items():
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
			idAndSymbolFromId[ensemblId] = '{}:{}'.format(ensemblId, listToUse[3])
	
	frameData = []
	for ensemblId,_list in filteredEnsemblData.items():
		# Insert orthologue info
		orth = ''
		if _list[1]=='MusMusculus':
			humanIds = humanIdsFromMouseId.get(ensemblId,[])	# list of orthologous human ensembl ids
			if humanIds:	# fetch matching human info
				orth = ','.join([idAndSymbolFromId[id] for id in humanIds if id in idAndSymbolFromId])
		elif _list[1]=='HomoSapiens':
			mouseIds = mouseIdsFromHumanId.get(ensemblId,[])
			if mouseIds:	# fetch matching mouse info
				orth = ','.join([idAndSymbolFromId[id] for id in mouseIds if id in idAndSymbolFromId])
		frameData.append(_list+[orth])

	df = pandas.DataFrame(frameData, columns=frameColumns).set_index('EnsemblId').fillna('')	# hdf file writing doesn't like nan
	print('\nCreating the final data frame - total data rows: {} (unique rows:{})'.format(df.shape[0], len(set(df.index))))
	
	# Write to file
	df.to_hdf(outfile, '/data')
	
	# Write description to file also
	desc = {'Ensembl Version':ensemblVersion, 'Entrez Data Source Date':entrezDataDate, 'Entrez Data Source URL':entrezDataUrl}
	pandas.Series(desc).to_hdf(outfile,'/metadata')


def create_full_geneset_20170526():
	create_full_geneset(outfile='/Users/jchoi/projects/dev/genedataset/python3/genedataset/data/Genes.h5',
						entrezFiles = ['/Users/jchoi/projects/received/Entrez/Mus_musculus.gene_info', 
									   '/Users/jchoi/projects/received/Entrez/Homo_sapiens.gene_info'],
						ensemblFiles = ['/Users/jchoi/projects/received/Ensembl/EnsemblGenes_MusMusculus.v88.GRCh38.txt',
										'/Users/jchoi/projects/received/Ensembl/EnsemblGenes_HomoSapiens.v88.GRCh38.txt'],
						transcriptLengthFiles = ['/Users/jchoi/projects/received/Ensembl/TranscriptLength_MusMusculus.v88.GRCh38.txt',
												 '/Users/jchoi/projects/received/Ensembl/TranscriptLength_HomoSapiens.v88.GRCh38.txt'],
						orthologueFile = '/Users/jchoi/projects/received/Ensembl/MouseHumanOrthologues.v88.GRCh38.txt',
						ensemblVersion = 'Ensembl Genes 88, GRCh38',
						entrezDataDate = '2017-05-26')


# ------------------------------------------------------------
# Tests - eg. nosetests geneset.py
# ------------------------------------------------------------
def test_Geneset():
	gs = Geneset()
	assert gs.size()>60000
	assert Geneset().subsetFromGeneIds(['ENSMUSG00000035448']).geneSymbols()[0]=='Ccr3'
	assert Geneset().subsetFromGeneSymbols(['Ccr3']).geneIds()[0]=='ENSMUSG00000035448'
	
def test_subset():
	gs = Geneset()
	
	# subset methods
	gs = gs.subset(queryStrings='ccr3')
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
	assert set(Geneset().species(ignoreMixed=False))==set(['HomoSapiens','MusMusculus'])
	gs = Geneset().subsetFromGeneIds(['ENSMUSG00000019982', 'ENSMUSG00000047591'])
	assert gs.species()=='MusMusculus'
	assert gs.species(ignoreMixed=False)==['MusMusculus']
	gs = Geneset().subset(queryStrings='ccr3')
	print(gs.species())