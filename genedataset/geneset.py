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
