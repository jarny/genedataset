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
		return self._dataframe.shape[0]
		
	def geneIds(self):
		return self._dataframe.index.tolist()
		
	def geneSymbols(self):
		return self._dataframe['GeneSymbol']
		
	def species(self):
		"""
		Return species string "MusMusculus" or "HomoSapiens". 
		Assumes that all genes in this set are of the same species, so just looks at the first gene.
		"""
		return self._dataframe['Species'].iloc[0] if self.size>0 else None
		
	def subset(self, *args, **kwargs):
		"""
		Return a copy of this instance, but with a subset of genes.
		
		Parameters
		----------
		Either a list of query strings such as subset(['myb','ccr3']) 
			or keyworded arguments:
		queryStrings: a list of query strings such as ['myb','ENSMUSG00000039601']
		species: one of ['MusMusculus','HomoSapiens'] to restrict search space; may be left out;
		caseSensitive: boolean; default False
		searchColumns: a list whose members may be any of ['EnsemblId','GeneSymbol','Synonyms','Description'].
			Leaving it out will use this full list
			
		Returns
		----------
		A new Geneset instance with the subset of genes.
			
		Examples
		----------
		>>> print Geneset().subset(['ENSMUSG00000039601','ccr3']).geneSymbols()
		EnsemblId
		ENSG00000183625       CCR3
		ENSMUSG00000035448    Ccr3
		"""
		gs = copy.copy(self)
		queryStrings = kwargs['queryStrings'] if 'queryStrings' in kwargs else args[0] if args else []
		if not queryStrings:
			return gs
			
		df = gs._dataframe
		
		if 'species' in kwargs and kwargs['species'] in ['MusMusculus','HomoSapiens']:
			df = df[df['Species']==species]
			
		caseSensitive = kwargs['caseSensitive'] if 'caseSensitive' in kwargs else False
		if not caseSensitive: queryStrings = [item.lower() for item in queryStrings]
		
		searchColumns = ['EnsemblId','GeneSymbol','Synonyms','Description']
		if 'searchColumns' in kwargs and kwargs['searchColumns'] and len(set(searchColumns).intersection(set(kwargs['searchColumns'])))>0:
			searchColumns = list(set(searchColumns).intersection(set(kwargs['searchColumns'])))

		rowsToKeep = set()
		for column in searchColumns:
			if column=='EnsemblId':	# match should be done on full string match rather than partial
				for rowIndex in df.index:
					if rowIndex in rowsToKeep: continue
					geneId = rowIndex.lower() if not caseSensitive else rowIndex
					if geneId in queryStrings:
						rowsToKeep.add(rowIndex)
			else:
				for rowIndex,value in df[column].iteritems():
					if rowIndex in rowsToKeep or not value: continue
					if not caseSensitive: value = value.lower()
					for queryString in queryStrings:
						if queryString in value: 
							rowsToKeep.add(rowIndex)
							break

		gs._dataframe = df.loc[list(rowsToKeep),:]
		return gs
		
	def to_json(self):
		"""
		Return an array of dictionaries for all genes in this geneset. Example:
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

		
# ------------------------------------------------------------
# Tests - eg. nosetests dataset.py
# ------------------------------------------------------------
def test_Geneset():
	"""
	Testing of Geneset class methods.
	"""
	gs = Geneset()
	assert gs.size>60000
	gs = gs.subset(queryStrings=['ccr3'])
	assert gs.geneIds()==['ENSG00000183625', 'ENSMUSG00000035448']
	assert gs.to_json().startswith('{"ENSG00000183625":')
	
	gs = Geneset(name='Search: %s' % 'myb', description='') \
				.subset(queryStrings=['myb'], searchColumns=None, species=None)
	assert gs.size()==44
	
	gs = Geneset().subset(queryStrings=['ENSMUSG00000019982', 'ENSMUSG00000047591'], searchColumns=["EnsemblId"])	
	assert gs.size()==2