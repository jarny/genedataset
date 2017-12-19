import os, sys
import pandas, numpy
import six

from . import dataDirectory

# ------------------------------------------------------------
# Utility functions 
# ------------------------------------------------------------
def createDatasetFile(destDir, **kwargs):
	"""
	Create a HDFStore file (.h5), which can be associated with a Dataset instance.
	This file contains multiple pandas DataFrame objects as well as 1 pandas Series object.
	The pandas Series object holds dataset metadata, such as name and description.
	The pandas DataFrame objects hold sample table and one or more expression matrices, as it
	is often useful to have raw counts as well as cpm matrices, for example.
		
	The expression DataFrame object should have sample ids as columns and feature ids as row indices.
	===================== ====== ======
	geneId                s01    s02 
	===================== ====== ======
	ENSMUSG000004353      0      324
	ENSMUSG00000039601    23     1342
	===================== ====== ======
	Often a dataset may contain multiple expression data frames and you can specify these by providing
	a list of data frames and specifying the keys used for these within 'expression_data_keys' key of
	the attributes dictionary. See example below.

	While the samples DataFrame object should have sample ids as rows and columns are sample groupings.
	See example below.

	Parameters:
		destDir: Destination directory for the file being created. Existing file of the same name will
			be removed before writing.
			
		attributes: dictionary describing the dataset, with following keys (all values are strings)
			name: usually the same as name attribute above, but make it different if different filename
				is desired to dataset name. Example: name='MyDataset', attributes={'name':'mydata',...}
				will create MyDataset.h5 file but the dataset name will be 'mydata'.
			fullname: A more descriptive name for the dataset.
			version: Some version string to associate with the dataset.
			description: Description of the dataset.
			expression_data_keys: list of keys describing the types of expression matrices, eg ['counts','cpm']
			pubmed_id: id string to pubmed if published. Leave empty or None if not applicable.
			species: one of ['MusMusculus','HomoSapiens']
			
		samples: pandas.DataFrame where sample ids must be the index matching the columns of the expression matrices.
			sampleId	celltype	tissue
			s01		LSK		bone marrow
			s02		T1		peripheral blood
			(sampleId must be the index of the DataFrame.)
	
		expressions: list of pandas.DataFrame's, where each DataFrame has sample ids as columns and feature ids as index
		Note that order of this list must match the order given by expression_data_keys.
			featureId		s01		s02
			ENSG000003535	4.35	8.42
			ENSG000004215	8.11	5.81

	Returns:
		A Dataset instance.

	Example:
		import pandas
		from genedataset import dataset
		attributes = {"name": "testdata",
					  "fullname": "Test Dataset",
					  "version": "1.0",
					  "description": "This dataset comes with the package for testing purposes.",
					  "expression_data_keys": ["counts","cpm"],
					  "pubmed_id": None,
					  "species": "MusMusculus"}
		samples = pandas.DataFrame([['B1', 'BM'], ['B1', 'BM'], ['B2', 'BM'], ['B2', 'BM']],
									index=['s01','s02','s03','s04'], columns=['celltype','tissue'])
		samples.index.name = "sampleId"
		counts = pandas.DataFrame([[35, 44, 21, 101], [50, 0, 14, 62], [0, 0, 39, 73]],
									index=['gene1', 'gene2', 'gene3'], columns=['s01', 's02', 's03', 's04'])
		counts.index.name = "geneId"
		cpm = pandas.DataFrame([[3.45, 4.65, 2.65, 8.23], [5.54, 0.0, 1.43, 6.43], [0.0, 0.0, 4.34, 5.44]],
								 index=['gene1', 'gene2', 'gene3'], columns=['s01', 's02', 's03', 's04'])
		cpm.index.name = "geneId"
		dataset.createDatasetFile("/datasets", attributes=attributes, samples=samples, expressions=[counts,cpm])
	"""
	# Parse input
	if destDir[-1]=='/':
		destDir = destDir[:-1]
	attributes = kwargs['attributes']
	samples = kwargs['samples']
	expressions = kwargs['expressions']

	# Check data
	error = None
	if not isinstance(samples, pandas.DataFrame):
		error = "samples table is not a pandas DataFrame object."
	elif not isinstance(expressions, list):
		error = "expressions is not a list object."
	elif not all(isinstance(item, pandas.DataFrame) for item in expressions):
		error = "one or more of the objects in expressions is not a pandas DataFrame object."
	elif not attributes:
		error = "attributes argument is missing or empty."
	elif "name" not in attributes:
		error = "'name' must be a key in attributes dictionary."
	elif "expression_data_keys" not in attributes:
		error = "'expression_data_keys' must be a key in attributes dictionary."
	elif not attributes['expression_data_keys']:
		error = "expression_data_keys in attributes dictionary must name the keys used for expressions."
	elif any(len([item for item in df.columns if item not in samples.index])>0 for df in expressions):
		error = "some sample ids in expression matrix missing in sample matrix"		

	if error is not None:
		print(error)
		return

	# Save all values to file
	filepath = os.path.join(destDir, '%s.%s.h5' % (attributes['name'], attributes['version']))
	if os.path.isfile(filepath): os.remove(filepath)

	# keys for expressions
	pandas.Series(attributes).to_hdf(filepath, '/series/attributes')
	samples.to_hdf(filepath, '/dataframe/samples')
	for i,key in enumerate(attributes['expression_data_keys']):
		expressions[i].to_hdf(filepath, '/dataframe/expression/%s' % key)

	# Return the Dataset instance
	return Dataset(filepath)
	
# ------------------------------------------------------------
# Dataset class 
# ------------------------------------------------------------
class Dataset(object):
	"""
	An instance of Dataset is associated with its HDF file, which contains all the information about the dataset,
	including expression matrices, sample information, and various metadata associated with the dataset.

	Attributes:
		name: (string) name of this dataset that can act as an identifier within an application
		filepath: (string) full path to the HDF file associated with this dataset instance
		attributes: (dict) complete metadata dictionary with keys 'description', 'expression_data_keys', 'fullname', 
			'name', 'pubmed_id', 'species', 'version'
		expression_data_keys: (list) keys used to retrieve different expression matrices, 
			built from self.attributes.expression_data_keys
		
	"""
	def __init__(self, pathToHdf):
		self.filepath = pathToHdf
		
		self.attributes = pandas.read_hdf(self.filepath, '/series/attributes').to_dict()
		self.name = self.attributes['name']
		self.expression_data_keys = self.attributes['expression_data_keys']#.split(',')
		
		# samples is a pandas DataFrame which describes samples and their groupings.
		# sample_id	level	name	value	colour
		# CAGRF9188-1460	0	sampleId	CAGRF9188-1460	#ffffff	
		# CAGRF9188-1333	0	sampleId	CAGRF9188-1333	#ffffff
		self.samples = pandas.read_hdf(self.filepath, '/dataframe/samples')

		# There may be multiple expression matrices, and we use attributes.expression_data_keys to find them
		self.expressions = dict([(key, pandas.read_hdf(self.filepath, '/dataframe/expression/%s' % key)) \
								for key in self.attributes['expression_data_keys']])
		
	def __repr__(self):
		return "<Dataset name:%s, %s samples>" % (self.name, len(self.samples))
					
	def hdfStore(self):
		"""Return the HDF store file.
		"""
		return pandas.HDFStore(self.filepath)
		
	# ------------------------------------------------------------
	# Expression matrix related methods 
	# ------------------------------------------------------------
	def expressionMatrix(self, expression_data_key=None, featureIds=None, sampleGroupForMean=None):
		"""Return pandas DataFrame of expression values matching featureIds.
	
		Parameters:
			expression_data_key: specifies which expression matrix to fetch. If None, it will fetch the first
				one specified by self.expression_data_keys
			featureIds: a list of row ids eg: ['ENSG000003535',...] or a string 'ENSG0000003535' or None, 
				in which case the full expression matrix is returned.
			sampleGroupForMean: a sample group name eg 'celltype' which will be used to return the mean
				over the sample ids.
				
		Returns:
			pandas.DataFrame instance	
		"""
		df = self.expressions[expression_data_key if expression_data_key in self.expression_data_keys else \
								self.expression_data_keys[0]]
		
		if isinstance(featureIds, six.string_types):	# assume single feature was specified
			featureIds = [featureIds]
		
		# subset rows if required
		if featureIds is not None and len(featureIds)>0:
			index = set(featureIds).intersection(set(df.index))
		else:
			index = df.index
		
		if len(index)>0:
			df = df.loc[index]
		else:
			return pandas.DataFrame()
					
		if sampleGroupForMean:
			sgi = self.sampleGroupItems(sampleGroup=sampleGroupForMean, duplicates=True)
			if ','.join([item if pandas.notnull(item) else '' for item in sgi])!=','.join(df.columns):	
				# it's possible for celltypes to be defined the same as sample ids, for example
				df = df.groupby(sgi, axis=1).mean()
			
		return df		
				
	def totalReads(self, expression_data_key=None):
		"""Return a dictionary that looks like {'CAGRF7126-1213':44831299, ...}
		which holds the total number of reads (basically sum of all raw reads) keyed on sample id.
		"""
		return self.expressionMatrix(expression_data_key=expression_data_key).sum(axis=0).to_dict()
		
	def valueRange(self, expression_data_key=None):
		"""Return a list of [min,max] value for the whole expression matrix.
		"""
		return [self.expressionMatrix(expression_data_key=expression_data_key).min().min(), self.expression.max().max()]
			
	def correlation(self, featureId, expression_data_key=None):
		"""
		Return a dictionary of correlation scores for each feature id in the dataset, eg: {'ENSMUSG0000034355':0.343, ...}
		Each value is the correlation between featureId and each feature in the expression matrix.
		featureId must be one of the row indices of the dataset. Otherwise None is returned.
		In this case only the keys contained in this dictionary are returned.
		"""
		from scipy.stats.stats import pearsonr
		
		# Use log2
		df = numpy.log2(self.expressionMatrix(expression_data_key=expression_data_key) + 1)
		if featureId not in df.index: return None
				
		# values for featureId
		values = df.loc[featureId]
			
		# loop through all features and calculate correlation
		score = {}
		for featureId,row in df.iterrows():
			if row.max()==0: continue
			corr = pearsonr(values, row)[0]   # only need correlation and not the p value
			if pandas.notnull(corr):   # null can happen if row is all zeros (or all the same values => zero covariance)
				score[featureId] = corr

		return score				
		
	# ------------------------------------------------------------
	# Sample related methods 
	# ------------------------------------------------------------
	def sampleTable(self):
		"""Return a pandas data frame of the samples table.
		"""
		return self.samples
		
	def sampleGroups(self):
		"""Return a list of sample group names eg: ["celltype","tissue"]
		"""
		return self.samples.columns.tolist()

	def sampleGroupItems(self, sampleGroup=None, groupBy=None, duplicates=False):
		"""
		Return sample group items belonging to sampleGroup, eg: ["B1","B2"].
		
		Parameters:
			sampleGroup: name of sample group, eg: 'celltype'
			groupBy: name of another sample group for grouping, eg: 'cell_lineage'
			duplicates: boolean to return a list of unique values in the list avoiding duplicates if False; 
				if True, it specifies a list of sample group items in the same order/position as columns of 
				expression matrix is requested; ignored if groupBy is specified.
		
		Returns:
			a list if groupBy is None, eg: ['B1','B2',...]. If duplicates is True, the list
				returned is the same length as dataset's columns in its expression matrix, and in the same
				corresponding position.
			a dictionary of list if groupBy is specified, eg: {'Stem Cell':['LSK','STHSC',...], ...}
				groupBy sorts the flat list which would have been returned without groupBy into
				appropriate groups.
			
		Note that this method does not make assumptions about the integrity of the data returned for 
		groupBy specification. So it's possible to return {'Stem Cell':['LSK','STHSC'], 'B Cells':['LSK','B1']},
		if there is a sample id which has been assigned to ('LSK','Stem Cell') and another to ('LSK','B Cells') by mistake.
		"""
		df = self.samples
	
		if sampleGroup in df.columns and groupBy in df.columns: # group each item by sample ids, then substitute items from sampleGroup
			sampleIdsFromGroupBy = dict([(item, df[df[groupBy]==item].index.tolist()) for item in set(df[groupBy])])
			# {'Stem Cell':['sample1','sample2',...], ... }
		
			# substitute items from sampleGroup for each sample id
			dictToReturn = {}
			for sampleGroupItem in sampleIdsFromGroupBy.keys():
				# this is the set of matching sample group items with duplicates removed, eg: ['LSK','CMP',...]
				dictToReturn[sampleGroupItem] = sorted(set([df.at[sampleId,sampleGroup] for sampleId in sampleIdsFromGroupBy[sampleGroupItem]]))                
			return dictToReturn
		
		elif sampleGroup in df.columns:
			if duplicates:
				sampleIds = self.expressionMatrix().columns
				return df.loc[sampleIds][sampleGroup].tolist()
			else:
				return sorted(set(df[sampleGroup]))

		else:
			return []
			
	def sampleIds(self, sampleGroup=None, sampleGroupItem=None):
		"""Return a list of sample ids, eg: ["sample1","sample2"].
		If either of sampleGroup or sampleGroupItem is None, returns all sampleIds.
		"""
		df = self.samples
		if sampleGroup and sampleGroupItem:
			return df[df[sampleGroup]==sampleGroupItem].index.tolist()
		else:
			return df.index.tolist()
	
	def sampleIdsFromSampleGroups(self):
		"""Return a dictionary of sample ids as lists, keyed on sample group items,
		eg: {'celltype': {'B':['s1','s2',...], ...}, ... }
		"""
		dict = {}
		df = self.samples
		for sampleGroup in df.columns:
			dict[sampleGroup] = {}
			for sampleId,value in df[sampleGroup].iteritems():
				if pandas.isnull(value): continue
				if value not in dict[sampleGroup]: dict[sampleGroup][value] = []
				dict[sampleGroup][value].append(sampleId)
		return dict
	
		
# ------------------------------------------------------------
# Tests - eg. nosetests dataset.py
# ------------------------------------------------------------
def test_createDatasetFile():
	attributes = {"name":"test",
				  "fullname": "Test Dataset",
				  "version": "1.0",
				  "description": "Created as part of daset test function",
				  "expression_data_keys": ['expression'],
				  "pubmed_id": None,
				  "species": "MusMusculus"}
	samples = pandas.DataFrame([['B1','B Cell Lineage'],['B2','B Cell Lineage']], index=['sample1','sample2'], columns=['celltype','cell_lineage'])
	expression = pandas.DataFrame([[3,0],[5,2],[0,4]], index=['gene1','gene2','gene3'], columns=['sample1','sample2'])
	createDatasetFile("/tmp", attributes=attributes, samples=samples, expressions=[expression])

def test_attributes():
	ds = Dataset("/tmp/test.1.0.h5")
	assert ds.attributes['species']=="MusMusculus"
	assert ds.name=="test"
	assert ds.attributes['pubmed_id'] is None

def test_expressionMatrix():
	ds = Dataset("/tmp/test.1.0.h5")
	assert ds.expressionMatrix()['sample1'].tolist()==[3,5,0]
	assert ds.expressionMatrix(sampleGroupForMean='celltype').at['gene2','B2']==2

def test_samples():
	ds = Dataset("%s/testdata.1.0.h5" % dataDirectory())	
	assert ds.sampleGroups()==['celltype','tissue']
	assert ds.sampleGroupItems(sampleGroup="celltype")==['B1','B2']
	assert ds.sampleGroupItems(sampleGroup='celltype', groupBy='tissue')=={'BM': ['B1', 'B2']}
	assert ds.sampleGroupItems(sampleGroup='celltype', duplicates=True)==['B1','B1','B2','B2']
	assert ds.sampleIds()==['s01','s02','s03','s04']
	assert ds.sampleIds(sampleGroup='celltype', sampleGroupItem='B1')==['s01','s02']
	assert ds.sampleIdsFromSampleGroups()['tissue']['BM']==['s01', 's02', 's03', 's04']

