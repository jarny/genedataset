import os, sys
import pandas, numpy

from . import dataDirectory

# ------------------------------------------------------------
# Utility functions 
# ------------------------------------------------------------
def rpkm(rawCount, totalReads, medianTranscriptLength):
	"""
	Return RPKM value for a gene. Example: rpkm(1022, 34119529, 2566)
	Formula is (10^9 * C)/(N * L). See https://www.biostars.org/p/55253/
	
    Parameters
    ----------
	rawCount: unnormalised summarised count for a gene
	totalReads: total reads for all genes in a sample - ie. library size
	medianTranscriptLength: median length of all transcripts in the gene.
		
	Returns
    ----------
	float
	"""
	return numpy.power(10,9)*rawCount/totalReads/medianTranscriptLength
	
def probeGeneMap(arrayType, probeIds=[]):
	"""
	Return dictionaries of mapping between probe ids and gene ids. Uses the matching file
	from data directory of this package to work out the mapping. Each Dataset instance should
	use its own methods there for dealing with probe to gene mapping, as it is more efficient.
	This function is mainly used when creating the data attached to each Dataset instance.

    Parameters
    ----------
	arrayType: one of ['IlluminaWG6','Affymetrix']
	probeIds: list of probe ids to restrict the search. Otherwise all probe ids will be searched.
		
	Returns
    ----------
	dictionary with keys 'geneIdsFromProbeId','probeIdsFromGeneId','nonMatchingProbeIds';
		(see example below)
		
	Examples
    ----------
	>>> print probeGeneMap("IlluminaWG6", probeIds=['ILMN_1212612'])
	{'geneIdsFromProbeId': {'ILMN_1212612': ['ENSMUSG00000039601']}, 'nonMatchingProbeIds': [], 'probeIdsFromGeneId': {'ENSMUSG00000039601': ['ILMN_1212612']}}
	
	"""
	filename = "%s/ProbeIdGeneIdMapping_%s.txt" % (dataDirectory(), arrayType)
	df = pandas.read_csv(filename, sep="\t", index_col=0)
	
	if len(probeIds)>0:  # only interested in this subset
		df = df.loc[probeIds]

	# Work out a list of probe ids for each gene and vice versa
	probeIdsFromGeneId = dict([(geneId,set()) for geneId in df['geneId']])
	geneIdsFromProbeId = dict([(probeId,set()) for probeId in df.index])
	nonMatchingProbeIds = [probeId for probeId in probeIds if probeId not in df.index] if len(probeIds)>0 else []
	
	# Because of multiple matches we have to return a list for each key
	for probeId,row in df.iterrows():
		probeIdsFromGeneId[row['geneId']].add(probeId)
		geneIdsFromProbeId[probeId].add(row['geneId'])
		
	return {'probeIdsFromGeneId': dict([(key,list(probeIdsFromGeneId[key])) for key in probeIdsFromGeneId]),
			'geneIdsFromProbeId': dict([(key,list(geneIdsFromProbeId[key])) for key in geneIdsFromProbeId]),
			'nonMatchingProbeIds': nonMatchingProbeIds}

def probeIdsFromGeneId(geneId):
	"""
	Return a list of probeIds given a geneId. This function will search through all arrayTypes,
	so probeGeneMap may be more efficient if arrayType is known.
	"""
	for filename in os.listdir(dataDirectory()):
		if filename.startswith("ProbeIdGeneIdMapping"):
			df = pandas.read_csv(filename, sep="\t", index_col=0)

def createDatasetFile(destDir, **kwargs):
	"""
	Create a HDFStore file (.h5), which can be associated with a Dataset instance.

	Parameters
	----------
	destDir: Destination directory for the file being created. Existing file of the same name will
		be removed before writing.
	
	name: Usually a short abbreviated name to associate with the dataset. Will be used as file name
		for .h5 file (so "MyDataset" will create MyDataset.h5 file), hence ensure it's file name compliant.
	attributes: A dictionary describing the dataset, with following keys (all values are strings)
		fullname: A more descriptive name for the dataset.
		version: Some version string to associate with the dataset.
		description: Description of the dataset.
		platform_type: one of ['microarray','rna-seq']
		platform_details: descriptive string for the platform, eg. "Illumina WG6 version 2"
		pubmed_id: id string to pubmed if published. Leave empty or None if not applicable.
		species: one of ['MusMusculus','HomoSapiens']
	samples: A pandas DataFrame which has sample ids matching the columns of the expression matrix. Example:
		sampleId	celltype	tissue
		s01		LSK		bone marrow
		s02		T1		peripheral blood
		(sampleId must be the index of the DataFrame.)
	
	If platform_type=='microarray', we need expression matrix with probe ids, and mappings between probes and genes:
	expression: A pandas DataFrame which has sample ids as columns and probe ids as rows (index)
		probeId		s01		s02
		ILMN_1234	4.35	8.42
		ILMN_4567	8.11	5.81
	probeIdsFromGeneId: A dictionary of form {'ENSMUSG000004353':['ILMN_1234',...],...}
	geneIdsFromProbeId: A dictionary of form {'ILMN_1234':['ENSMUSG000004353',...],...}
		Both of these dictionaries allow for the possibility of multiple matches.
		
	If platform_type=='rna-seq', we need summarised counts, cpm and rpkm values:
	'counts', 'cpm', 'rpkm': 3 pandas DataFrame instances, one for each type of value.
	All DataFrames should have sample ids as columns and gene ids as row indices.
		geneId		s01		s02
		ENSG0000324	0	324
		ENSG0000132	23	1342
	
	Returns
	----------
	A Dataset instance.

	Example
	----------
	attributes = {"fullname": "Haemopedia",
				  "version": "1.0",
				  "description": "Microarray expression of murine haematopoietic cells",
				  "platform_type": "microarray",
				  "platform_details": "Illumina WG6 version 2",
				  "pubmed_id": None,
				  "species": "MusMusculus"}
	# If samples, expression, probeIdsFromGeneId, geneIdsFromProbeId have been set as above:
	createDatasetFile("/datasets", name="haemopedia", attributes=attributes, samples=samples,
		expression=expression, probeIdsFromGeneId=probeIdsFromGeneId, geneIdsFromProbeId=geneIdsFromProbeId)
	"""
	# Parse input
	if destDir[-1]=='/':
		destDir = destDir[:-1]
	name = kwargs['name']
	attributes = kwargs['attributes']
	samples = kwargs['samples']
	expression = kwargs['expression']
	probeIdsFromGeneId = kwargs['probeIdsFromGeneId']	
	geneIdsFromProbeId = kwargs['geneIdsFromProbeId']
	
	# convert a dictionary of lists into a pandas Series with repeated indices
	values, index = [], []
	for geneId,probeIds in probeIdsFromGeneId.iteritems():
		for i in range(len(probeIds)):
			index.append(geneId)
			values.append(probeIds[i])
	probeIdsFromGeneId = pandas.Series(values, index=index)
	
	values, index = [], []
	for probeId,geneIds in geneIdsFromProbeId.iteritems():
		for i in range(len(geneIds)):
			index.append(probeId)
			values.append(geneIds[i])
	geneIdsFromProbeId = pandas.Series(values, index=index)
	
	# Save all values to file
	filepath = '%s/%s.h5' % (destDir, name)
	if os.path.isfile(filepath): os.remove(filepath)
	pandas.Series(attributes).to_hdf(filepath, '/series/attributes')
	samples.to_hdf(filepath, '/dataframe/samples')
	expression.to_hdf(filepath, '/dataframe/expression')
	probeIdsFromGeneId.to_hdf(filepath, '/series/probeIdsFromGeneId')
	geneIdsFromProbeId.to_hdf(filepath, '/series/geneIdsFromProbeId')

	return Dataset(filepath)
	
# ------------------------------------------------------------
# Dataset class 
# ------------------------------------------------------------
class Dataset(object):
	"""
	An instance of Dataset is associated with its HDF file, which contains all the information about the dataset,
	including expression matrix, sample information, and various metadata associated with the dataset.

	Attributes
	----------
	filepath: full path to the HDFStore file associated with this dataset instance
	name: name of this dataset, derived by basename of filepath (eg: Dataset('/path/to/mydataset.h5').name wil be 'mydataset')
	fullname: a more descriptive name stored within the HDFStore file.
	species: one of ['MusMusculus','HomoSapiens']
	platform_type: one of ['microarray','rna-seq']
	isRnaSeqData: boolean value determined only by platform_type=='rna-seq'
		
	"""
	def __init__(self, pathToHdf):
		self.filepath = pathToHdf
		
		self._attributes = pandas.read_hdf(self.filepath, '/series/attributes')
		self.name = os.path.splitext(os.path.basename(self.filepath))[0]
		self.fullname = self._attributes['fullname']
		self.species = self._attributes['species']
		self.platform_type = self._attributes['platform_type']
		self.isRnaSeqData = self.platform_type=='rna-seq'
		
		# samples is a pandas DataFrame which describes samples and their groupings.
		# sample_id	level	name	value	colour
		# CAGRF9188-1460	0	sampleId	CAGRF9188-1460	#ffffff	
		# CAGRF9188-1333	0	sampleId	CAGRF9188-1333	#ffffff
		self._samples = pandas.read_hdf(self.filepath, '/dataframe/samples')		
			
		# expression data is held differently based on data type
		if self.isRnaSeqData:
			self._expression = {'counts': pandas.read_hdf(self.filepath, '/dataframe/counts'),
								'cpm': pandas.read_hdf(self.filepath, '/dataframe/cpm'),
								'rpkm': pandas.read_hdf(self.filepath, '/dataframe/rpkm')}
		else:
			self._expression = pandas.read_hdf(self.filepath, '/dataframe/expression')
			self._probeIdsFromGeneId = pandas.read_hdf(self.filepath, '/series/probeIdsFromGeneId')
			self._geneIdsFromProbeId = pandas.read_hdf(self.filepath, '/series/geneIdsFromProbeId')
		
	def __repr__(self):
		return "<Dataset name:%s species:%s, platform_type:%s>" % (self.name, self.species, self.platform_type)

	def attributes(self):
		"""Return a dictionary of attributes. Keys are
		['fullname','version','description','platform_type','platform_details','pubmed_id','species']
		"""
		return self._attributes.to_dict()
					
	def hdfStore(self):
		"""Return the HDF store file.
		"""
		return pandas.HDFStore(self.filepath)
		
	# ------------------------------------------------------------
	# Expression matrix related methods 
	# ------------------------------------------------------------
	def expressionMatrix(self, geneIds=None, datatype='rpkm', sampleGroupForMean=None):
		"""Return pandas DataFrame of expression values matching geneIds.
	
		Parameters
		----------
		geneIds: a list of gene ids eg: ['ENSG000003535',...] or a string 'ENSG0000003535' or None, 
			in which case the full expression matrix is returned.
		datatype: one of ['counts','rpkm','cpm'] if self.isRnaSeqData is True, ignored if False.
		sampleGroupForMean: a sample group name eg 'celltype' which will be used to return the mean
			over the sample ids.
				
		Returns
		----------
		pandas.DataFrame instance
				
		The method works by subsetting the dataframe by index for rnaseq data or by first working out the index
		to use based on self._probeIdsFromGeneId for microarray data.
		"""
		df = self._expression[datatype] if self.isRnaSeqData else self._expression
		
		if isinstance(geneIds, str) or isinstance(geneIds, unicode):	# assume single gene was specified
			geneIds = [geneIds]
				
		if self.isRnaSeqData:
			try:
				if geneIds:
					df = df.loc[set(geneIds).intersection(set(df.index))]   # otherwise all geneIds are returned with NaN values for no matches
			except KeyError:	# no geneIds found in df
				return pandas.DataFrame()
		else:
			if geneIds:
				probeIds = []
				for key,val in self._probeIdsFromGeneId.iteritems():
					if key in geneIds:
						for probeId in val:
							probeIds.append(probeId)
				df = df.loc[probeIds]  # will be empty DataFrame (with columns though) if probeIds is empty
		
		if sampleGroupForMean:
			sgi = self.sampleGroupItems(sampleGroupForMean)
			if ','.join(sgi)!=','.join(df.columns):	# it's possible for celltypes to be defined the same as sample ids, for example
				df = df.groupby(sgi, axis=1).mean()
			
		return df
		
	def expressionValues(self, geneIds=[]):
		"""
		Return a dictionary of expression values matching geneIds.
		
		Example
		----------
		expressionValues(geneIds=['ENSMUSG00000019982','ENSMUSG00000005672'])
			{'values': 
				{
					'ILMN_2683910':{'B.1':1.2, ...},
					'ILMN_2752817':{'B.1':2.2, ...}
				},
			 'featureGroups': {'Gene1':['ILMN_2683910'], ...},
			 'cpm': {}
			}
		"""
		if isinstance(geneIds, str) or isinstance(geneIds, unicode):	# assume single gene was specified
			geneIds = [geneIds]
		cpm = {}

		if self.isRnaSeqData:
			rowIds = set(self._expression['rpkm'].index).intersection(set(geneIds))
			df = self._expression['rpkm'].loc[rowIds]
			cpm = dict([(rowId,row.to_dict()) for rowId,row in self._expression['cpm'].loc[rowIds].iterrows()])
			featureGroups = dict([(index,[index]) for index in df.index])
		else:	# get matching probeIds to use as row index
			df = self._expression.loc[self.probeIdsFromGeneIds(geneIds=geneIds)]
			featureGroups = self.probeIdsFromGeneIds(geneIds=geneIds, returnType="dict")
		
		# return a dictionary
		return {'values':dict([(rowId,row.to_dict()) for rowId,row in df.iterrows()]), 'featureGroups':featureGroups, 'cpm':cpm}
		
		
	def probeIdsFromGeneIds(self, geneIds=[], returnType="list"):
		"""Return a pandas Series with geneIds as index and matching probe ids as values.
		"""
		# The following is really really slow
		#index = set(self._probeIdsFromGeneId.index).intersection(set(geneIds))
		#return self._probeIdsFromGeneId.loc[index] if len(index)>0 else pandas.Series()
		
		pfgi = self._probeIdsFromGeneId
		if returnType=="list":
			return [value for index,value in pfgi.iteritems() if index in geneIds]
		elif returnType=="dict":
			return dict([(geneId, [value for index,value in pfgi.iteritems() if index==geneId]) for geneId in geneIds])
		
	def totalReads(self):
		"""Return a dictionary that looks like {'CAGRF7126-1213':44831299, ...}
		which holds the total number of reads (basically sum of all raw reads) keyed on sample id.
		Return an empty dictionary for microarray data
		"""
		if not self.isRnaSeqData: return {}
		return self._expression['/dataframe/counts'].sum(axis=0).to_dict()
		
	def valueRange(self):
		"""Return a list of [min,max] value for the whole dataset if it's microarray data.
		"""
		df = self._expression['rpkm'] if self.isRnaSeqData else self._expression
		return [df.min().min(), df.max().max()]
			
	# ------------------------------------------------------------
	# Sample related methods 
	# ------------------------------------------------------------
	def sampleGroupNames(self):
		"""Return a list of sample group names eg: ["celltype","tissue"]
		"""
		return self._samples.columns.tolist()


# ------------------------------------------------------------
# Tests - eg. nosetests dataset.py
# ------------------------------------------------------------
def test_utilityFunctions():
	"""Testing of utility functions.
	"""
	# rpkm
	assert rpkm(1022, 34119529, 2566)==11
	
	# probeGeneMap
	pgm = probeGeneMap("IlluminaWG6", probeIds=['ILMN_1212612','ILMN_1213657'])
	assert pgm['probeIdsFromGeneId']['ENSMUSG00000039601']==['ILMN_1212612']
	assert len(pgm['geneIdsFromProbeId']['ILMN_1213657'])==15

def test_Dataset():
	"""Test Dataset class methods.
	"""
	ds = Dataset("%s/testdata.h5" % dataDirectory())
	
	# metadata test
	assert ds.platform_type=="microarray"
	assert not ds.isRnaSeqData
	assert ds.species=="MusMusculus"
	
	# probeIdsFromGeneIds
	assert ds.probeIdsFromGeneIds(geneIds=['gene1'])==['probe1']
	assert ds.probeIdsFromGeneIds()==[]
	assert ds.probeIdsFromGeneIds(geneIds=['gene2'], returnType="dict")=={'gene2':['probe2','probe3']}
	
	# expressionValues
	ev = ds.expressionValues(geneIds=['gene2'])
	assert ev['values']['probe2']['s01']==5.54
	assert ev['featureGroups']['gene2']==['probe2','probe3']

	# samples
	assert ds.sampleGroupNames()==['celltype','tissue']
	