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
	
	Parameters:
		rawCount: unnormalised summarised count for a gene
		totalReads: total reads for all genes in a sample - ie. library size
		medianTranscriptLength: median length of all transcripts in the gene.
		
	Returns:
		float
	"""
	return numpy.power(10,9)*rawCount/totalReads/medianTranscriptLength
	
def probeGeneMap(arrayType, probeIds=[]):
	"""
	Return dictionaries of mapping between probe ids and gene ids. Uses the matching file
	from data directory of this package to work out the mapping. Each Dataset instance should
	use its own methods there for dealing with probe to gene mapping, as it is more efficient.
	This function is mainly used when creating the data attached to each Dataset instance.

	Parameters:
		arrayType: {'IlluminaWG6','Affymetrix'}
		probeIds: list of probe ids to restrict the search. Otherwise all probe ids will be searched.
		
	Returns:
		dictionary with keys 'geneIdsFromProbeId','probeIdsFromGeneId','nonMatchingProbeIds' (see example below)
		
	Examples:
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
		if row['geneId']: probeIdsFromGeneId[row['geneId']].add(probeId)
		if probeId: geneIdsFromProbeId[probeId].add(row['geneId'])
		
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
	If platform_type=='microarray', we need expression matrix with probe ids, and mappings between probes and genes,
	so supply expression, probeIdsFromGeneId and geneIdsFromProbeId.
	If platform_type=='rna-seq', we need summarised counts, cpm and rpkm values.
	
	All DataFrames should have sample ids as columns and gene ids as row indices.
	===================== ====== ======
	geneId                s01    s02 
	===================== ====== ======
	ENSMUSG000004353      0      324
	ENSMUSG00000039601    23     1342
	===================== ====== ======

	Parameters:
		destDir: Destination directory for the file being created. Existing file of the same name will
			be removed before writing.

		name: Usually a short abbreviated name to associate with the dataset. Will be used as file name
			for .h5 file (so "MyDataset" will create MyDataset.h5 file), hence ensure it's file name compliant.
		attributes: dictionary describing the dataset, with following keys (all values are strings)
			fullname: A more descriptive name for the dataset.
			version: Some version string to associate with the dataset.
			description: Description of the dataset.
			platform_type: one of ['microarray','rna-seq']
			platform_details: descriptive string for the platform, eg. "Illumina WG6 version 2"
			pubmed_id: id string to pubmed if published. Leave empty or None if not applicable.
			species: one of ['MusMusculus','HomoSapiens']
		samples: pandas DataFrame which has sample ids matching the columns of the expression matrix. Example:
			sampleId	celltype	tissue
			s01		LSK		bone marrow
			s02		T1		peripheral blood
			(sampleId must be the index of the DataFrame.)
	
		expression: pandas DataFrame which has sample ids as columns and probe ids as index (microarray only)
			probeId		s01		s02
			ILMN_1234	4.35	8.42
			ILMN_4567	8.11	5.81
		probeIdsFromGeneId: A dictionary of form {'ENSMUSG000004353':['ILMN_1234',...],...} (microarray only)
		geneIdsFromProbeId: A dictionary of form {'ILMN_1234':['ENSMUSG000004353',...],...} (microarray only)
			Both of these dictionaries allow for the possibility of multiple matches.
		
		counts: pandas DataFrame of raw counts summarised at genes (rna-seq only)
		cpm: pandas DataFrame of counts per million values for each gene and sample (rna-seq only)
		rpkm: pandas DataFrame of RPKM values for each gene and sample (rna-seq only)

	Returns:
		A Dataset instance.

	Example:
		>>> attributes = {"fullname": "Haemopedia",
					  "version": "1.0",
					  "description": "Microarray expression of murine haematopoietic cells",
					  "platform_type": "microarray",
					  "platform_details": "Illumina WG6 version 2",
					  "pubmed_id": None,
					  "species": "MusMusculus"}
		# If samples, expression, probeIdsFromGeneId, geneIdsFromProbeId have been set as above:
		>>> createDatasetFile("/datasets", name="haemopedia", attributes=attributes, samples=samples,
			expression=expression, probeIdsFromGeneId=probeIdsFromGeneId, geneIdsFromProbeId=geneIdsFromProbeId)
	"""
	# Parse input
	if destDir[-1]=='/':
		destDir = destDir[:-1]
	name = kwargs['name']
	attributes = kwargs['attributes']
	samples = kwargs['samples']
	
	expression = kwargs.get('expression')
	probeIdsFromGeneId = kwargs.get('probeIdsFromGeneId',{})
	geneIdsFromProbeId = kwargs.get('geneIdsFromProbeId',{})
	
	counts = kwargs.get('counts')
	cpm = kwargs.get('cpm')
	rpkm = kwargs.get('rpkm')

	# Check data
	
	# Save all values to file
	filepath = '%s/%s.h5' % (destDir, name)
	if os.path.isfile(filepath): os.remove(filepath)
	pandas.Series(attributes).to_hdf(filepath, '/series/attributes')
	samples.to_hdf(filepath, '/dataframe/samples')

	if attributes['platform_type']=='microarray':
		# Only keep probeIds which occur in expression matrix
		pgi = dict([(geneId, probeIdsFromGeneId[geneId]) for geneId in probeIdsFromGeneId.keys() \
			if pandas.notnull(geneId) and len([probeId for probeId in probeIdsFromGeneId[geneId] if probeId in expression.index])>0])
		gpi = dict([(probeId, geneIdsFromProbeId[probeId]) for probeId in geneIdsFromProbeId.keys() \
			if probeId in expression.index and len([geneId for geneId in geneIdsFromProbeId[probeId] if pandas.notnull(geneId)])>0])		
		expression.to_hdf(filepath, '/dataframe/expression')
		pandas.Series(pgi).to_hdf(filepath, '/series/probeIdsFromGeneId')
		pandas.Series(gpi).to_hdf(filepath, '/series/geneIdsFromProbeId')
	else:
		counts.to_hdf(filepath, '/dataframe/counts')
		cpm.to_hdf(filepath, '/dataframe/cpm')
		rpkm.to_hdf(filepath, '/dataframe/rpkm')

	return Dataset(filepath)
	
# ------------------------------------------------------------
# Dataset class 
# ------------------------------------------------------------
class Dataset(object):
	"""
	An instance of Dataset is associated with its HDF file, which contains all the information about the dataset,
	including expression matrix, sample information, and various metadata associated with the dataset.

	Attributes:
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
	
		Parameters:
			geneIds: a list of gene ids eg: ['ENSG000003535',...] or a string 'ENSG0000003535' or None, 
				in which case the full expression matrix is returned.
			datatype: one of ['counts','rpkm','cpm'] if self.isRnaSeqData is True, ignored if False.
			sampleGroupForMean: a sample group name eg 'celltype' which will be used to return the mean
				over the sample ids.
				
		Returns:
			pandas.DataFrame instance	
		"""
		df = self._expression[datatype] if self.isRnaSeqData else self._expression
		
		if isinstance(geneIds, str) or isinstance(geneIds, unicode):	# assume single gene was specified
			geneIds = [geneIds]
		
		# work out which row index to use based on platform type
		if geneIds:
			index = set(geneIds).intersection(set(df.index)) if self.isRnaSeqData else \
				set(self.probeIdsFromGeneIds(geneIds=geneIds)).intersection(set(df.index))
				#[probeId for geneId,probeId in self._probeIdsFromGeneId.iteritems() if geneId in geneIds]
		else:
			index = df.index
		
		if len(index)>0:
			df = df.loc[index]
		else:
			return pandas.DataFrame()
					
		if sampleGroupForMean:
			sgi = self.sampleGroupItems(sampleGroup=sampleGroupForMean, duplicates=True)
			if ','.join(sgi)!=','.join(df.columns):	# it's possible for celltypes to be defined the same as sample ids, for example
				df = df.groupby(sgi, axis=1).mean()
			
		return df
		
	def expressionValues(self, geneIds=[]):
		"""
		Return a dictionary of expression values matching geneIds.
		
		Example:
			>>> expressionValues(geneIds=['ENSMUSG00000019982','ENSMUSG00000005672'])
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
			
	def correlation(self, featureId):
		"""
		Return a dictionary of correlation scores for each feature id in the dataset, eg: {'ENSMUSG0000034355':0.343, ...}
		Each value is the correlation between featureId and each feature in the appropriate expression matrix.
		featureId must be a member of the dataset (ie. probe id for microarray data). Otherwise None is returned.
		In this case only the keys contained in this dictionary are returned.
		"""
		from scipy.stats.stats import pearsonr
		
		# The methods are quite different for rna-seq vs microarray datasets
		df = self._expression['rpkm'] if self.isRnaSeqData else self._expression
		if featureId not in df.index: return None
		
		if self.isRnaSeqData:	# use log2 rpkm values
			df = numpy.log2(df)		
		
		# values for featureId
		values = df.loc[featureId]
			
		# loop through all features and calculate correlation
		score = {}
		for featureId,row in df.iterrows():
			corr = pearsonr(values, row)[0]   # only need correlation and not the p value
			if pandas.notnull(corr):   # null can happen if row is all zeros (or all the same values => zero covariance)
				score[featureId] = corr

		if self.isRnaSeqData:
			return score
		else:	# return gene ids as keys instead of probe ids, using highest value probe per gene
			scoreModified = {}
			for geneId,probeIds in self._probeIdsFromGeneId.iteritems():
				scorelist = [score[probeId] for probeId in probeIds if probeId in score]
				if scorelist:
					scoreModified[geneId] = max(scorelist)
			return scoreModified
		
	# ------------------------------------------------------------
	# Probe and genes 
	# ------------------------------------------------------------		
	def probeIdsFromGeneIds(self, geneIds=[], returnType="list"):
		"""Return probe ids matching geneIds.
		
		Parameters:
			geneIds: a list of gene ids
			returnType: {'list', 'dict'}
		
		Returns:
			If returnType=='list': returns a flat list of probe ids
			If returnType=='dict': returns a dictionary of lists, eg: {'gene1':['probe1'],...}
		"""
		# The following is really really slow
		#index = set(self._probeIdsFromGeneId.index).intersection(set(geneIds))
		#return self._probeIdsFromGeneId.loc[index] if len(index)>0 else pandas.Series()
		
		pfgi = self._probeIdsFromGeneId
		if returnType=="list":
			return sum([pfgi[geneId] for geneId in geneIds if geneId in pfgi], [])
		elif returnType=="dict":
			return dict([(geneId, pfgi[geneId]) for geneId in geneIds if geneId in pfgi])
		
	def geneIdsFromProbeIds(self, probeIds=[], returnType="list"):
		"""Return gene ids matching probeIds.

		Parameters:
			probeIds: a list of probe ids
			returnType: {'list', 'dict'}
		
		Returns:
			If returnType=='list': returns a flat list of gene ids
			If returnType=='dict': returns a dictionary of lists, eg: {'probe1':['gene1'],...}
		"""
		gfpi = self._geneIdsFromProbeId
		if returnType=="list":
			return sum([gfpi[probeId] for probeId in probeIds if probeId in gfpi], [])
		elif returnType=="dict":
			return dict([(probeId, gfpi[probeId]) for probeId in probeIds if probeId in gfpi])
				
		
	# ------------------------------------------------------------
	# Sample related methods 
	# ------------------------------------------------------------
	def sampleGroups(self):
		"""Return a list of sample group names eg: ["celltype","tissue"]
		"""
		return self._samples.columns.tolist()

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
		df = self._samples
	
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
		df = self._samples
		if sampleGroup and sampleGroupItem:
			return df[df[sampleGroup]==sampleGroupItem].index.tolist()
		else:
			return df.index.tolist()
	
	def sampleIdsFromSampleGroups(self):
		"""Return a dictionary of sample ids as lists, keyed on sample group items,
		eg: {'celltype': {'B':['s1','s2',...], ...}, ... }
		"""
		dict = {}
		df = self._samples
		for sampleGroup in df.columns:
			dict[sampleGroup] = {}
			for sampleId,value in df[sampleGroup].iteritems():
				if pandas.isnull(value): continue
				if value not in dict[sampleGroup]: dict[sampleGroup][value] = []
				dict[sampleGroup][value].append(sampleId)
		return dict
	
	def sampleTable(self):
		"""Return a pandas DataFrame of the entire sample table.
		"""
		return self._samples
		
# ------------------------------------------------------------
# Tests - eg. nosetests dataset.py
# ------------------------------------------------------------
def test_utilityFunctions():
	# rpkm
	assert rpkm(1022, 34119529, 2566)==11
	
	# probeGeneMap
	pgm = probeGeneMap("IlluminaWG6", probeIds=['ILMN_1212612','ILMN_1213657'])
	assert pgm['probeIdsFromGeneId']['ENSMUSG00000039601']==['ILMN_1212612']
	assert len(pgm['geneIdsFromProbeId']['ILMN_1213657'])==15

def test_metadata():
	ds = Dataset("%s/testdata.h5" % dataDirectory())	
	assert ds.platform_type=="microarray"
	assert not ds.isRnaSeqData
	assert ds.species=="MusMusculus"
	
def test_probeIdGeneIdMapping():
	ds = Dataset("%s/testdata.h5" % dataDirectory())	
	assert ds.probeIdsFromGeneIds(geneIds=['gene1'])==['probe1']
	assert ds.probeIdsFromGeneIds()==[]
	assert ds.probeIdsFromGeneIds(geneIds=['gene2'], returnType="dict")=={'gene2':['probe2','probe3']}
	
def test_expressionMatrix():
	ds = Dataset("%s/testdata.h5" % dataDirectory())
	assert ds.expressionMatrix()['s01'].tolist()==[3.45,5.54,0]
	assert ds.expressionMatrix(sampleGroupForMean='celltype').at['probe2','B1']==2.77
	
def test_expressionValues():
	ds = Dataset("%s/testdata.h5" % dataDirectory())	
	ev = ds.expressionValues(geneIds=['gene2'])
	assert ev['values']['probe2']['s01']==5.54
	assert ev['featureGroups']['gene2']==['probe2','probe3']
	assert ds.expressionValues(geneIds=['nonsense'])['values']=={}
	assert str(ds.correlation('probe1')['gene2'])[:4]=='0.53'

def test_samples():
	ds = Dataset("%s/testdata.h5" % dataDirectory())	
	assert ds.sampleGroups()==['celltype','tissue']
	assert ds.sampleGroupItems(sampleGroup="celltype")==['B1','B2']
	assert ds.sampleGroupItems(sampleGroup='celltype', groupBy='tissue')=={'BM': ['B1', 'B2']}
	assert ds.sampleGroupItems(sampleGroup='celltype', duplicates=True)==['B1','B1','B2','B2']
	assert ds.sampleIds()==['s01','s02','s03','s04']
	assert ds.sampleIds(sampleGroup='celltype', sampleGroupItem='B1')==['s01','s02']
	assert ds.sampleIdsFromSampleGroups()['tissue']['BM']==['s01', 's02', 's03', 's04']
	