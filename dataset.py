import os, sys
import pandas, numpy

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
	return
	
	ds = Dataset("/Users/jchoi/projects/Hematlas/process_dataset/dmap/v1.1/DataOutput/dmap.h5")
	
	# metadata test
	assert ds.name=="dmap"
	assert not ds.isRnaSeqData
	assert ds.species=="HomoSapiens"
	
	# geneIdsFromFeatureIds


# ------------------------------------------------------------
# Utility functions 
# ------------------------------------------------------------
def rpkm(rawCount, totalReads, medianTranscriptLength):
	"""Return RPKM value for a gene. Example: rpkm(1022, 34119529, 2566)
	Formula is (10^9 * C)/(N * L). See https://www.biostars.org/p/55253/
	
	Args:
		rawCount: unnormalised summarised count for a gene
		totalReads: total reads for all genes in a sample - ie. library size
		medianTranscriptLength: median length of all transcripts in the gene.
		
	Returns:
		float
	"""
	return numpy.power(10,9)*rawCount/totalReads/medianTranscriptLength

def dataDirectory():
	"""Return the full path to the directory where data used by this package is kept.
	"""
	return os.path.dirname(os.path.abspath(__file__)) + "/data"
	
def probeGeneMap(arrayType, probeIds=[]):
	"""Return dictionaries of mapping between probe ids and gene ids.

	Args:
		arrayType: one of ['IlluminaWG6']
		probeIds: list of probe ids to restrict the search. Otherwise all probe ids will be searched.
		
	Returns:
		dictionary; see example below
		
	>>> probeGeneMap("IlluminaWG6", probeIds=['ILMN_1212612']
	{'geneIdsFromProbeId': {'ILMN_1212612': ['ENSMUSG00000039601']}, 'nonMatchingProbeIds': [], 'probeIdsFromGeneId': {'ENSMUSG00000039601': ['ILMN_1212612']}}
	
	"""
	filename = "%s/ProbeIdGeneIdMapping_%s.txt" % (dataDirectory(), arrayType)
	df = pandas.read_csv(filename, sep="\t", index_col=0)
	
	if probeIds:  # only interested in this subset
		df = df.loc[probeIds]

	# Work out a list of probe ids for each gene and vice versa
	probeIdsFromGeneId = dict([(geneId,set()) for geneId in df['geneId']])
	geneIdsFromProbeId = dict([(probeId,set()) for probeId in df.index])
	nonMatchingProbeIds = [probeId for probeId in probeIds if probeId not in df.index] if probeIds else []
	
	# Because of multiple matches we have to return a list for each key
	for probeId,row in df.iterrows():
		probeIdsFromGeneId[row['geneId']].add(probeId)
		geneIdsFromProbeId[probeId].add(row['geneId'])
		
	return {'probeIdsFromGeneId': dict([(key,list(probeIdsFromGeneId[key])) for key in probeIdsFromGeneId]),
			'geneIdsFromProbeId': dict([(key,list(geneIdsFromProbeId[key])) for key in geneIdsFromProbeId]),
			'nonMatchingProbeIds': nonMatchingProbeIds}
