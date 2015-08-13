import os, sys
import pandas, numpy

# ------------------------------------------------------------
# tests 
# ------------------------------------------------------------
def test_Dataset():
	"""
	nosetests datasets.py
	"""
	# rpkm
	assert rpkm(1022, 34119529, 2566)==11

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
