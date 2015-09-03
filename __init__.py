"""
"""

__all__ = ['dataset','geneset']

def dataDirectory():
	"""Return the full path to the directory where data used by this package is kept.
	"""
	import os
	return os.path.dirname(os.path.abspath(__file__)) + "/data"
