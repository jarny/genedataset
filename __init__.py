"""
[jchoi@:~/projects/Hematlas/Haemosphere-new]% pip install git+file:/Users/jchoi/projects/Hematlas/
"""

__all__ = ['dataset','geneset']

def dataDirectory():
	"""Return the full path to the directory where data used by this package is kept.
	"""
	import os
	return os.path.dirname(os.path.abspath(__file__)) + "/data"
