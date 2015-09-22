import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.txt')) as f:
    README = f.read()
with open(os.path.join(here, 'CHANGES.txt')) as f:
    CHANGES = f.read()

requires = [
    'pandas',
    ]

setup(name='genedataset',
      version='0.1',
      description='Store and access gene expression datasets and gene definitions.',
      long_description=README + '\n\n' + CHANGES,
      classifiers=[
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
      author='Jarny Choi',
      author_email='jchoi@wehi.edu.au',
      url='https://github.com/jarny/genedataset',
      keywords='gene microarray rna-seq',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
      install_requires=requires,
      tests_require=requires,
      )