from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from scorpio import __version__, _program

setup(name='scorpio',
      version=__version__,
      packages=find_packages(),
      scripts=["scorpio/scripts/scorpiofunks.py",
               "scorpio/scripts/type_constellations.py"],
      package_data={"scorpio":["data/*"]},
      install_requires=[
            "biopython>=1.70",
            "pytools>=2020.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4',
            "matplotlib>=3.2.1",
            "scipy>=1.4.1",
            "numpy>=1.13.3",
            "geopandas>=0.7.0",
            "descartes>=1.1.0",
            "adjustText>=0.7.3",
            "tabulate>=0.8.7",
            "snipit>=1.0.3",
            "seaborn>=0.10.1",
            "epiweeks>=2.1.2"
        ],

      description='serious constellations of reoccurring phylogenetically-independent origin',
      url='https://github.com/cov-lineages/scorpio',
      author='Ben Jackson & Rambaut Group',
      author_email='bjackso4@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = scorpio.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
