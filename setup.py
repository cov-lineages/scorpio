from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os

from scorpio import __version__, _program

setup(name='scorpio',
      version=__version__,
      packages=find_packages(),
      scripts=["scorpio/scripts/type_constellations.py",
               "scorpio/scripts/extract_definitions.py",
               "scorpio/scripts/definitions.py"],
      install_requires=[
            "biopython>=1.70"
        ],
      description='serious constellations of reoccurring phylogenetically-independent origin',
      url='https://github.com/cov-lineages/scorpio',
      author='Rachel Colquhoun, Ben Jackson & Rambaut Group',
      author_email='rachel.colquhoun@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = scorpio.__main__:main
      """.format(program=_program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
