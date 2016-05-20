#!/usr/bin/env python
from distutils.core import setup
import sys


# Check that Python 3.2+ is installed
if sys.version_info < (3.0):
    print('ERROR: SparCC is only supported up to python 2.x. Python %d.%d detected'  % sys.version_info[:2])
    print('):')
    sys.exit(1)
# Installed modules
PACKAGES = ['sparcc']
# Required packages
REQUIRED = ['numpy']
# Get and set version from package __init__.py
__version__ = 'Undefined'
for line in open('sparcc/__init__.py'):
    if line.startswith('__version__'):
        exec(line.strip())
# Call setup
setup(name='sparcc',
      version=__version__,
      description='Correlations from sparse compositional data',
      author='yonatanf (modified by Sptehen Watts)',
      licence='mit',
      url='https://github.com/scwatts/sparcc',
      scripts=['scripts/MakeBootstraps.py', 'scripts/PseudoPvals.py',
               'scripts/SampleDist.py', 'scripts/SparCC.py'],
      install_requires=REQUIRED,
      packages=PACKAGES,)
