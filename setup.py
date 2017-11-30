"""
Setup for clermontpcr

A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
import re
from codecs import open
from os import path
import sys
from pip.req import parse_requirements

here = path.abspath(path.dirname(__file__))

VERSIONFILE = "clermontpcr/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: clermontpcr requires Python 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

setup(
    name='clermontpcr',
    version=verstr,

    description='clermontpcr: phylotype your strains, in silico',
    # long_description=long_description,
    long_description="""
    check out the GitHub
    repo for the real README.md file
    """,

    url='https://github.com/nickp60/clermontpcr',

    # Author details
    author='Nick Waters',
    author_email='nickp60@gmail.com',
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics evolution genomics development',
    packages=find_packages(),
    install_requires=[
        'Biopython==1.68'
    ],
    include_package_data=True,
    package_data={
       '': [path.join(__name__, "tests", "refs/*")],
    },
    scripts=['clermontpcr/clermontpcr.py'],
)
