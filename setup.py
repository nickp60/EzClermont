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

# Get the long description from the README file
# with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#     long_description = f.read()

# # from pyani
# with open(path.join(here, 'clermontpcr', '__init__.py'), 'r') as fh:
#     for line in fh:
#         m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
#         if m:
#             init_version = m.group('version')
#             break
# from http://stackoverflow.com/questions/458550/
#         standard-way-to-embed-version-into-python-package
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

## parse requirements file
# install_reqs = parse_requirements(path.join(here, "requirements.txt"),
#                                  session=False)
# requirements = [str(ir.req) for ir in install_reqs]
setup(
    name='clermontpcr',
    # dependency_links=['https://github.com/WarrenWeckesser/heatmapcluster.git#egg=heatmapcluster'],
    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=verstr,

    description='clermontpcr: phylotype your strains, in silico',
    # long_description=long_description,
    long_description="""
    check out the GitHub
    repo for the real README.md file
    """,

    # The project's main homepage.
    url='https://github.com/nickp60/clermontpcr',

    # Author details
    author='Nick Waters',
    author_email='nickp60@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        #'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.6',
        #'Programming Language :: Python :: 2.7',
        #'Programming Language :: Python :: 3',
        #'Programming Language :: Python :: 3.3',
        #'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # What does your project relate to?
    keywords='bioinformatics evolution genomics development',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    # packages=['clermontpcr'],
    packages=find_packages(),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    # install_requires=requirements,
    install_requires=[
        'Biopython==1.68'
    ],
    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    #extras_require={
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},
    include_package_data=True,
    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
       '': [path.join(__name__, "tests", "refs/*")],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
    # Using the scripts keyword here for speed, as the entry points are not
    # yet defined
    scripts=['clermontpcr/clermontpcr.py'],
)
