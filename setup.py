#!/usr/bin/env python
#
# setup for TreeFix library packages
#
# use the following to install:
#   python setup.py build
#   python setup.py install
#

import os,sys
from distutils.core import setup, Extension
import numpy

sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), "lib")))

USE_CYTHON = True

try:
    from Cython.Distutils import build_ext
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False
        
cmdclass = { }
recon_module = []


VERSION = "1.0.1rc"
DESC = "GATC (Genetic Algorithm for Tree Construction/Correction) find the best tree from a list of candidate trees according to sequence likelihood and reconciliation with a species tree."

extra_link_args = ['-lm']
if sys.platform != 'darwin':
    extra_link_args.append('-s')

srcs = [os.path.join('src/raxml',fn) for fn in os.listdir('src/raxml')
        if (not os.path.isdir(fn)) and fn.endswith('.c')]
raxml_module = [ Extension('lib.raxmlib._raxml',
                         sources=['lib/raxmlib/raxml.i'] + srcs,
                         extra_link_args=extra_link_args
                         )]
if USE_CYTHON:
    recon_module += [
        Extension("lib.reclkl.computeLKL", sources=[ "src/recon/computeLKL.pyx"], include_dirs=[numpy.get_include()]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    recon_module += [
        Extension("lib.reclkl.computeLKL", sources=[ "src/recon/computeLKL.c" ], include_dirs=[numpy.get_include()]),
    ]

modules = raxml_module + recon_module

setup(
    name='GATC',
    version = VERSION, 
    description=DESC,

    author='Emmanuel Noutahi',
    author_email='noutahie@iro.umontreal.ca',

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Education',
        ],

    packages=['lib', 'lib.TreeLib', 'lib.raxmlib', 'lib.ga', 'lib.ga.evolve', 'lib.reclkl', 'lib.PolytomySolver'],
    py_modules=[],
    scripts=['bin/gatc'],
    install_requires=['numpy', 'scipy', 'matplotlib', 'ete3', 'biopython'],
    cmdclass=cmdclass,
    ext_modules=modules
    )
