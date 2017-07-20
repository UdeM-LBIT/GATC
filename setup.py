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

from lib import VERSION, DESC

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
        Extension("lib.reclkl.computeLKL", sources=[ "src/recon/computeLKL.pyx" ]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    recon_module += [
        Extension("lib.reclkl.computeLKL", sources=[ "src/recon/computeLKL.c" ]),
    ]

modules = raxml_module + recon_module
print modules

setup(
    name='gaperm',
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

    packages=['lib', 'lib.TreeLib', 'lib.raxmlib', 'lib.ga', 'lib.reclkl'],
    py_modules=[],
    scripts=['bin/gaperm'],
    install_requires=['scipy', 'numpy', 'ete3', 'biopython'],
    cmdclass=cmdclass,
    ext_modules=modules
    )
