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

from lib import VERSION, DESC

extra_link_args = ['-lm']
if sys.platform != 'darwin':
    extra_link_args.append('-s')

srcs = [os.path.join('src/raxml',fn) for fn in os.listdir('src/raxml')
        if (not os.path.isdir(fn)) and fn.endswith('.c')]
raxml_module = Extension('raxmlib._raxml',
                         sources=['lib/raxmlib/raxml.i'] + srcs,
                         extra_link_args=extra_link_args
                         )

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

    package_dir = {'': 'lib'},
    packages=['TreeLib',
              'raxmlib',
              'ga'],
    py_modules=[],
    scripts=['bin/gaperm'],
    ext_modules=[raxml_module]
    )
