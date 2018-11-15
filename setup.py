# -*- coding: utf-8 -*-
# Copyright (C) 2018  Thibault Hallouin
from setuptools import setup, Extension
import numpy


with open("README.md", "r") as fh:
    long_desc = fh.read()

with open('smartcpp/version.py') as fv:
    exec(fv.read())

setup(
    name='smartcpp',

    version=__version__,

    description='SMARTcpp: a C++ extension of the rainfall-runoff SMART for Python',
    long_description=long_desc,
    long_description_content_type="text/markdown",

    url='https://github.com/ThibHlln/smartcpp',

    author='Thibault Hallouin, Michael Bruen, and Eva Mockler',
    author_email='thibault.hallouin@ucdconnect.ie',

    license='GPLv3',

    classifiers=[
        'Development Status :: 4 - Beta',

        'Natural Language :: English',

        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering',

        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows :: Windows 10',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: C++',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython'
    ],

    packages=['smartcpp'],

    ext_modules=[Extension('smartcpp', ['smartcpp/smartcpp.cpp'],
                           include_dirs=[numpy.get_include()])],

    python_requires='>=2.7'
)
