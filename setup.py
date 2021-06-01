#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

import mcroni

from setuptools import setup

setup(
    name='mcroni',
    version=flanker.__version__,
    description=' mcr-1 analysis ',
    url='https://github.com/liampshaw/mcroni',
    license='LICENSE',
    python_requires='>=3.6',
    packages=['mcroni'],
    install_requires=['pandas', 'biopython', 'numpy'],
    entry_points={'console_scripts':['mcroni=mcroni.mcroni:main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English'])