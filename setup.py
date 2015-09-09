#!/usr/bin/env python
"""GFam -- installer script"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010-2015, Tamas Nepusz"
__license__ = "GPL"

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

params = {}
params["name"] = "gfam"
params["version"] = "1.3"
params["description"] = "Genome Families"

params["packages"] = find_packages(exclude='tests')
params["scripts"] = ["bin/gfam"]

setup(**params)

