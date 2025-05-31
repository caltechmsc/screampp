#!/usr/bin/env python

"""
setup.py for scream
"""

from distutils.core import setup, Extension
import os, glob

# list of src files
libfiles = glob.glob("src/*.cpp")+["swig/py_scream_ee.i"]

scream_library = Extension("_py_scream_ee",
                           sources = libfiles,
                           include_dirs = ["include"],
                           swig_opts=['-c++', '-Iinclude'] 
                           )

setup(name = "scream",
      version = "7Jul11",
      author = "Vaclav Cvicek",
      author_email = "vcvicek@caltech.edu",
      url = "http://www.caltech.edu",
      description = """scream from Victor Kam""",
      py_modules = ['swig.py_scream_ee','python.scream','python.scream_wrap','python.timing'],
      ext_modules = [scream_library]
      )

