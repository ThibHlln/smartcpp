from setuptools import setup, Extension


setup(name = 'SMARTc', version = '1.1', \
      ext_modules = [Extension('SMARTc', ['SMARTc.cpp'])])
