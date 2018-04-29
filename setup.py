from setuptools import setup, Extension


setup(name = 'SMARTc', version = '1.0', \
      ext_modules = [Extension('SMARTc', ['C:/PycharmProjects/C/SMARTc/source/SMARTc.cpp'])])
