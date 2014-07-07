__author__ = 'aleeee'

import os
from distutils.core import setup


setup(name='dnavec',
      version='0.1',
      description='DNA vectoring in Python.',
      long_description=open('README.md').read(),
      author='Fule Liu',
      author_email='liufule12@gmail.com',
      url='https://github.com/liufule12/bio_package',
      packages=['dnavec', 'dnavec.kmer', 'dnavec.psednc', 'dnavec.pseknc', 'dnavec.inucpseknc'],)