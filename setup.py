__author__ = 'aleeee'

import os
from distutils.core import setup

setup(name='dnavec',
      version='0.1',
      description='DNA vectoring in Python.',
      long_description=open('README.md').read(),
      author='Fule Liu',
      author_email='liufule12@gmail.com',
      packages=['dnavec', 'dnavec.nac', 'dnavec.psenac', 'dnavec.ac'],)