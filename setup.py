__author__ = 'aleeee'

import os
from distutils.core import setup


setup(name='bioHelper',
      version='0.1',
      description='bioHelper in Python.',
      long_description=open('README.md').read(),
      author='Fule Liu',
      author_email='liufule12@gmail.com',
      url='https://github.com/liufule12/bio_package',
      packages=['bioHelper', 'bioHelper.kmer', 'bioHelper.pseudo'],)

