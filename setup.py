__author__ = 'aleeee'

from distutils.core import setup

setup(name='repDNA',
      version='0.1',
      description='A Python package to generate various representations for DNA sequences',
      long_description=open('README.md').read(),
      author='Fule Liu',
      author_email='liufule12@gmail.com',
      packages=['repDNA'],
      package_data={'repDNA': ['data/*.data', 'example/*.fasta', 'example/*.py']},
      )