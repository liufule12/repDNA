__author__ = 'Fule Liu'

from distutils.core import setup

setup(name='repDNA',
      version='1.0',
      author='Bin Liu and Fule Liu',
      author_email='bliu@insun.hit.edu.cn and liufule12@gmail.com',
      maintainer='Fule Liu',
      maintainer_email='liufule12@gmail.com',
      url='bioinformatics.hitsz.edu.cn/repDNA/',
      description='a Python package to generate various feature vectors of DNA sequences incorporating physicochemical properties and sequence-order effects',
      download_url='bioinformatics.hitsz.edu.cn/repDNA/download',
      platforms=['MS Windows', 'Mac OS', 'Unix/Linux'],
      license='GPL',
      packages=['repDNA'],
      package_data={'repDNA': ['data/*.data', 'example/*.fasta', 'example/*.py']},
)