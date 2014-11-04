__author__ = 'Fule Liu'

from distutils.core import setup

setup(name='repDNA',
      version='1.1.1',
      author='Bin Liu and Fule Liu',
      author_email='bliu@insun.hit.edu.cn and liufule12@gmail.com',
      maintainer='Fule Liu',
      maintainer_email='liufule12@gmail.com',
      url='bioinformatics.hitsz.edu.cn/repDNA/',
      description='a Python package to generate various modes of feature vectors for DNA sequences '
                  'by incorporating user-defined physicochemical properties and sequence-order effects',
      download_url='bioinformatics.hitsz.edu.cn/repDNA/download',
      platforms=['MS Windows', 'Mac X', 'Unix/Linux'],
      license='GPL',
      packages=['repDNA'],
      package_data={
          'repDNA': ['data/*.data', 'doc/*.pdf', 'example/*.*', 'test/*.*',
                     'data/12_trinucleotide_physicochemical_indices/*.txt',
                     'data/38_dinucleotide_physicochemical_indices/*.txt',
                     'data/6_dinucleotide_physicochemical_indices/*.txt']},
      data_files = [('', ['LICENSE', 'README.md', 'News.md'])]
)