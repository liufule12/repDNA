## 1. Introduction

### 1.1 repDNA description

In the field of bioinformatics, more and more computational methods employed well-known machine learning techniques to construct their classifiers, such as Support Vector Machines (SVMs), Random Forest (RF), Artificial Neural Networks (ANN), etc. These methods require fixed length feature vectors as inputs. Therefore, many computational tools have been proposed to generate various representations for protein sequences (Cao et al. 2013, propy: a tool to generate various modes of Chou’s PseAAC), but little work has been done on developing useful computational tools to calculate the representations of DNA sequences. To facilitate the studies of DNA and nucleic acids, repDNA was proposed, which is a freely available Python package to generate various representations for DNA sequences. To our best knowledge, this is the first Python tool focusing on DNA representations.

### 1.2 Variou DNA representations in repDNA
repDNA computes three feature groups composed of 12 features, which will described in following part of this manual.
**Basic function**
* Read sequence data from FASTA files

**Nucleic acid Composition**
* Basic kmer
* Reverse compliment kmer.
* Increment of Diversity

**Autocorrelation**
* Auto covariance (AC)
* Cross covariance (CC)
* Auto-cross covariance (ACC)

**Pseudo nucleotide Composition**
* Pseudo dinucleotide composition (PseDNC)
* Pseudo k-tupler nucleotide composition (PseKNC)
* Parallel correlation pseudo dinucleotide composition (PC-PseDNC)
* Parallel correlation pseudo trinucleotide composition (PC-PseTNC)
* Series correlation pseudo dinucleotide composition (SC- PseDNC)
* Series correlation pseudo trinucleotide composition (SC- PseTNC)

### 1.3 Installing repDNA

#### 1.3.1 Requirements

We currently recommend using Python 2.7 from http://www.python.org which is the final version of Python 2.
repDNA is currently supported and tested on the Python 2.7.

#### 1.3.2 For Linux:

To install repDNA, download and unzip the repDNA package, go to the directory at the command line, and type:

    sudo python setup.py install or python setup.py install

#### 1.3.3 For Windows:
To install repDNA, download and unzip the repDNA package, go to the directory at the command line, and type:

    python setup.py install
