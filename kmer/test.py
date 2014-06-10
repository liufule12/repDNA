__author__ = 'aleeee'

import sys
from kmer import Kmer

k, upto, revcomp, normalize = 1, False, False, False
if len(sys.argv[1]) != 0:
    k = int(sys.argv[1])
if sys.argv[2] == '-upto':
    upto = True
if sys.argv[3] == '-revcomp':
    revcomp = True
if sys.argv[4] == '-normalize':
    normalize = True

kmer = Kmer(k, upto, revcomp)
print "k:", kmer.k_value, "upto:", kmer.upto, "revcomp:", kmer.revcomp, "normalize:", kmer.normalize, "alphabet:", kmer.alphabet
vector = kmer.make_kmer_vector("hs.txt")
for e in vector:
    print e

