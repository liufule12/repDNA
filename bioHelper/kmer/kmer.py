__author__ = 'aleeee'

from bioHelper.kmer.util_kmer import make_upto_kmer_list, make_revcomp_kmer_list, make_kmer_vector
from bioHelper.util import read_fasta_sequence


class Kmer():
    def __init__(self,
                 k_value=1,
                 upto=False,
                 revcomp=False,
                 normalize=False,
                 alphabet="ACGT"):
        self.k_value = k_value
        self.upto = upto
        self.revcomp = revcomp
        self.normalize = normalize
        self.alphabet = alphabet

    def make_kmer_vector(self, fasta_file):
        sequence_list = read_fasta_sequence(open(fasta_file, 'r'), self.alphabet)

        if self.upto:
            k_values = range(1, self.k_value + 1)
        else:
            k_values = range(self.k_value, self.k_value + 1)
        kmer_list = make_upto_kmer_list(k_values, self.alphabet)

        # Use lexicographically first version of {kmer, revcomp(kmer)}.
        rev_kmer_list = []
        if self.revcomp:
            rev_kmer_list = make_revcomp_kmer_list(kmer_list)

        vector = make_kmer_vector(sequence_list, kmer_list, rev_kmer_list, self.k_value, self.upto, self.revcomp,
                                  self.normalize)
        return vector



