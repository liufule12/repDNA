__author__ = 'aleeee'

import sys

sys.path.append('..')
from util import read_fasta_sequence
from util_pseDNC import make_psednc_vector


class PseDNC():
    def __init__(self, lamada, w, alphabet="ACGT"):
        self.lamada = lamada
        self.w = w
        self.alphabet = alphabet

    def make_vector(self, fasta_file):
        sequence_list = read_fasta_sequence(open(fasta_file, 'r'), self.alphabet)

        vector = make_psednc_vector(sequence_list, self.lamada, self.w)
        for e in vector:
            print e

        return vector