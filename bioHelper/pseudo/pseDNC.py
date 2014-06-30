__author__ = 'aleeee'

import sys
sys.path.append('..')
from util import get_sequence_check_dna
from util_pseDNC import make_psednc_vector


class PseDNC():
    def __init__(self, lamada, w, alphabet="ACGT"):
        self.lamada = lamada
        self.w = w
        self.alphabet = alphabet

    def make_vector(self, fasta_file):
        sequence_list = get_sequence_check_dna(open(fasta_file, 'r'))

        vector = make_psednc_vector(sequence_list, self.lamada, self.w)
        for e in vector:
            print e

        return vector


if __name__ == '__main__':
    psednc = PseDNC(1, 0.3)
    vector = psednc.make_vector('hs.txt')
    for e in vector:
        print e
