__author__ = 'aleeee'

import sys

sys.path.append('..')
from util_PseKNC import make_type1_vector, make_type2_vector
from util import get_sequence_check_dna


class PseKNC():
    def __init__(self, k, lamada, w, alphabet='ACGT'):
        self.k = k
        self.lamada = lamada
        self.w = w
        self.alphabet = alphabet

    def make_vector(self, vec_type):
        sequence_list = get_sequence_check_dna(open('hs.txt', 'r'))

        if self.k == 2:
            phyche_list = ['Tilt', 'Twist']
        elif self.k == 3:
            phyche_list = ['Dnase I']
        else:
            sys.stderr.write("The value of k must be 2 or 3.")
            sys.exit(0)

        if vec_type == 1:
            vector = make_type1_vector(sequence_list, self.k, self.lamada, self.w, self.alphabet, phyche_list)
        elif vec_type == 2:
            vector = make_type2_vector(sequence_list, self.k, self.lamada, self.w, self.alphabet, phyche_list)

        return vector


if __name__ == '__main__':
    # k = int(sys.argv[1])
    # lamada = int(sys.argv[2])
    # w = float(sys.argv[3])
    # vec_type = int(sys.argv[4])
    k = 2
    lamada = 1
    w = 0.05
    vec_type = 1
    pseKNC = PseKNC(k, lamada, w)
    vector = pseKNC.make_vector(vec_type)
    for e in vector:
        print e
        print len(e)
    print len(vector)