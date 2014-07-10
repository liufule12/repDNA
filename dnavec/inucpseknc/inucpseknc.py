__author__ = 'aleeee'

from dnavec.util import get_data
from dnavec.inucpseknc.inucpsekncutil import make_pseknc_vector


class INucPseKNC():
    """This class should be used to make iNucPseKNC vector."""

    def __init__(self, k=4, lamada=6, w=0.5):
        """Initialize a new INucPseKNC object.

        :param k: k-tuple.
        :param lamada: The number of the total counted ranks(or tiers) of the correlations along a DNA sequence.
        :param w: The weight factor.
        """
        self.k = k
        self.lamada = lamada
        self.w = w

    def make_vector(self, data, upto=False, revcomp=False):
        """Make iNucPseKNC vector.

        :param data: The fasta file path or single DNA sequence or DNA sequence list.
        :return: vector: The iNucPseKNC vector.
        """

        sequence_list = get_data(data)

        vector = make_pseknc_vector(sequence_list, self.k, self.lamada, self.w, upto, revcomp)

        return vector


if __name__ == '__main__':
    iNuPseKNC = INucPseKNC(3, 1, 0.05)
    # res = iNuPseKNC.make_vector(open('hs.txt'), upto=False, revcomp=True)
    res = iNuPseKNC.make_vector(open('hs.txt'))
    for e in res:
        print e
        print len(e)