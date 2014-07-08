__author__ = 'aleeee'

from dnavec.kmer.kmerutil import make_upto_kmer_list, make_revcomp_kmer_list, make_kmer_vector
from dnavec.util import get_data


class Kmer():
    def __init__(self,
                 k=1,
                 upto=False,
                 revcomp=False,
                 normalize=False,
                 alphabet="ACGT"):
        self.k = k
        self.upto = upto
        self.revcomp = revcomp
        self.normalize = normalize
        self.alphabet = alphabet

    def make_kmer_vector(self, data):
        """Make a kmer vector with options k, upto, revcomp, normalize.

        :param data: file object or sequence list.
        :return: kmer vector.

        Examples
        --------
        >>> from dnavec.kmer.kmer import Kmer

        >>> kmer = Kmer(k=2, upto=False, revcomp=False, normalize=False)
        >>> vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
        The kmer is ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        >>> vec
        [[1, 3, 0, 3, 2, 0, 0, 4, 2, 2, 1, 1, 2, 2, 4, 7]]

        >>> kmer = Kmer(k=2, upto=True, revcomp=False, normalize=False)
        >>> vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
        The kmer is ['A', 'C', 'G', 'T', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA',
        'TC', 'TG', 'TT']
        >>> vec
        [[7, 7, 6, 15, 1, 3, 0, 3, 2, 0, 0, 4, 2, 2, 1, 1, 2, 2, 4, 7]]

        >>> kmer = Kmer(k=2, upto=False, revcomp=True, normalize=False)
        >>> vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
        Reduced to 10 kmers.
        The kmer is ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        >>> vec
        [[8, 4, 4, 3, 6, 1, 0, 4, 2, 2]]

        >>> kmer = Kmer(k=2, upto=False, revcomp=False, normalize=True)
        >>> vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
        The kmer is ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
        >>> vec
        [[0.029411764705882353, 0.08823529411764706, 0.0, 0.08823529411764706, 0.058823529411764705, 0.0, 0.0,
        0.11764705882352941, 0.058823529411764705, 0.058823529411764705, 0.029411764705882353, 0.029411764705882353,
         0.058823529411764705, 0.058823529411764705, 0.11764705882352941, 0.20588235294117646]]

        """
        sequence_list = get_data(data)

        if self.upto:
            k_list = range(1, self.k + 1)
        else:
            k_list = range(self.k, self.k + 1)
        kmer_list = make_upto_kmer_list(k_list, self.alphabet)

        # Use lexicographically first version of {kmer, revcomp(kmer)}.
        rev_kmer_list = []
        if self.revcomp:
            rev_kmer_list = make_revcomp_kmer_list(kmer_list)

        vec = make_kmer_vector(sequence_list, kmer_list, rev_kmer_list, self.k, self.upto, self.revcomp,
                               self.normalize)
        return vec


if __name__ == '__main__':
    from dnavec.kmer.kmer import Kmer

    kmer = Kmer(k=2, upto=False, revcomp=False, normalize=False)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, upto=True, revcomp=False, normalize=False)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, upto=False, revcomp=True, normalize=False)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, upto=False, revcomp=False, normalize=True)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    print get_data(open('hs.txt'))