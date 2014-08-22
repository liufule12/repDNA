__author__ = 'aleeee'

import sys

from repDNA.nacutil import make_upto_kmer_list, make_revcomp_kmer_list, make_kmer_vector
from repDNA.util import get_data


class Kmer():
    def __init__(self, k=1, normalize=False, upto=False, alphabet="ACGT"):
        self.k = k
        self.upto = upto
        self.normalize = normalize
        self.alphabet = alphabet

    def make_kmer_vec(self, data):
        """Make a kmer vector with options k, upto, revcomp, normalize.

        :param data: file object or sequence list.
        :return: kmer vector.
        """
        sequence_list = get_data(data)

        if self.upto:
            k_list = range(1, self.k + 1)
        else:
            k_list = range(self.k, self.k + 1)
        kmer_list = make_upto_kmer_list(k_list, self.alphabet)

        rev_kmer_list = []
        revcomp = False
        vec = make_kmer_vector(sequence_list, kmer_list, rev_kmer_list, self.k, self.upto, revcomp, self.normalize)
        return vec


class RevcKmer(Kmer):
    def make_revckmer_vec(self, data):
        """Make a reverse compliment kmer vector with options k, upto, normalize.

        :param data: file object or sequence list.
        :return: reverse compliment kmer vector.
        """
        sequence_list = get_data(data)

        if self.upto:
            k_list = range(1, self.k + 1)
        else:
            k_list = range(self.k, self.k + 1)
        kmer_list = make_upto_kmer_list(k_list, self.alphabet)

        # Use lexicographically first version of {kmer, revcomp(kmer)}.
        rev_kmer_list = make_revcomp_kmer_list(kmer_list)
        revcomp = True
        vec = make_kmer_vector(sequence_list, kmer_list, rev_kmer_list, self.k, self.upto, revcomp, self.normalize)
        return vec


class IDkmer():
    def __init__(self, k=6, upto=True):
        if k <= 0:
            error_info = 'Sorry, the k value must be an integer larger than 0.'
            sys.stdout.write(error_info)
            sys.exit(0)
        if not isinstance(upto, bool):
            error_info = 'Sorry, the upto type must be bool.'
            sys.stdout.write(error_info)
            sys.exit(0)
        self.k = k
        self.upto = upto

    def make_idkmer_vec(self, data, hs, non_hs):
        from repDNA.nacutil import make_kmer_list
        from repDNA.nacutil import diversity
        from repDNA.nacutil import id_x_s

        alphabet, rev_kmer_list, upto, revcomp, normalize = 'ACGT', [], False, False, False

        pos_s_list = get_data(hs)
        neg_s_list = get_data(non_hs)
        print self.k
        if self.upto is False:
            k_list = [self.k]
        else:
            k_list = range(1, self.k+1)

        print 'k_list =', k_list

        # Get all kmer ID from 1-kmer to 6-kmer.
        # Calculate standard source S vector.
        pos_s_vec, neg_s_vec = [], []
        diversity_pos_s, diversity_neg_s = [], []
        for k in k_list:
            kmer_list = make_kmer_list(k, alphabet)

            temp_pos_s_vec = make_kmer_vector(pos_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
            temp_neg_s_vec = make_kmer_vector(neg_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)

            temp_pos_s_vec = [sum(e) for e in zip(*[e for e in temp_pos_s_vec])]
            temp_neg_s_vec = [sum(e) for e in zip(*[e for e in temp_neg_s_vec])]

            pos_s_vec.append(temp_pos_s_vec)
            neg_s_vec.append(temp_neg_s_vec)

            diversity_pos_s.append(diversity(temp_pos_s_vec))
            diversity_neg_s.append(diversity(temp_neg_s_vec))

        # Calculate Diversity(X) and ID(X, S).
        sequence_list = get_data(data)
        vec = []

        for seq in sequence_list:
            # print seq
            temp_vec = []
            for k in k_list:
                kmer_list = make_kmer_list(k, alphabet)
                seq_list = [seq]
                kmer_vec = make_kmer_vector(seq_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
                print 'k', k
                print 'kmer_vec', kmer_vec

                # print diversity_pos_s
                if upto is False:
                    k = 1

                print 'pos_vec', pos_s_vec
                print 'neg_vec', neg_s_vec
                print 'diversity_pos_s', diversity_pos_s

                temp_vec.append(round(id_x_s(kmer_vec[0], pos_s_vec[k-1], diversity_pos_s[k-1]), 3))
                temp_vec.append(round(id_x_s(kmer_vec[0], neg_s_vec[k-1], diversity_neg_s[k-1]), 3))

            vec.append(temp_vec)

        return vec


if __name__ == '__main__':
    from repDNA.nac import Kmer

    kmer = Kmer(k=2)
    vec = kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, normalize=True)
    vec = kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, normalize=False, upto=True)
    vec = kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    from repDNA.nac import RevcKmer

    revckmer = RevcKmer(k=2, normalize=False, upto=False)
    vec = revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    revckmer = RevcKmer(k=2, normalize=True, upto=False)
    vec = revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    revckmer = RevcKmer(k=2, normalize=True, upto=True)
    vec = revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    print 'Begin IDkmer.'
    from repDNA.nac import IDkmer

    print 'Test: default mod.'
    idkmer = IDkmer()
    vec = idkmer.make_idkmer_vec(open('test.txt'), open('hs.txt'), open('non-hs.txt'))
    print vec
    print

    print 'Test: k=2.'
    idkmer = IDkmer(k=2)
    vec = idkmer.make_idkmer_vec(open('test.txt'), open('hs.txt'), open('non-hs.txt'))
    print vec
    print

    print 'Test: k=2, upto=False'
    idkmer = IDkmer(k=2, upto=False)
    vec = idkmer.make_idkmer_vec(open('test.txt'), open('hs.txt'), open('non-hs.txt'))
    print vec
    print


    # x = [11, 20, 13, 9, 27, 17, 1, 16, 9, 9, 14, 10, 6, 16, 13, 41]
    # s = [68, 67, 67, 58, 87, 52, 7, 76, 56, 46, 69, 48, 49, 58, 73, 90]
    # d_s = 3797.6619762268665
    # from repDNA.nacutil import id_x_s
    # print id_x_s(x, s, d_s)