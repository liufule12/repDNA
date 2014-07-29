__author__ = 'aleeee'

from dnavec.nac.kmerutil import make_upto_kmer_list, make_revcomp_kmer_list, make_kmer_vector

from dnavec.util import get_data


class Kmer():
    def __init__(self,
                 k=1,
                 normalize=False,
                 upto=False,
                 alphabet="ACGT"):
        self.k = k
        self.upto = upto
        self.normalize = normalize
        self.alphabet = alphabet

    def make_kmer_vector(self, data):
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
        vec = make_kmer_vector(sequence_list, kmer_list, rev_kmer_list, self.k, self.upto, revcomp,
                               self.normalize)
        return vec


class RevcKmer(Kmer):
    def make_revckmer_vector(self, data):
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
        print rev_kmer_list
        revcomp = True
        vec = make_kmer_vector(sequence_list, kmer_list, rev_kmer_list, self.k, self.upto, revcomp,
                               self.normalize)
        return vec


class IDkmer():
    def __init__(self, k=0):
        self.k = k

    def make_idkmer_vector(self, data, hs, non_hs):
        from dnavec.nac.kmerutil import make_kmer_list
        from dnavec.nac.kmerutil import diversity
        from dnavec.nac.kmerutil import id_x_s

        alphabet, rev_kmer_list, upto, revcomp, normalize = 'ACGT', [], False, False, False

        pos_s_list = get_data(hs)
        neg_s_list = get_data(non_hs)

        if self.k > 0:
            # Get a kmer ID.
            # Calculate standard source S vector.
            kmer_list = make_kmer_list(self.k, alphabet)
            pos_s_vec = make_kmer_vector(pos_s_list, kmer_list, rev_kmer_list, self.k, upto, revcomp, normalize)
            neg_s_vec = make_kmer_vector(neg_s_list, kmer_list, rev_kmer_list, self.k, upto, revcomp, normalize)
            pos_s_vec = [sum(e) for e in zip(*[e for e in pos_s_vec])]
            neg_s_vec = [sum(e) for e in zip(*[e for e in neg_s_vec])]
            diversity_pos_s = diversity(pos_s_vec)
            diversity_neg_s = diversity(neg_s_vec)

            # Calculate Diversity(X) and ID(X, S).
            sequence_list = get_data(data)
            vec = []
            for seq in sequence_list:
                # print seq
                temp_vec = []
                kmer_list = make_kmer_list(self.k, alphabet)
                seq_list = [seq]
                kmer_vec = make_kmer_vector(seq_list, kmer_list, rev_kmer_list, self.k, upto, revcomp, normalize)
                temp_vec.append(id_x_s(kmer_vec[0], pos_s_vec, diversity_pos_s))
                temp_vec.append(id_x_s(kmer_vec[0], neg_s_vec, diversity_neg_s))
                # print 'id_x_s', id_x_s(kmer_vec[0], pos_s_vec[k - 1])
                vec.append(temp_vec)
        elif self.k == 0:
            # Get all kmer ID from 1-kmer to 6-kmer.
            # Calculate standard source S vector.
            pos_s_vec, neg_s_vec = [], []
            diversity_pos_s, diversity_neg_s = [], []
            for k in range(1, 7):
                kmer_list = make_kmer_list(k, alphabet)
                temp_pos_s_vec = make_kmer_vector(pos_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
                # print k, '--------------------------------'
                # print temp_pos_s_vec
                temp_neg_s_vec = make_kmer_vector(neg_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
                temp_pos_s_vec = [sum(e) for e in zip(*[e for e in temp_pos_s_vec])]
                # print '+++++++++++++++++++'
                # print temp_pos_s_vec
                temp_neg_s_vec = [sum(e) for e in zip(*[e for e in temp_neg_s_vec])]
                pos_s_vec.append(temp_pos_s_vec)
                neg_s_vec.append(temp_neg_s_vec)

                diversity_pos_s.append(diversity(temp_pos_s_vec))
                diversity_neg_s.append(diversity(temp_neg_s_vec))

            # print 'diversity_pos_s', diversity_pos_s
            # print 'diversity_neg_s', diversity_neg_s

            # Calculate Diversity(X) and ID(X, S).
            sequence_list = get_data(data)
            vec = []
            for seq in sequence_list:
                # print seq
                temp_vec = []
                for k in range(1, 7):
                    kmer_list = make_kmer_list(k, alphabet)
                    seq_list = [seq]
                    kmer_vec = make_kmer_vector(seq_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
                    # print k, kmer_vec
                    # if k == 1:
                    #     print kmer_vec
                    temp_vec.append(id_x_s(kmer_vec[0], pos_s_vec[k - 1], diversity_pos_s[k-1]))
                    temp_vec.append(id_x_s(kmer_vec[0], neg_s_vec[k - 1], diversity_neg_s[k-1]))
                    # print 'id_x_s', id_x_s(kmer_vec[0], pos_s_vec[k - 1])
                vec.append(temp_vec)
        else:
            import sys
            error_info = 'Sorry, the k value must be a integer large than 0.'
            sys.stdout.write(error_info)
            sys.exit(0)

        return vec


if __name__ == '__main__':
    from dnavec.nac.kmer import Kmer

    kmer = Kmer(k=2)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, normalize=True)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    kmer = Kmer(k=2, normalize=False, upto=True)
    vec = kmer.make_kmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    from dnavec.nac.kmer import RevcKmer

    revckmer = RevcKmer(k=2, normalize=False, upto=False)
    vec = revckmer.make_revckmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    revckmer = RevcKmer(k=2, normalize=True, upto=False)
    vec = revckmer.make_revckmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    revckmer = RevcKmer(k=2, normalize=True, upto=True)
    vec = revckmer.make_revckmer_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    print "The vector is ", vec
    print

    print 'Begin IDkmer.'
    from dnavec.nac.kmer import IDkmer

    idkmer = IDkmer()
    vec = idkmer.make_idkmer_vector(open('test.txt'), open('hs.txt'), open('non-hs.txt'))
    print vec

    idkmer = IDkmer(k=2)
    vec = idkmer.make_idkmer_vector(open('test.txt'), open('hs.txt'), open('non-hs.txt'))
    print vec