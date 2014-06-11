__author__ = 'aleeee'
import sys

sys.path.append("..")
import util

import iMcRNA
import time


# The file path that we generate.
BRACKET_FILE_PATH = 'bracket.txt'
MATCHED_FILE_PATH = 'matched.txt'
P_VALUES_PATH = "p_values.txt"
N_GRAM_PATH = "N_gram.txt"
VECTOR_PATH = "vector.txt"

# The file path that we have.
N_GRAM_PERL_PATH = "64_coding_triplet_pri-sequence.pl"

if __name__ == '__main__':
    start_time = time.time()
    test_file = sys.argv[1]
    way_choice = sys.argv[2]
    if util.read_fasta_seq(open(test_file, 'r'), 'ACGU'):
        iMcRNA.generate_bracket_seq(test_file, BRACKET_FILE_PATH)
        iMcRNA.match_2st_has_name(BRACKET_FILE_PATH, MATCHED_FILE_PATH)
        p_values, n_gram, mfe = [], [], []
        if way_choice == 'ExSSC':
            p_values = iMcRNA.generate_p_values(test_file, P_VALUES_PATH)
            n_gram = iMcRNA.generate_n_gram(test_file, N_GRAM_PERL_PATH, N_GRAM_PATH)
            mfe = iMcRNA.generate_mfe(BRACKET_FILE_PATH)
        elif way_choice != 'SSC':
            print 'The way_choice error!'
        k, lamada, w, label = 1, 1, 0.1, '+1'
        iMcRNA.make_test_vector_libsvm_format(MATCHED_FILE_PATH, VECTOR_PATH, k, lamada, w, way_choice, p_values,
                                              n_gram, mfe, label)

    print time.time() - start_time, 's'