__author__ = 'aleeee'

import sys
import cPickle
from math import pow

from bioHelper.util import frequency
from bioHelper.kmer.util_kmer import make_kmer_list


def get_phyche_factor_dic(k):
    if 2 == k:
        f = open('mmc3.data', 'rb')
    elif 3 == k:
        f = open('mmc4.data', 'rb')
    else:
        sys.stderr.write("The k can just be 2 or 3.")
        sys.exit(0)

    phyche_factor_dic = cPickle.load(f)
    f.close()

    return phyche_factor_dic


def parallel_cor_function(k, nucleotide1, nucleotide2, phyche_list):
    """Get the parallel correlation Factor(Type 1)."""
    phyche_factor_dic = get_phyche_factor_dic(k)
    sum = 0.0

    for phyche in phyche_list:
        h1 = None
        h2 = None

        for phyche_tripe in phyche_factor_dic[nucleotide1]:
            if phyche_tripe[0] == phyche:
                h1 = phyche_tripe[1]
                print nucleotide1, phyche_tripe
        if None == h1:
            error_info = "Do not find the Physicochemical properties " + str(phyche)
            sys.stderr.write(error_info)
            sys.exit(0)

        for phyche_tripe in phyche_factor_dic[nucleotide2]:
            if phyche_tripe[0] == phyche:
                h2 = phyche_tripe[1]
                print nucleotide2, phyche_tripe
        if None == h2:
            error_info = "Do not find the Physicochemical properties " + str(phyche)
            sys.stderr.write(error_info)
            sys.exit(0)

        sum += pow(h1-h2, 2)

    return sum / len(phyche_list)


def get_cor_factor1(k, lamada, sequence, phyche_list):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)
    for i in range(1, lamada+1):
        sum = 0.0
        for j in range(0, l-k-i+1):
            nucleotide1 = ''
            nucleotide2 = ''
            for u in range(0, k):
                nucleotide1 += sequence[j+u]
                nucleotide2 += sequence[j+i+u]
            print nucleotide1, nucleotide2
            sum += parallel_cor_function(k, nucleotide1, nucleotide2, phyche_list)
            print sum

        theta.append(sum/(l-k-i+1))

    return theta


def make_type1_vector(sequence_list, k, lamada, w, alphabet, phyche_list):
    """Generate the parallel(Type 1) vector."""
    vector = []
    kmer = make_kmer_list(k, alphabet)

    for sequence in sequence_list:
        # Get the dinucleotide frequency in the DNA sequence.
        fre_list = []
        fre_sum = 0.0
        # print '-----------------------------'
        for key in kmer:
            temp_fre = frequency(sequence, str(key))
            fre_list.append(temp_fre)
            fre_sum += temp_fre
            # print key, temp_fre
            if temp_fre != 0:
                print key, temp_fre
        # Get the normalized occurrence frequency of dinucleotide in the DNA sequence.
        fre_len = len(fre_list)
        for i in range(0, fre_len):
            fre_list[i] /= fre_sum

        # Get the theta_list according the Equation 5.
        theta_list = get_cor_factor1(k, lamada, sequence, phyche_list)

        # Generate the vector according the Equation 9.
        theta_sum = 0.0
        for theta in theta_list:
            theta_sum += theta
        print theta_sum
        denominator = 1 + w*theta_sum

        temp_vec = []
        for f in fre_list:
            # print f
            temp_vec.append(round(f/denominator, 3))
        print denominator
        # print theta_list
        for theta in theta_list:
            temp_vec.append(round(w*theta/denominator, 3))

        vector.append(temp_vec)

    return vector


def series_cor_function(k, nucleotide1, nucleotide2, phyche):
    """Get the series correlation Factor(Type 2)."""
    phyche_factor_dic = get_phyche_factor_dic(k)

    h1 = 0
    h2 = 0
    for phyche_tripe in phyche_factor_dic[nucleotide1]:
        # print phyche_tripe, phyche
        if phyche_tripe[0] == phyche:
            h1 = phyche_tripe[1]
    for phyche_tripe in phyche_factor_dic[nucleotide2]:
        if phyche_tripe[0] == phyche:
            h2 = phyche_tripe[1]

    return h1*h2


def get_cor_factor2(k, lamada, sequence, phyche_list):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)
    for lamada1 in range(1, lamada+1):
        for phyche in phyche_list:
            sum = 0.0
            for i in range(0, l-k-lamada1):
                nucleotide1 = ''
                nucleotide2 = ''
                for u in range(0, k):
                    nucleotide1 += sequence[i+u]
                    nucleotide2 += sequence[i+lamada1+u]
                # print nucleotide1, nucleotide2
                sum += series_cor_function(k, nucleotide1, nucleotide2, phyche)

            theta.append(sum/(l-k-lamada1))

    return theta


def make_type2_vector(sequence_list, k, lamada, w, alphabet, phyche_list):
    """Generate the series(Type 2) vector."""
    vector = []
    kmer = make_kmer_list(k, alphabet)

    for sequence in sequence_list:
        # Get the nucleotide frequency in the DNA sequence.
        fre_list = []
        fre_sum = 0.0
        for key in kmer:
            temp_fre = frequency(sequence, str(key))
            fre_list.append(temp_fre)
            fre_sum += temp_fre

        # Get the normalized occurrence frequency of dinucleotide in the DNA sequence.
        fre_len = len(fre_list)
        for i in range(0, fre_len):
            fre_list[i] /= fre_sum

        # Get the theta_list according the Equation 13.
        theta_list = get_cor_factor2(k, lamada, sequence, phyche_list)

        # Generate the vector according the Equation 16.
        theta_sum = 0.0
        for theta in theta_list:
            theta_sum += theta
        denominator = 1 + w*theta_sum

        temp_vec = []
        for f in fre_list:
            temp_vec.append(round(f/denominator, 3))
        for theta in theta_list:
            temp_vec.append(round(w*theta/denominator, 3))

        vector.append(temp_vec)

    return vector