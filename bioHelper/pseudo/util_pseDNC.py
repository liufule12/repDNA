__author__ = 'aleeee'

from math import pow
import sys
sys.path.append('..')

from util import frequency


U = 6
pu_ri_rj = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
            'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
            'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
            'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17],
            'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
            'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
            'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39],
            'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
            'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
            'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59],
            'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
            'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
            'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39],
            'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
            'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
            'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11]}


def cor_function(dinuleotide1, dinuleotide2):
    """Get the cFactor."""
    sum = 0.0
    for u in range(0, U):
        sum += pow(pu_ri_rj[dinuleotide1][u] - pu_ri_rj[dinuleotide2][u], 2)
    return sum / U


def get_cor_factor(lamada, sequence):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)
    for i in range(1, lamada+1):
        sum = 0.0
        for j in range(0, l-1-lamada):
            rirj1 = sequence[j] + sequence[j+1]
            rirj2 = sequence[j+i] + sequence[j+i+1]
            sum += cor_function(rirj1, rirj2)

        theta.append(sum/(l-i-1))

    return theta


def make_psednc_vector(sequence_list, lamada, w):
    """Generate the psednc vector."""
    vector = []
    for sequence in sequence_list:

        # Get the dinucleotide frequency in the DNA sequence.
        fre_list = []
        fre_sum = 0.0
        for key in pu_ri_rj.keys():
            temp_fre = frequency(sequence, str(key))
            fre_list.append(temp_fre)
            fre_sum += temp_fre

        # Get the normalized occurrence frequency of dinucleotide in the DNA sequence.
        fre_len = len(fre_list)
        for i in range(0, fre_len):
            fre_list[i] /= fre_sum

        # Get the theta_list according the Equation 5.
        theta_list = get_cor_factor(lamada, sequence)

        # Generate the vector according the Equation 9.
        theta_sum = 0.0
        for theta in theta_list:
            theta_sum += theta
        denominator = 1 + w*theta_sum

        temp_vec = []
        for f in fre_list:
            temp_vec.append(f/denominator)
        for theta in theta_list:
            temp_vec.append(w*theta/denominator)

        vector.append(temp_vec)
        print len(temp_vec)

    return vector








