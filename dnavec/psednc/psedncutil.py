__author__ = 'aleeee'

import sys
from math import pow

from dnavec.util import frequency


U = 6
PU_RI_RJ = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
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
    for u in range(U):
        sum += pow(PU_RI_RJ[dinuleotide1][u] - PU_RI_RJ[dinuleotide2][u], 2)
    return sum / U


def get_cor_factor(lamada, sequence):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)
    for i in range(1, lamada+1):
        sum = 0.0
        for j in range(0, l-1-lamada):
            nucleotide1 = sequence[j] + sequence[j+1]
            nucleotide2 = sequence[j+i] + sequence[j+i+1]
            sum += cor_function(nucleotide1, nucleotide2)

        theta.append(sum/(l-i-1))

    return theta


def make_psednc_vector(sequence_list, lamada, w):
    """Generate the pseknc vector."""
    vector = []

    for sequence in sequence_list:
        if len(sequence) < 2+lamada:
            error_info = "Sorry, the sequence length must be larger than " + str(lamada+2)
            sys.stderr.write(error_info)
            sys.exit(0)

        # Get the dinucleotide frequency in the DNA sequence.
        fre_list = [frequency(sequence, str(key)) for key in PU_RI_RJ.keys()]
        fre_sum = float(sum(fre_list))

        # Get the normalized occurrence frequency of dinucleotide in the DNA sequence.
        fre_list = [e / fre_sum for e in fre_list]

        # Get the theta_list according the Equation 5.
        theta_list = get_cor_factor(lamada, sequence)
        theta_sum = sum(theta_list)

        # Generate the vector according the Equation 9.
        denominator = 1 + w*theta_sum

        temp_vec = [round(f/denominator, 3) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w*theta/denominator, 4))

        vector.append(temp_vec)

    return vector