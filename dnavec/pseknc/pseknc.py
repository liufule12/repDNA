__author__ = 'aleeee'

import sys

from dnavec.util import get_data

DIPYHYCHE_LIST = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content',
                  'A-philicity', 'Propeller twist', 'Duplex stability:(freeenergy)', 'Duplex tability(disruptenergy)',
                  'DNA denaturation', 'Bending stiffness', 'Protein DNA twist', 'Stabilising energy of Z-DNA',
                  'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Breslauer_dS', 'Electron_interaction',
                  'Hartman_trans_free_energy', 'Helix-Coil_transition', 'Ivanov_BA_transition', 'Lisser_BZ_transition',
                  'Polar_interaction', 'SantaLucia_dG', 'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility',
                  'Stability', 'Stacking_energy', 'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS',
                  'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise']
TRIPHYCHE_LIST = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                  'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                  'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']


class PseKNC():
    def __init__(self, k, lamada, w, alphabet='ACGT'):
        self.k = k
        self.lamada = lamada
        self.w = w
        self.alphabet = alphabet

    def make_vector(self, data, phyche_list, vec_type=1, all_property=False):
        """Make a PseKNC vector.

        :param data: file object or sequence list.
        :param vec_type: 1 or 2.
        :param phyche_list: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :return: PseKNC vector.

        Examples
        --------

        >>> from dnavec.pseknc.pseknc import PseKNC

        >>> pseKNC = PseKNC(k=2, lamada=1, w=0.05)
        >>> pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Twist', 'Tilt'])
        [[0.027, 0.081, 0.0, 0.081, 0.054, 0.0, 0.0, 0.108, 0.054, 0.054, 0.027, 0.027, 0.054, 0.054, 0.108, 0.19, 0.078]]

        >>> pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], vec_type=2, phyche_list=[], all_property=True)
        [[0.026, 0.079, 0.0, 0.079, 0.052, 0.0, 0.0, 0.105, 0.052, 0.052, 0.026, 0.026, 0.052, 0.052, 0.105, 0.183, 0.001,
         -0.005, -0.005, 0.023, 0.008, 0.016, 0.017, 0.004, 0.001, 0.017, -0.015, -0.01, -0.004, 0.004, -0.008, -0.009,
         0.023, -0.015, 0.004, 0.008, -0.015, 0.023, 0.017, -0.006, -0.007, -0.003, 0.006, 0.001, 0.017, 0.008, 0.004,
         0.023, -0.008, -0.009, 0.036, 0.009, -0.022, -0.022]]

        >>> pseKNC = PseKNC(k=3, lamada=1, w=0.05)
        >>> pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Dnase I'])
        [[0.0, 0.027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.082, 0.0, 0.0, 0.0, 0.0, 0.027, 0.0, 0.0, 0.054, 0.0, 0.027, 0.0, 0.027,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.027, 0.054, 0.027, 0.027, 0.027, 0.0, 0.0, 0.027, 0.0, 0.0, 0.027,
        0.0, 0.0, 0.0, 0.027, 0.0, 0.0, 0.0, 0.027, 0.0, 0.0, 0.0, 0.054, 0.027, 0.0, 0.0, 0.0, 0.027, 0.054, 0.027, 0.0,
        0.027, 0.027, 0.054, 0.082, 0.102]]

        >>> pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], vec_type=2, phyche_list=['Dnase I'])
        [[0.0, 0.031, 0.0, 0.0, 0.0, 0.0, 0.0, 0.092, 0.0, 0.0, 0.0, 0.0, 0.031, 0.0, 0.0, 0.061, 0.0, 0.031, 0.0, 0.031,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.031, 0.061, 0.031, 0.031, 0.031, 0.0, 0.0, 0.031, 0.0, 0.0, 0.031,
        0.0, 0.0, 0.0, 0.031, 0.0, 0.0, 0.0, 0.031, 0.0, 0.0, 0.0, 0.061, 0.031, 0.0, 0.0, 0.0, 0.031, 0.061, 0.031, 0.0,
        0.031, 0.031, 0.061, 0.092, -0.009]]
        """
        sequence_list = get_data(data)

        # Set and check physicochemical properties.
        if self.k == 2:
            if all_property is True:
                phyche_list = DIPYHYCHE_LIST
            else:
                for e in phyche_list:
                    if e not in DIPYHYCHE_LIST:
                        error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                        sys.stderr.write(error_info)
                        sys.exit(0)
        elif self.k == 3:
            if all_property is True:
                phyche_list = TRIPHYCHE_LIST
            else:
                for e in phyche_list:
                    if e not in TRIPHYCHE_LIST:
                        error_info = 'Sorry, the physicochemical properties ' + e + 'is not exit.'
                        sys.stderr.write(error_info)
                        sys.exit(0)
        else:
            sys.stderr.write("The value of k must be 2 or 3.")
            sys.exit(0)

        # Set type and make vector.
        if vec_type == 1:
            from dnavec.pseknc.psekncutil import make_type1_vector
            vector = make_type1_vector(sequence_list, self.k, self.lamada, self.w, self.alphabet, phyche_list)
        elif vec_type == 2:
            from dnavec.pseknc.psekncutil import make_type2_vector
            vector = make_type2_vector(sequence_list, self.k, self.lamada, self.w, self.alphabet, phyche_list)
        else:
            error_info = 'Sorry, the type can only be 1 or 2, but your type is ' + str(vec_type) + '.'
            sys.stderr.write(error_info)
            sys.exit(0)

        return vector


if __name__ == '__main__':
    pseKNC = PseKNC(k=2, lamada=1, w=0.05)
    vector = pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Twist', 'Tilt'])
    for e in vector:
        print e
        print len(e)
    print len(vector)

    vector = pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], vec_type=2, phyche_list=[], all_property=True)
    for e in vector:
        print e
        print len(e)
    print len(vector)

    pseKNC = PseKNC(k=3, lamada=1, w=0.05)
    vector = pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Dnase I'])
    for e in vector:
        print e
        print len(e)
    print len(vector)

    vector = pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], vec_type=2, phyche_list=['Dnase I'])
    for e in vector:
        print e
        print len(e)
    print len(vector)