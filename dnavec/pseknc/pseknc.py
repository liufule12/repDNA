__author__ = 'aleeee'

import cPickle
import sys

from dnavec.util import get_data, write_libsvm

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

    def make_vector(self, data, vec_type, phyche_list, added='none'):
        """Make a PseKNC vector.

        :param data: file object or sequence list.
        :param vec_type: 1 or 2.
        :param phyche_list: physicochemical properties list.
        :param added: choose all physicochemical properties or not.
        :return: PseKNC vector.

        Examples
        --------

        >>> from dnavec.pseknc.pseknc import PseKNC
        >>> pseKNC = PseKNC(k=2, lamada=1, w=0.05)
        >>> pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], vec_type=1, phyche_list=['Twist', 'Tilt'])
        [[0.027, 0.08, 0.0, 0.08, 0.053, 0.0, 0.0, 0.106, 0.053, 0.053, 0.027, 0.027, 0.053, 0.053, 0.106, 0.186, 0.095]]

        >>> pseKNC.make_vector(data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], vec_type=1, phyche_list=[], added='all')
        [[0.027, 0.081, 0.0, 0.081, 0.054, 0.0, 0.0, 0.108, 0.054, 0.054, 0.027, 0.027, 0.054, 0.054, 0.108, 0.19, 0.078]]

        """
        sequence_list = get_data(data)

        if self.k == 2:
            if added == 'all':
                phyche_list = DIPYHYCHE_LIST
            else:
                for e in phyche_list:
                    if e not in DIPYHYCHE_LIST:
                        error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                        sys.stderr.write(error_info)
                        sys.exit(0)
        elif self.k == 3:
            if added == 'all':
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

        if vec_type == 1:
            from dnavec.pseknc.psekncutil import make_type1_vector
            vector = make_type1_vector(sequence_list, self.k, self.lamada, self.w, self.alphabet, phyche_list)
        elif vec_type == 2:
            from dnavec.pseknc.psekncutil import make_type2_vector
            vector = make_type2_vector(sequence_list, self.k, self.lamada, self.w, self.alphabet, phyche_list)

        return vector


if __name__ == '__main__':
    k = 3
    lamada = 1
    w = 0.05
    vec_type = 1
    pseKNC = PseKNC(k, lamada, w)
    phyche_list = ['Bendability (DNAse)']
    vector = pseKNC.make_vector(open('hs.txt'), vec_type, phyche_list)
    for e in vector:
        print e
        print len(e)
    print len(vector)