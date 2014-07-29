__author__ = 'aleeee'

from dnavec.util import get_data
from dnavec.psenac.psenacutil import extend_phyche_index


class PseDNC():
    def __init__(self, lamada=3, w=0.05):
        self.lamada = lamada
        self.w = w

    def make_psednc_vector(self, input_data, extra_phyche_value={}):
        from dnavec.psenac.psenacutil import make_pseknc_vector

        original_phyche_value = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
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

        sequence_list = get_data(input_data)
        phyche_value = extend_phyche_index(original_phyche_value, extra_phyche_value)
        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, 2, phyche_value, theta_type=1)

        return vector


class PseKNC():
    """This class should be used to make PseKNC vector."""

    def __init__(self, k=3, lamada=1, w=0.5):
        """
        :param k: k-tuple.
        """
        self.k = k
        self.lamada = lamada
        self.w = w

    def make_pseknc_vector(self, input_data, extra_phyche_value={}):
        """Make PseKNC vector.

        :param input_data: The fasta file path or single DNA sequence or DNA sequence list.
        :return: vector: The iNucPseKNC vector.
        """
        from dnavec.psenac.psenacutil import make_old_pseknc_vector

        original_phyche_value = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
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

        sequence_list = get_data(input_data)
        phyche_value = extend_phyche_index(original_phyche_value, extra_phyche_value)
        return make_old_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)


class PCPseDNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 2

    def make_pcpsednc_vector(self, input_data, phyche_list, all_property=False, extra_phyche_index={}):
        """Make a PseDNC type1 vector.

        :param input_data: file object or sequence list.
        :param phyche_list: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :return: PseDNC type1 vector.
        """
        diphyche_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content',
                         'A-philicity', 'Propeller twist', 'Duplex stability:(freeenergy)',
                         'Duplex tability(disruptenergy)',
                         'DNA denaturation', 'Bending stiffness', 'Protein DNA twist', 'Stabilising energy of Z-DNA',
                         'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Breslauer_dS', 'Electron_interaction',
                         'Hartman_trans_free_energy', 'Helix-Coil_transition', 'Ivanov_BA_transition',
                         'Lisser_BZ_transition',
                         'Polar_interaction', 'SantaLucia_dG', 'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility',
                         'Stability', 'Stacking_energy', 'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS',
                         'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise']

        sequence_list = get_data(input_data)

        # Set and check physicochemical properties.
        if all_property is True:
            phyche_list = diphyche_list
        else:
            for e in phyche_list:
                if e not in diphyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    import sys
                    sys.stderr.write(error_info)
                    sys.exit(0)

        # Generate phyche_value.
        from dnavec.psenac.psenacutil import make_pseknc_vector
        from dnavec.psenac.psenacutil import get_phyche_index

        phyche_value = extend_phyche_index(get_phyche_index(self.k, phyche_list), extra_phyche_index)
        # print phyche_value
        # Make vector.
        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)

        return vector


class PCPseTNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 3

    def make_pcpsetnc_vector(self, input_data, phyche_list, all_property=False, extra_phyche_index={}):
        """Make a PseDNC type1 vector.

        :param input_data: file object or sequence list.
        :param phyche_list: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :return: PseDNC type1 vector.
        """
        triphyche_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                          'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                          'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

        sequence_list = get_data(input_data)

        # Set and check physicochemical properties.
        if all_property is True:
            phyche_list = triphyche_list
        else:
            for e in phyche_list:
                if e not in triphyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    import sys

                    sys.stderr.write(error_info)
                    sys.exit(0)

        # Generate phyche_value.
        from dnavec.psenac.psenacutil import make_pseknc_vector
        from dnavec.psenac.psenacutil import get_phyche_index

        phyche_value = extend_phyche_index(get_phyche_index(self.k, phyche_list), extra_phyche_index)
        print phyche_value
        # Make vector.
        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)

        return vector


class SCPseDNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 2

    def make_scpsednc_vector(self, input_data, phyche_list, all_property=False, extra_phyche_index={}):
        """Make a PseDNC type2 vector.

        :param input_data: file object or sequence list.
        :param phyche_list: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :return: PseDNC type2 vector.
        """
        diphyche_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content',
                         'A-philicity', 'Propeller twist', 'Duplex stability:(freeenergy)',
                         'Duplex tability(disruptenergy)',
                         'DNA denaturation', 'Bending stiffness', 'Protein DNA twist', 'Stabilising energy of Z-DNA',
                         'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Breslauer_dS', 'Electron_interaction',
                         'Hartman_trans_free_energy', 'Helix-Coil_transition', 'Ivanov_BA_transition',
                         'Lisser_BZ_transition',
                         'Polar_interaction', 'SantaLucia_dG', 'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility',
                         'Stability', 'Stacking_energy', 'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS',
                         'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise']

        sequence_list = get_data(input_data)

        # Set and check physicochemical properties.
        if all_property is True:
            phyche_list = diphyche_list
        else:
            for e in phyche_list:
                if e not in diphyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    import sys

                    sys.stderr.write(error_info)
                    sys.exit(0)

        # Generate phyche_value.
        from dnavec.psenac.psenacutil import make_pseknc_vector
        from dnavec.psenac.psenacutil import get_phyche_index

        phyche_value = extend_phyche_index(get_phyche_index(self.k, phyche_list), extra_phyche_index)
        print phyche_value
        # Make vector.
        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=2)

        return vector


class SCPseTNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 3

    def make_scpsetnc_vector(self, input_data, phyche_list, all_property=False, extra_phyche_index={}):
        """Make a PseDNC type2 vector.

        :param input_data: file object or sequence list.
        :param phyche_list: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :return: PseDNC type2 vector.
        """
        triphyche_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                          'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                          'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

        sequence_list = get_data(input_data)

        # Set and check physicochemical properties.
        if all_property is True:
            phyche_list = triphyche_list
        else:
            for e in phyche_list:
                if e not in triphyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    import sys

                    sys.stderr.write(error_info)
                    sys.exit(0)

        # Generate phyche_value.
        from dnavec.psenac.psenacutil import make_pseknc_vector
        from dnavec.psenac.psenacutil import get_phyche_index

        phyche_value = extend_phyche_index(get_phyche_index(self.k, phyche_list), extra_phyche_index)
        # print phyche_value
        # Make vector.
        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=2)

        return vector


if __name__ == '__main__':
    import time

    start_time = time.time()

    psednc = PseDNC(lamada=1, w=0.05)
    res = psednc.make_psednc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    for e in res:
        print e
    print len(e)
    print '-----------------------------------------------------------------'

    extra_phyche_index = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1],
    'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
    'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
                          'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17, 1],
                          'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
                          'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
                          'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39, 1],
                          'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
                          'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
                          'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59, 1],
                          'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
                          'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
                          'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39, 1],
                          'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
                          'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
                          'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1]}

    extra_phyche_index = {'AA': [1.02, -0.64, 111], 'AC': [-0.92, -0.83, 1], 'AG': [0.49, -0.89, 1],
                          'AT': [0.57, -1.05, 1],
                          'CA': [0.57, 1.51, 1], 'CC': [-0.07, 0.36, 1], 'CG': [-0.58, 2.23, 1], 'CT': [0.49, -0.89, 1],
                          'GA': [-0.65, -0.14, 1], 'GC': [-2.46, -0.30, 1], 'GG': [-0.07, 0.36, 1],
                          'GT': [-0.92, -0.83, 1],
                          'TA': [1.60, 0.42, 1], 'TC': [-0.65, -0.14, 1], 'TG': [0.57, 1.51, 1], 'TT': [1.02, -0.64, 1]}

    psednc = PseDNC(lamada=1, w=0.05)
    vector = psednc.make_psednc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], extra_phyche_index)
    for e in vector:
        print e
        print len(e)

    pseknc = PseKNC()
    res = pseknc.make_pseknc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    for e in res:
        print e
        print len(e)

    pseknc = PseKNC(k=2, lamada=1, w=0.05)
    res = pseknc.make_pseknc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    for e in res:
        print e
        print len(e)

    pseknc = PseKNC(k=2, lamada=1, w=0.05)
    res = pseknc.make_pseknc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], extra_phyche_index)
    for e in res:
        print e
        print len(e)

    print 'Time:', time.time() - start_time

    print

    extra_phyche_index = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1],
                          'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
                          'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
                          'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17, 1],
                          'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
                          'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
                          'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39, 1],
                          'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
                          'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
                          'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59, 1],
                          'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
                          'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
                          'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39, 1],
                          'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
                          'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
                          'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1]}

    pc_psednc = PCPseDNC(1, 0.05)
    print 'This is PC-PseDNC', pc_psednc.make_pcpsednc_vector(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'],
                                                              phyche_list=['Base stacking'],
                                                              all_property=False)
    print 'This is extra PC-PseDNC', pc_psednc.make_pcpsednc_vector(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'],
                                                                    phyche_list=['Base stacking'], all_property=False,
                                                                    extra_phyche_index=extra_phyche_index)

    pc_psetnc = PCPseTNC(1, 0.05)
    print 'This is PC-PseTNC'
    print pc_psetnc.make_pcpsetnc_vector(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Dnase I'])

    sc_psednc = SCPseDNC(lamada=1, w=0.05)
    print 'This is SC-PseDNC'
    print sc_psednc.make_scpsednc_vector(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Base stacking'])

    sc_psetnc = SCPseTNC(lamada=1, w=0.05)
    print 'This is SC-PseTNC'
    print sc_psetnc.make_scpsetnc_vector(input_data=['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_list=['Dnase I'])

