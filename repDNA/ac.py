__author__ = 'liufule12'

from repDNA.util import get_data, generate_phyche_value


class DAC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2

    def make_dac_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        if phyche_index is None:
            phyche_index = []
        if extra_phyche_index is None:
            extra_phyche_index = {}

        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_index, all_property, extra_phyche_index)
        # print phyche_value

        from repDNA.acutil import make_ac_vector

        return make_ac_vector(sequence_list, self.lag, phyche_value, self.k)


class DCC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2

    def make_dcc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        if phyche_index is None:
            phyche_index = []
        if extra_phyche_index is None:
            extra_phyche_index = {}

        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_index, all_property, extra_phyche_index)

        from repDNA.acutil import make_cc_vector

        return make_cc_vector(sequence_list, self.lag, phyche_value, self.k)


class DACC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2

    def make_dacc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        if phyche_index is None:
            phyche_index = []
        if extra_phyche_index is None:
            extra_phyche_index = {}

        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_index, all_property, extra_phyche_index)
        # print phyche_value

        from repDNA.acutil import make_ac_vector, make_cc_vector

        zipped = zip(make_ac_vector(sequence_list, self.lag, phyche_value, self.k),
                     make_cc_vector(sequence_list, self.lag, phyche_value, self.k))
        vector = [reduce(lambda x, y: x + y, e) for e in zipped]

        return vector


class TAC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3

    def make_tac_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        if phyche_index is None:
            phyche_index = []
        if extra_phyche_index is None:
            extra_phyche_index = {}

        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_index, all_property, extra_phyche_index)
        # print phyche_value

        from repDNA.acutil import make_ac_vector

        return make_ac_vector(sequence_list, self.lag, phyche_value, self.k)


class TCC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3

    def make_tcc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        if phyche_index is None:
            phyche_index = []
        if extra_phyche_index is None:
            extra_phyche_index = {}

        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_index, all_property, extra_phyche_index)
        # print phyche_value

        from repDNA.acutil import make_cc_vector

        return make_cc_vector(sequence_list, self.lag, phyche_value, self.k)


class TACC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3

    def make_tacc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        if phyche_index is None:
            phyche_index = []
        if extra_phyche_index is None:
            extra_phyche_index = {}

        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_index, all_property, extra_phyche_index)
        # print phyche_value

        from repDNA.acutil import make_ac_vector, make_cc_vector

        zipped = zip(make_ac_vector(sequence_list, self.lag, phyche_value, self.k),
                     make_cc_vector(sequence_list, self.lag, phyche_value, self.k))
        vector = [reduce(lambda x, y: x + y, e) for e in zipped]

        return vector


if __name__ == '__main__':
    extra_phyche_value = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
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
    # phyche_index = \
    # [[0.026, 0.036, 0.031, 0.033, 0.016, 0.026, 0.014, 0.031, 0.025, 0.025, 0.026, 0.036, 0.017, 0.025, 0.016,
    # 0.026],
    # [0.038, 0.038, 0.037, 0.036, 0.025, 0.042, 0.026, 0.037, 0.038, 0.036, 0.042, 0.038, 0.018, 0.038, 0.025,
    # 0.038]]

    phyche_index = \
        [[2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
         [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]

    from repDNA.util import normalize_index

    dac = DAC(2)
    vec = dac.make_dac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'])
    print vec
    print len(vec[0])

    vec = dac.make_dac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    print(vec)
    print len(vec[0])

    # print normalize_index(phyche_index)

    vec = dac.make_dac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                           extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    print vec
    print len(vec[0])
    print

    dcc = DCC(2)
    vec = dcc.make_dcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'])
    print(vec)
    print len(vec[0])

    vec = dcc.make_dcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    print(vec)
    print len(vec[0])

    vec = dcc.make_dcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                           extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    print vec
    print len(vec[0])
    print

    print 'DACC'
    dacc = DACC(2)
    vec = dacc.make_dacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'])
    print vec
    print len(vec[0])

    vec = dacc.make_dacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    print vec
    print len(vec[0])

    vec = dacc.make_dacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                             extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    print vec
    print len(vec[0])
    print

    phyche_index = [
        # [6.882, 5.260, 3.995, 6.698, 3.516, 3.619, 3.625, 4.471, 3.879, 2.683, 3.782, 4.471, 3.047, 4.153, 2.185, 6.698,
        #  3.958, 2.832, 2.671, 2.185, 5.000, 3.311, 4.502, 3.782, 2.570, 3.275, 4.502, 3.625, 3.813, 3.221, 2.671, 3.995,
        #  4.385, 3.498, 3.221, 4.153, 2.754, 1.387, 3.275, 2.683, 3.819, 1.387, 3.311, 3.619, 3.770, 3.498, 2.832, 5.260,
        #  4.013, 3.770, 3.813, 3.047, 2.197, 3.819, 2.570, 3.879, 10.000, 2.754, 5.000, 3.516, 4.013, 4.385, 3.958,
        #  0.100],
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581, 7.176]
    ]
    print normalize_index(phyche_index, is_convert_dict=True)

    print 'Begin TAC'
    tac = TAC(2)
    vec = tac.make_tac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'])
    print vec
    print len(vec[0])

    vec = tac.make_tac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    print vec
    print len(vec[0])

    vec = tac.make_tac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                           extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    print vec
    print len(vec[0])
    print

    print 'Begin TCC'
    tcc = TCC(2)
    vec = tcc.make_tcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'])
    print vec
    print len(vec[0])

    vec = tcc.make_tcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    print vec
    print len(vec[0])

    vec = tcc.make_tcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                           extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    print vec
    print len(vec[0])
    print

    print 'Bengin TACC'
    tacc = TACC(2)
    vec = tacc.make_tacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'])
    print vec
    print len(vec[0])

    vec = tacc.make_tacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    print vec
    print len(vec[0])

    vec = tacc.make_tacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                             extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    print vec
    print len(vec[0])
    print