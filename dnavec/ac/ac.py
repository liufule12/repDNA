__author__ = 'liufule12'

from dnavec.util import get_data, generate_phyche_value


class DAC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2

    def make_ac_vector(self, input_data, phyche_list=[], all_property=False, extra_phyche_index={}):
        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_list, all_property, extra_phyche_index)
        print phyche_value

        from dnavec.ac.acutil import make_ac_vector

        return make_ac_vector(sequence_list, self.lag, phyche_value, self.k)


class DCC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2

    def make_cc_vector(self, input_data, phyche_list=[], all_property=False, extra_phyche_index={}):
        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_list, all_property, extra_phyche_index)
        print phyche_value

        from dnavec.ac.acutil import make_cc_vector

        return make_cc_vector(sequence_list, self.lag, phyche_value, self.k)


class DACC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2

    def make_acc_vector(self, input_data, phyche_list=[], all_property=False, extra_phyche_index={}):
        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_list, all_property, extra_phyche_index)
        print phyche_value

        from dnavec.ac.acutil import make_ac_vector, make_cc_vector

        zipped = zip(make_ac_vector(sequence_list, self.lag, phyche_value, self.k),
                     make_cc_vector(sequence_list, self.lag, phyche_value, self.k))
        vector = [reduce(lambda x, y: x + y, e) for e in zipped]

        return vector


class TAC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3

    def make_tac_vector(self, input_data, phyche_list=[], all_property=False, extra_phyche_index={}):
        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_list, all_property, extra_phyche_index)
        print phyche_value

        from dnavec.ac.acutil import make_ac_vector

        return make_ac_vector(sequence_list, self.lag, phyche_value, self.k)


class TCC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3

    def make_cc_vector(self, input_data, phyche_list=[], all_property=False, extra_phyche_index={}):
        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_list, all_property, extra_phyche_index)
        print phyche_value

        from dnavec.ac.acutil import make_cc_vector

        return make_cc_vector(sequence_list, self.lag, phyche_value, self.k)


class TACC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3

    def make_acc_vector(self, input_data, phyche_list=[], all_property=False, extra_phyche_index={}):
        sequence_list = get_data(input_data)

        phyche_value = generate_phyche_value(self.k, phyche_list, all_property, extra_phyche_index)
        print phyche_value

        from dnavec.ac.acutil import make_ac_vector, make_cc_vector

        zipped = zip(make_ac_vector(sequence_list, self.lag, phyche_value, self.k),
                     make_cc_vector(sequence_list, self.lag, phyche_value, self.k))
        vector = [reduce(lambda x, y: x + y, e) for e in zipped]

        return vector


if __name__ == '__main__':
    from pprint import pprint

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

    ac = DAC(1)
    vec = ac.make_ac_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], [], all_property=True)
    print(vec)
    print len(vec[0])

    cc = DCC(2)
    vec = cc.make_cc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], [],
                            all_property=True)
    print(vec)
    print len(vec[0]), len(vec[1])

    acc = DACC(1)
    vec = acc.make_acc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], [],
                              all_property=True)
    for e in vec:
        print e
    print len(vec[0]), len(vec[1])

    acc = DACC(2)
    vec = acc.make_acc_vector(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], [],
                              all_property=True, extra_phyche_index=extra_phyche_value)
    for e in vec:
        print e
    print len(vec[0]), len(vec[1])