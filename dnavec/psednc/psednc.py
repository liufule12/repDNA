__author__ = 'aleeee'

from dnavec.util import get_data
from dnavec.psednc.psedncutil import make_psednc_vector


class PseDNC():
    def __init__(self, lamada, w, alphabet="ACGT"):
        self.lamada = lamada
        self.w = w
        self.alphabet = alphabet

    def make_vector(self, data):
        sequence_list = get_data(data)

        vector = make_psednc_vector(sequence_list, self.lamada, self.w)

        return vector


if __name__ == '__main__':
    import time
    start_time = time.time()
    psednc = PseDNC(1, 0.05)
    vector = psednc.make_vector(open('hs.txt'))
    for e in vector:
        print e
    print 'Time:', time.time() - start_time
