__author__ = 'aleeee'

import sys
from pseDNC import PseDNC

# k, upto, revcomp, normalize = 1, False, False, False
# if len(sys.argv[1]) != 0:
#     k = int(sys.argv[1])
# if sys.argv[2] == '-upto':
#     upto = True
# if sys.argv[3] == '-revcomp':
#     revcomp = True
# if sys.argv[4] == '-normalize':
#     normalize = True

psednc = PseDNC(1, 0.3)
vector = psednc.make_vector('hs.txt')
for e in vector:
    print e

