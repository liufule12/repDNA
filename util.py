__author__ = 'aleeee'

import sys


class Seq:
    def __init__(self, name, seq, no):
        self.name = name
        self.seq = seq.upper()
        self.length = len(seq)
        self.no = no

    def __str__(self):
        """Output seq when 'print' method is called."""
        return ">%s\tNo:%s\tlength:%s\n%s" % (self.name, str(self.no), str(self.length), self.seq)


def is_format_legal(seq, alphabet):
    """
    Judge the seq is in alphabet or null.
    Four situation:
    1. No seq name.
    2. Seq name is illegal.
    3. No sequence.
    4. Sequence is illegal.
    """
    if not seq.name:
        print 'Error, sequence', seq.no, 'has no sequence name.'
        return False
    if -1 != seq.name.find('>'):
        print 'Error, sequence', seq.no, 'name has > character.'
        return Fasta
    if 0 == seq.length:
        print 'Error, sequence', seq.no, 'is null.'
        return False
    for e in seq.seq:
        if e not in alphabet:
            print 'Error, sequence', seq.no, 'has non-alphabet character.'
            return False

    return True


def read_seq(f):
    """
    Yields a Seq object.
    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    """
    name, seq = '', ''
    count = 0
    while True:
        line = f.readline()
        if not line:
            break

        if '>' == line[0]:
            if 0 != count or (0 == count and seq != ''):
                yield Seq(name, seq, count)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    count += 1
    yield Seq(name, seq, count)


def judge_fasta(f, alphabet_choice):
    """
    Judge the input is Fasta or not.
    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    Return True or False.
    """
    if alphabet_choice == 'rna':
        alphabet = ['A', 'C', 'G', 'U']
    elif alphabet_choice == 'amino':
        alphabet = ['A', 'R', 'D', 'C', 'Q', 'E', 'H', 'I', 'G', 'N', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    for e in read_seq(f):
        # print e
        if not is_format_legal(e, alphabet):
            return False
    return True

if __name__ == '__main__':
    test_file = sys.argv[1]
    judge_fasta(open(test_file, 'r'), 'rna')