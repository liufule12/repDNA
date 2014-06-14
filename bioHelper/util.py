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


def frequency(tol_str, tar_str):
    """
    Generate the frequency of tar_str in tol_str.
    Return the count.
    """
    i, j, tar_count = 0, 0, 0
    len_tol_str = len(tol_str)
    len_tar_str = len(tar_str)
    while i < len_tol_str and j < len_tar_str:
        if tol_str[i] == tar_str[j]:
            i += 1
            j += 1
            if j >= len_tar_str:
                tar_count += 1
                i = i - j + 1
                j = 0
        else:
            i = i - j + 1
            j = 0

    return tar_count


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
        return False
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


def read_fasta_seq(f, alphabet):
    """
    Read the fasta file.
    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    Return the seq list.
    """
    seq = []
    for e in read_seq(f):
        # print e
        if not is_format_legal(e, alphabet):
            return None
        seq.append(e)

    return seq


def read_fasta_sequence(f, alphabet):
    """
    Read the fasta file.
    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    Return the seq list.
    """
    sequence_list = []
    for e in read_seq(f):
        # print e
        if not is_format_legal(e, alphabet):
            return None
        sequence_list.append(e.seq)

    return sequence_list


if __name__ == '__main__':
    test_file = sys.argv[1]
    temp_seq = read_fasta_seq(open(test_file, 'r'), 'ACGU')
    for e in temp_seq:
        print e