__author__ = 'aleeee'

import sys

ALPHABET = 'ACGT'


class Seq:
    def __init__(self, name, seq, no):
        self.name = name
        self.seq = seq.upper()
        self.no = no
        self.length = len(seq)

    def __str__(self):
        """Output seq when 'print' method is called."""
        return ">%s\tNo:%s\tlength:%s\n%s" % (self.name, str(self.no), str(self.length), self.seq)


def frequency(tol_str, tar_str):
    """Generate the frequency of tar_str in tol_str.

    :param tol_str: mother string.
    :param tar_str: substring.
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


def write_libsvm(vector_list, label_list, write_file):
    """Write the vector into disk in livSVM format."""
    len_vector_list = len(vector_list)
    len_label_list = len(label_list)
    if len_vector_list == 0:
        sys.stderr.write("The vector is none.")
        sys.exit(1)
    if len_label_list == 0:
        sys.stderr.write("The label is none.")
        sys.exit(1)
    if len_vector_list != len_label_list:
        sys.stderr.write("The length of vector and label is different.")
        sys.exit(1)

    f = open(write_file, 'w')
    len_vector = len(vector_list[0])
    for i in range(0, len_vector_list):
        temp_write = str(label_list[i])
        for j in range(0, len_vector):
            temp_write += ' ' + str(j + 1) + ':' + str(vector_list[i][j])
        f.write(temp_write)
        f.write('\n')
    f.close()


def is_under_alphabet(s, alphabet):
    """Judge the string is within the scope of the alphabet or not.

    :param s: The string.
    :param alphabet: alphabet.

    Return True or the error character.
    """
    for e in s:
        if e not in alphabet:
            return e

    return True


def is_fasta(seq):
    """Judge the Seq object is in FASTA format.
    Two situation:
    1. No seq name.
    2. Seq name is illegal.
    3. No sequence.

    :param seq: Seq object.
    """
    if not seq.name:
        error_info = 'Error, sequence ' + str(seq.no) + ' has no sequence name.'
        sys.stderr.write(error_info)
        return False
    if -1 != seq.name.find('>'):
        error_info = 'Error, sequence ' + str(seq.no) + ' name has > character.'
        sys.stderr.write(error_info)
        return False
    if 0 == seq.length:
        error_info = 'Error, sequence ' + str(seq.no) + ' is null.'
        sys.stderr.write(error_info)
        return False

    return True


def read_fasta(f):
    """Read a fasta file.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return Seq obj list.
    """
    name, seq = '', ''
    count = 0
    seq_list = []
    lines = f.readlines()
    for line in lines:
        if not line:
            break

        if '>' == line[0]:
            if 0 != count or (0 == count and seq != ''):
                if is_fasta(Seq(name, seq, count)):
                    seq_list.append(Seq(name, seq, count))
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    count += 1
    if is_fasta(Seq(name, seq, count)):
        seq_list.append(Seq(name, seq, count))
    else:
        sys.exit(0)

    return seq_list


def read_fasta_yield(f):
    """Yields a Seq object.

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
                if is_fasta(Seq(name, seq, count)):
                    yield Seq(name, seq, count)
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    if is_fasta(Seq(name, seq, count)):
        yield Seq(name, seq, count)
    else:
        sys.exit(0)


def read_fasta_check_dna(f):
    """Read the fasta file, and check its legality.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return the seq list.
    """
    seq_list = []
    for e in read_fasta_yield(f):
        # print e
        res = is_under_alphabet(e.seq, ALPHABET)
        if res:
            seq_list.append(e)
        else:
            error_info = 'Sorry, sequence ' + str(e.no) \
                         + ' has character ' + res + '.(The character must be A or C or G or T)'
            sys.stderr(error_info)
            sys.exit(0)

    return seq_list


def get_sequence_check_dna(f):
    """Read the fasta file.

    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return the sequence list.
    """
    sequence_list = []
    for e in read_fasta_yield(f):
        # print e
        res = is_under_alphabet(e.seq, ALPHABET)
        if res is not True:
            error_info = 'Sorry, sequence ' + str(e.no) \
                         + ' has character ' + str(res) + '.(The character must be A, C, G or T)'
            sys.stderr(error_info)
            sys.exit(0)
        else:
            sequence_list.append(e.seq)

    return sequence_list


def is_sequence_list(sequence_list):
    """Judge the sequence list is within the scope of alphabet."""
    count = 0
    new_sequence_list = []

    for e in sequence_list:
        count += 1
        res = is_under_alphabet(e, ALPHABET)
        if res is not True:
            error_info = 'Sorry, sequence ' + str(count) \
                         + ' has character ' + str(res) + '.(The character must be A, C, G or T)'
            sys.stderr.write(error_info)
            return False
        else:
            new_sequence_list.append(e.upper())

    return True


def get_data(data):
    if isinstance(data, str):
        import os

        if os.path.isfile(data):
            sequence_list = get_sequence_check_dna(open(data, 'r'))
            return sequence_list
        else:
            sys.stderr.write("Sorry, the input file is not exist.")
            sys.exit(0)
    elif isinstance(data, list):
        if is_sequence_list(data):
            return data
        else:
            sys.exit(0)


if __name__ == '__main__':
    # test_file = sys.argv[1]
    # temp_seq = read_fasta_seq(open(test_file, 'r'), 'ACGT')
    # temp_seq = read_fasta(open(test_file, 'r'))

    # temp_seq = read_fasta_check_dna(open(test_file))
    test_file = ['AAAAAAAAAAAAAAAAAAA', 'CCCCCCCCCCCCCCCCCCCCCCCCC']
    temp_seq = get_data(test_file)
    for e in temp_seq:
        print e