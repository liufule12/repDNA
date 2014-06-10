__author__ = 'aleeee'

import subprocess
import shlex


def generate_bracket_seq(test_file, bracket_file):
    """
    This is a system command to generate bracket_seq file according receive_file.
    """
    cmd = "RNAfold --noPS"
    args = shlex.split(cmd)
    subprocess.Popen(args, stdin=open(test_file),
                     stdout=open(bracket_file, 'w')).wait()


def match_2st_has_name(bracket_file, matched_file):
    """
    Get matched sequence from 2st structure sequence and brackets.
    File per line is seq_name, seq_old and seq_matched included MFE.
    Input read file name, write file name.
    Return seq, seq is a list of triple(seq_name, seq_old, seq_bracket, seq_matched).
    """
    f_open = open(bracket_file)
    f_write = open(matched_file, 'w')
    lines = f_open.readlines()
    len_lines = len(lines)
    seq = []

    for i in range(0, len_lines, 3):
        # Get seq_name.
        seq_name = ''
        for j in range(1, len(lines[i])):
            if ' ' != lines[i][j]:
                seq_name += lines[i][j]

        # Get the sequence_match sequence.
        sequence_old = lines[i + 1].strip()
        sequence_bracket = lines[i + 2].strip()
        stack_sequence = []
        stack_bracket = []
        dict_match = {}
        len_line = len(sequence_old)
        for j in range(0, len_line):

            if '.' == sequence_bracket[j]:
                dict_match[j] = '.'
            elif '(' == sequence_bracket[j]:
                stack_bracket.append('(')
                stack_sequence.append(j)
                dict_match[j] = '('
            elif ')' == sequence_bracket[j]:
                stack_bracket.pop()
                temp_order = stack_sequence.pop()
                dict_match[temp_order] = sequence_old[j]
                dict_match[j] = sequence_old[temp_order]
            else:
                print j
                print "match_2st error!!!!!!!!!!"
                return

        # Write seq_name, seq_old, seq_matched and include MFE.
        f_write.write(seq_name)
        f_write.write(lines[i + 1])
        str_match_values = ''.join(dict_match.values())
        f_write.write(str_match_values)
        f_write.write(lines[i + 2][len_line:])
        seq.append((seq_name.strip(), sequence_old, lines[i + 2][:len_line], str_match_values))
    return seq


def frequency_pair(sequence_old, sequence_match, pair):
    # Get the frequency of pair.
    frequency = 0
    len_sequence = len(sequence_old)
    for i in range(0, len_sequence):
        if (sequence_old[i] == pair[0] and sequence_match[i] == pair[1]) or \
           (sequence_match[i] == pair[0] and sequence_old[i] == pair[1]):
            frequency += 1
    return frequency


def free_energy(pair):
    # Get the free energy in a pair.
    if '.' == pair[0] or '.' == pair[1]:
        return 0
    elif ('G' == pair[0] and 'C' == pair[1]) or ('C' == pair[0] and 'G' == pair[1]):
        return -3
    elif ('A' == pair[0] and 'U' == pair[1]) or ('U' == pair[0] and 'A' == pair[1]):
        return -2
    elif ('G' == pair[0] and 'U' == pair[1]) or ('U' == pair[0] and 'G' == pair[1]):
        return -1
    else:
        print pair
        print "free energy error!!!!!"
        return -9999999


def seta_rirj(pair_ri, pair_rj):
    # Get seta_RiRj result.
    return (free_energy(pair_ri) - free_energy(pair_rj)) * (free_energy(pair_ri) - free_energy(pair_rj))


def seta(sequence_old, sequence_match, lamada):
    # Get seta_lamada.
    l = len(sequence_old)
    r_sum = 0.0
    for i in range(0, l - lamada):
        ri = [sequence_old[i], sequence_match[i]]
        rj = [sequence_old[i + lamada], sequence_match[i + lamada]]
        r_sum += seta_rirj(ri, rj)

    return r_sum / (l - lamada)


def frequency_kmer(sequence_old, sequence_match, k):
    # Get kmer's frequency.
    sequence_frequency = [0] * (10 ** k)
    len_sequence = len(sequence_old)
    for i in range(0, len_sequence - k + 1):
        temp_index = 0
        j = 0
        for j in range(i, i + k):
            if ('.' == sequence_match[j] and 'A' == sequence_old[j]) or \
                    ('A' == sequence_match[j] and '.' == sequence_old[j]):
                temp_index += 0
            elif ('.' == sequence_match[j] and 'U' == sequence_old[j]) or \
                    ('U' == sequence_match[j] and '.' == sequence_old[j]):
                temp_index += 1 * (10 ** (k - j + i - 1))
            elif ('.' == sequence_match[j] and 'G' == sequence_old[j]) or \
                    ('G' == sequence_match[j] and '.' == sequence_old[j]):
                temp_index += 2 * (10 ** (k - j + i - 1))
            elif ('.' == sequence_match[j] and 'C' == sequence_old[j]) or \
                    ('C' == sequence_match[j] and '.' == sequence_old[j]):
                temp_index += 3 * (10 ** (k - j + i - 1))
            elif 'A' == sequence_match[j] and 'U' == sequence_old[j]:
                temp_index += 4 * (10 ** (k - j + i - 1))
            elif 'U' == sequence_match[j] and 'A' == sequence_old[j]:
                temp_index += 5 * (10 ** (k - j + i - 1))
            elif 'G' == sequence_match[j] and 'C' == sequence_old[j]:
                temp_index += 6 * (10 ** (k - j + i - 1))
            elif 'C' == sequence_match[j] and 'G' == sequence_old[j]:
                temp_index += 7 * (10 ** (k - j + i - 1))
            elif 'G' == sequence_match[j] and 'U' == sequence_old[j]:
                temp_index += 8 * (10 ** (k - j + i - 1))
            elif 'U' == sequence_match[j] and 'G' == sequence_old[j]:
                temp_index += 9 * (10 ** (k - j + i - 1))
            else:
                print 'frequencey_kmer error!!!!!'
                temp_index = -999999
                break
        if j == i + k - 1:
            sequence_frequency[temp_index] += 1

    return sequence_frequency


def calculate_d(d, sequence_old, sequence_match, k, lamada, w):
    """
    Calculate D.
    """
    theta = [0] * lamada

    # Calculate frequency.
    f = frequency_kmer(sequence_old, sequence_match, k)
    temp_f_sum = 0
    for temp_f in f:
        temp_f_sum += temp_f

    # Normalization.
    for j in range(0, len(f)):
        f[j] /= (temp_f_sum * 1.0)

    # Calculate theta.
    sum_theta = 0
    for j in range(1, lamada + 1):
        theta[j - 1] = seta(sequence_old, sequence_match, j)
        sum_theta += theta[j - 1]

    denominator = 1.0 + w * sum_theta

    # Calculate dk.
    for j in range(0, (10 ** k) + lamada):
        if 0 <= j < 10 ** k:
            d[j] = f[j] / denominator
        elif j >= 10 ** k and k < (10 ** k) + lamada:
            d[j] = w * theta[j - 10 ** k] / denominator
        else:
            print 'error!!!!'
    return d


def read_p_value(read_name):
    """
    Read p-value file, return all seq p-vaule.
    """
    f_open = open(read_name)
    p_values = []
    lines = f_open.readlines()
    for e in lines:
        e = e.split('\t')
        p_values.append(float(e[-1]))
    f_open.close()
    return p_values


def read_ngram(read_name):
    """
    Read 64-Ngram, return 64-vector.
    """
    f_open = open(read_name)
    ngram = []
    lines = f_open.readlines()
    for e in lines:
        # print e
        e = e.strip().split(' ')
        ngram.append(e)
    f_open.close()
    return ngram


def generate_p_values(test_file, p_values_file):
    cmd = "randfold -d " + test_file + " 1000"
    args = shlex.split(cmd)
    subprocess.Popen(args, stdout=open(p_values_file, 'w')).wait()
    p_value = read_p_value(p_values_file)
    print 'step 2.1 is completed, p_value ok!'
    return p_value


def generate_n_gram(test_file, n_gram_perl, n_gram):
    cmd = "perl " + n_gram_perl + " " + test_file + " " + n_gram
    args = shlex.split(cmd)
    subprocess.Popen(args).wait()
    n_gram = read_ngram(n_gram)
    print 'step 2.2 is completed, Ngram ok!'
    return n_gram


def generate_mfe(bracket_file):
    f = open(bracket_file)
    lines = f.readlines()
    len_lines = len(lines)
    mfe = []
    for i in range(0, len_lines, 3):
        temp_mfe = ''
        len_line = len(lines[i + 1])
        for j in range(len_line + 2, len(lines[i + 2])):
            if lines[i + 2][j] != ')':
                temp_mfe += lines[i + 2][j]
            else:
                break
        mfe.append(float(temp_mfe))
    print "step 2.3 is completed, mfe ok!"
    return mfe


def make_test_vector_libsvm_format(matched_file, vector_file, k, lamada, w, way_choice, p_value, ngram, mfe, label):
    """
    Generate vector according way_choice.
    Return MFE_list.
    """
    f_open = open(matched_file)
    f_write = open(vector_file, 'w')

    lines = f_open.readlines()
    len_lines = len(lines)
    # vector = {}
    for i in range(0, len_lines, 3):
        sequence_old = lines[i + 1].strip()
        sequence_match = lines[i + 2].strip()

        if way_choice == 'SSC':
            d = [0] * (10 ** k + lamada)
        elif way_choice == 'ExSSC':
            d = [0] * (10 ** k + lamada + 2 + 64)
        else:
            print 'way_choice -------------------------------error!'

        # Calculate D.
        d = calculate_d(d, sequence_old, sequence_match, k, lamada, w)

        # MFE + p_value + Ngram.
        if way_choice == 'ExSSC':
            d[(10 ** k) + lamada] = mfe[i / 3]
            d[(10 ** k) + lamada + 1] = float(p_value[i / 3])
            for j in range(0, len(ngram[i / 3])):
                d[(10 ** k) + lamada + 2 + j] = float(ngram[i / 3][j])

        # Write the livSVM format.
        f_write.write(label)
        len_d = len(d)
        for j in range(0, len_d):
            temp_write = ' ' + str(j + 1) + ':' + str(d[j])
            f_write.write(temp_write)
        f_write.write('\n')

    print 'End of Vector...'