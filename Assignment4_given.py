# -*- coding: utf-8 -*-
"""DA259Assignment4.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1aBXWwGh7C4vRyhtv8rtJNF_n-CMJ0-U7
"""
import numpy as np
import math

# Index for binary array of all four characters
index = {'A': 0,
         'C': 1,
         'G': 2,
         'T': 3,
         '$': 4}

# Offset to select in first column, this must be added to rand to run the select
offset = {'A': -1,
          'C': 45648951,
          'G': 45648952 + 29813353 - 1,
          'T': 45648952 + 29813353 + 29865831 - 1,
          '$': 151100559}

r_map = {'A': 'T',
         'C': 'G',
         'G': 'C',
         'T': 'A',
         '$': '$'}

delta = 1000

file_map = np.loadtxt("../data/chrX_bwt/chrX_map.txt", dtype=int)


def load_reference_sequence():
    text_file = open("../data/chrX_bwt/chrX_last_col.txt", "r")
    last_column = text_file.read().replace('\n', "")
    print(f"Total characters in last column: {len(last_column)}")
    print(
        f"Total of A,C,G,T {last_column.count('A') + last_column.count('C') + last_column.count('G') + last_column.count('T')}")
    print(
        f"Total of A,C,G,T {last_column.count('A')},  {last_column.count('C')}, {last_column.count('G')}, {last_column.count('T')}")
    print(f"Distinct characters {set(last_column)}")
    print()
    return last_column


# Get o-1 array of all characters
def get_binary_array(reference_seq):
    binary_array = np.zeros([len(reference_seq), 5])
    for i, ch in enumerate(reference_seq):
        binary_array[i][index[ch]] = 1
    return binary_array


# Calculates delta sum matrix
def get_delta_sum_metrics(binary_arr):
    size = math.ceil(len(binary_arr) / delta + 1)
    sum_arr = np.zeros([size + 1, 5])
    for i in range(1, size + 1):
        sum_arr[i] = sum_arr[i - 1] + np.sum(binary_arr[(i - 1) * delta: min(len(binary_arr), i * delta)], axis=0)
    return sum_arr


# Returns first and last rank of a character in the range start - end
def rank(start, end, binary_arr, sum_arr, delta, ch):
    (i, j) = (int(start / delta), int(end / delta))
    if sum_arr[i][index[ch]] == sum_arr[j][index[ch]] and \
            np.sum(binary_arr[start:end + 1], axis=0)[index[ch]] == 0:
        return -1, -1

    first = sum_arr[i][index[ch]] + np.sum(binary_arr[int(i * delta): start], axis=0)[index[ch]] + 1
    last = sum_arr[j][index[ch]] + np.sum(binary_arr[j * delta: end + 1], axis=0)[index[ch]]

    return first, last


# load reads from file
def load_reads():
    text_file = open("../data/chrX_bwt/reads", "r")
    reads = text_file.readlines()
    reads = [read.replace('\n', "").replace('N', "A") for read in reads]
    return reads


def get_exons():
    red = [(149249757, 149249868), (149256127, 149256423), (149258412, 149258580), (149260048, 149260213),
           (149261768, 149262007), (149264290, 149264400)]

    green = [(149288166, 149288277), (149295542, 149295710), (149293258, 149293554), (149297178, 149297343),
             (149298898, 149299137), (149301420, 149301530)]
    return red, green


def get_gene_score(start, end, red, green, file_map, l):
    r = False
    g = False
    for i in range(start, end + 1):
        for gene in red:
            if file_map[i] >= gene[0] and file_map[i] + l <= gene[1]:
                r = True
        for gene in green:
            if file_map[i] >= gene[0] and file_map[i] + l <= gene[1]:
                g = True
    if r and g:
        return 0.5, 0.5
    elif r:
        return 1, 0
    elif g:
        return 0, 1
    else:
        return 0, 0


def match_to_reference_string(read, binary_arr, sum_arr, delta, mis_matches=2):
    start = 0
    end = len(binary_arr) - 1
    for ch in reversed(read):
        (start_rank, end_rank) = rank(start, end, binary_arr, sum_arr, delta, ch)
        if start_rank == -1:
            if mis_matches > 0:
                mis_matches -= 1
                continue
            return start_rank, end_rank
        start = int(offset[ch] + start_rank)
        end = int(offset[ch] + end_rank)
    return start, end


def reversed_read(read):
    for i, r in enumerate(read):
        read = read[:i] + r_map[r] + read[i + 1:]
    return read[::-1]


def match(read, arr, counts):
    (start, end) = match_to_reference_string(read, arr, counts, delta)
    if start == -1:
        (start, end) = match_to_reference_string(reversed_read(read), arr, counts, delta)
    return start, end


def execute():
    (red, green) = get_exons()
    reads = load_reads()
    p_red, p_green = 0, 0
    arr = get_binary_array(load_reference_sequence())
    counts = get_delta_sum_metrics(arr)
    for read in reads:
        (start, end) = match(read, arr, counts)
        (r, g) = get_gene_score(start, end, red, green, file_map, len(read))
        p_red += r
        p_green += g
    return p_red, p_green


def main():
    (p, g) = execute()
    print(f"The red to green percentage: {p * 100 / g}")


if __name__ == "__main__":
    main()