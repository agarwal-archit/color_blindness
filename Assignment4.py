#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math


# In[16]:


# Indexing all four characters 
index = {'A': 0,
         'C': 1,
         'G': 2,
         'T': 3,
         '$': 4}

# Offsets observed to select in first column, this must be added to rand to run the select
offset = {'A': -1,
          'C': 45648951,
          'G': 45648952 + 29813353 - 1,
          'T': 45648952 + 29813353 + 29865831 - 1,
          '$': 151100559}

rev_map = {'A': 'T',
           'C': 'G',
           'G': 'C',
           'T': 'A',
           '$': '$'}


# In[17]:


delta = 1000

data_map = np.loadtxt("../data/chrX_bwt/chrX_map.txt", dtype=int)


# In[18]:


def load_sequence():
    text_file = open("../data/chrX_bwt/chrX_last_col.txt", "r")
    last_column = text_file.read().replace('\n', "")
    print(f"characters in last column: {len(last_column)}")
    print(f"Total count of A,C,G,T {last_column.count('A') + last_column.count('C') + last_column.count('G') + last_column.count('T')}")
    print(f"count of A,C,G,T individually {last_column.count('A')},  {last_column.count('C')}, {last_column.count('G')}, {last_column.count('T')}")
    print(f"character set observed {set(last_column)}")
    print()
    return last_column


# In[19]:


# Get 0-1 array of all characters
def get_bin_array(reference_seq):
    binary_array = np.zeros([len(reference_seq), 5])
    for i, ch in enumerate(reference_seq):
        binary_array[i][index[ch]] = 1
    return binary_array


# In[20]:


# Calculates delta sum matrix
def get_delta_sum(binary_arr):
    size = math.ceil(len(binary_arr) / delta + 1)
    sum_arr = np.zeros([size + 1, 5])
    for i in range(1, size + 1):
        sum_arr[i] = sum_arr[i - 1] + np.sum(binary_arr[(i - 1) * delta: min(len(binary_arr), i * delta)], axis=0)
    return sum_arr


# In[21]:


# Returns first and last rank of a character in the range start - end
def rank(start, end, binary_arr, sum_arr, delta, ch):
    (i, j) = (int(start / delta), int(end / delta))
    if sum_arr[i][index[ch]] == sum_arr[j][index[ch]] and \
            np.sum(binary_arr[start:end + 1], axis=0)[index[ch]] == 0:
        return -1, -1

    first = sum_arr[i][index[ch]] + np.sum(binary_arr[int(i * delta): start], axis=0)[index[ch]] + 1
    last = sum_arr[j][index[ch]] + np.sum(binary_arr[j * delta: end + 1], axis=0)[index[ch]]

    return first, last


# In[22]:


# load reads from file
def load_reads():
    text_file = open("../data/chrX_bwt/reads", "r")
    reads = text_file.readlines()
    reads = [read.replace('\n', "").replace('N', "A") for read in reads]
    return reads


# In[23]:


def get_exons():
    red = [(149249757, 149249868), (149256127, 149256423), (149258412, 149258580), (149260048, 149260213),
           (149261768, 149262007), (149264290, 149264400)]

    green = [(149288166, 149288277), (149295542, 149295710), (149293258, 149293554), (149297178, 149297343),
             (149298898, 149299137), (149301420, 149301530)]
    return red, green


# In[24]:


def get_gene_score(start, end, red, green, data_map, l):
    g = False
    r = False
    for i in range(start, end + 1):
        for gene in red:
            if data_map[i] >= gene[0] and data_map[i] + l <= gene[1]:
                r = True
        for gene in green:
            if data_map[i] >= gene[0] and data_map[i] + l <= gene[1]:
                g = True
    if r and g:
        return 0.5, 0.5
    elif g:
        return 0, 1
    elif r:
        return 1, 0
    else:
        return 0, 0


# In[25]:


def match_to_reference_string(read1, binary_arr, sum_arr, delta):
    start = 0
    end = len(binary_arr) - 1
    mis_matches=2
    for ch in reversed(read1):
        (start_rank, end_rank) = rank(start, end, binary_arr, sum_arr, delta, ch)
        if start_rank == -1:
            if mis_matches > 0:
                mis_matches -= 1
                continue
            return start_rank, end_rank
        start = int(offset[ch] + start_rank)
        end = int(offset[ch] + end_rank)
    return start, end


# In[28]:


def read_reverse(read):
    for i, r in enumerate(read):
        read = read[:i] + rev_map[r] + read[i + 1:]
    return read[::-1]


def match(read, arr, counts):
    (start, end) = match_to_reference_string(read, arr, counts, delta)
    if start == -1:
        (start, end) = match_to_reference_string(read_reverse(read), arr, counts, delta)
    return start, end


# In[ ]:





# In[ ]:





# In[ ]:





# In[31]:


(red, green) = get_exons()
reads = load_reads()
p_red, p_green = 0, 0
arr = get_bin_array(load_sequence())


# In[32]:


counts = get_delta_sum(arr)


# In[ ]:





# In[47]:


for read in reads:
    (start, end) = match(read, arr, counts)
    (r, g) = get_gene_score(start, end, red, green, data_map, len(read))
    p_red += r
    p_green += g
    if(p_green!=0):
        print(f"current red - green percentage: {p_red * 100 / p_green}")
print(f"Final red - green percentage: {p_red * 100 / p_green}")


# In[ ]:





# In[ ]:




