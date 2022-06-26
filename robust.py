import numpy as np
import edlib
from Bio.Seq import Seq
from joblib import Memory
from sys import argv
import hdbscan
import time


name, final_decomposition, cen = argv


def edlib_edit_distance(str1, str2):
    return edlib.align(str1, str2)['editDistance']


def parse_centromerefa(p_centromere_file):
    centromere_array = p_centromere_file.readlines()
    del centromere_array[0]
    centromere_str = ""
    for s in centromere_array:
        centromere_str += s[:-1]
    return centromere_str


def my_hdbscan(blocks_char_array):
    max_size = -1
    for str1 in blocks_char_array:
        max_size = max(max_size, len(str1))

    for str1 in blocks_char_array:
        str1 += ['N' for _ in range(max_size - len(str1))]

    blocks_str_np_array = np.array(blocks_char_array)
    blocks_str_np_array_num = blocks_str_np_array.view(np.int32)

    start = time.time()
    # clusterer = hdbscan.HDBSCAN()
    clusterer = hdbscan.RobustSingleLinkage(algorithm='best', alpha=1.0, metric=edlib_edit_distance, cut=0.125, k=7)
    clusterer.fit(blocks_str_np_array_num)
    end = time.time()
    with open('out_robust_1000.txt', 'w') as f:
        for s in clusterer.labels_:
            f.write(str(s) + '\n')

    return end - start


tsv_file = open(final_decomposition)
centromere_file = open(cen)
tsv_array = [s.split() for s in tsv_file.readlines()]  # absolutely right
centromere_str = parse_centromerefa(centromere_file)  # absolutely right
blocks_list_char_array = []
for s in tsv_array:
    tmp = centromere_str[int(s[2]):int(s[3]) + 1]
    if "'" in s[1]:
        blocks_list_char_array.append(list(str(Seq(tmp).reverse_complement())))
    else:
        blocks_list_char_array.append(list(tmp))

print("robust")
print("robust:\n", my_hdbscan(blocks_list_char_array[:1000]))
