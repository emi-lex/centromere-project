import numpy as np
import edlib
from Bio.Seq import Seq
from sklearn.cluster import DBSCAN
from sys import argv

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


def dbscan(blocks_char_array):
    max_size = -1
    for str1 in blocks_char_array:
        max_size = max(max_size, len(str1))

    for str1 in blocks_char_array:
        str1 += ['N' for _ in range(max_size - len(str1))]

    blocks_str_np_array = np.array(blocks_char_array)
    blocks_str_np_array_num = blocks_str_np_array.view(np.int32)

    clustering = DBSCAN(eps=3, metric=edlib_edit_distance, n_jobs=2).fit(blocks_str_np_array_num)
    return clustering.labels_


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

print("meow")
print("DBSCAN:\n", dbscan(blocks_list_char_array))
