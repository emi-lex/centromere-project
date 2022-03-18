import numpy as np
import sklearn.metrics
import edlib


def my_edit_distance(str1, str2):
    if len(str1) > len(str2):
        str1, str2 = str2, str1

    distances = range(len(str1) + 1)
    for i2, c2 in enumerate(str2):
        distances_ = [i2 + 1]
        for i1, c1 in enumerate(str1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]


def edlib_edit_distance(str1, str2):
    return edlib.align(str1.astype('U'), str2.astype('U'))['editDistance']


def parse_centromerefa(centromere_file):
    centromere_array = centromere_file.readlines()
    del centromere_array[0]
    centromere_str = ""
    for s in centromere_array:
        centromere_str += s[:-1]
    return centromere_str


def parse_monomersfa(monomers_file):
    monomers_str_array = []
    monomers_array = monomers_file.readlines()
    del monomers_array[0]
    strg = ""
    for s in monomers_array:
        if '>' not in s:
            strg += s[:-1]
        else:
            monomers_str_array.append(strg)
            strg = ""
    if strg != "":
        monomers_str_array.append(strg)
    return monomers_str_array


centromere_file = open("centromere.fa")
tsv_file = open("final_decomposition.tsv")
monomers_file = open("monomers.fa")

centromere_str = parse_centromerefa(centromere_file)
monomers_str_array = parse_monomersfa(monomers_file)

tsv_array = [s.split() for s in tsv_file.readlines()]
blocks_str_array = [list(centromere_str[int(s[2]):int(s[3]) + 1]) for s in tsv_array]
labels_array = [int(s[1][3:].split("_")[0].replace("'", "")) for s in tsv_array]

# blocks_and_labels_map = dict(zip(str(blocks_str_array), labels_array))

max_size = -1
for str in blocks_str_array:
    max_size = max(max_size, len(str))

for str in blocks_str_array:
    str += ['N' for _ in range(max_size - len(str))]

blocks_str_np_array = np.array(blocks_str_array)
blocks_str_np_array_num = blocks_str_np_array.view(np.int32)

labels_np_array = np.array(labels_array)

num = 500

print("Silhouette score: ",
      sklearn.metrics.silhouette_score(X=blocks_str_np_array_num[:num], labels=labels_np_array[:num],
                                       metric=edlib_edit_distance))

# need dict label -> centre (monomers.fa)
labels_and_centres_map = {}


def distortion(blocks, blocks_and_centres_map):
    coef = 1 / len(blocks)
    res = 0
    for block in blocks:
        res += edlib_edit_distance(block, blocks_and_centres_map[block]) ** 2
    return coef * res


def DBI(monomers,
        centres_and_mblocks_map):  # monomers -- iterable smth of monomers; centres_and_mblocks -- map:str -> list
    coef = 1 / len(monomers)
    res = 0
    for m in monomers:
        maximum = -1
        for m1 in monomers:
            if m != m1:
                dist = edlib_edit_distance(m, m1)
                maximum = max(maximum, (R(m, centres_and_mblocks_map[m]) + R(m1, centres_and_mblocks_map[m1])) / dist)
        res += maximum
    return coef * res


def R(m, m_blocks):
    sum_distance = 0
    for b in m_blocks[m]:
        sum_distance += edlib_edit_distance(m, b)
    return sum_distance / len(m_blocks)


centromere_file.close()
tsv_file.close()
monomers_file.close()
