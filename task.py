import numpy as np
import sklearn.metrics
import edlib
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from sklearn.cluster import DBSCAN, OPTICS


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
    return edlib.align(str1, str2)['editDistance']


def parse_centromerefa(p_centromere_file):
    centromere_array = p_centromere_file.readlines()
    del centromere_array[0]
    centromere_str = ""
    for s in centromere_array:
        centromere_str += s[:-1]
    return centromere_str


def parse_monomersfa(p_monomers_file):
    monomers_str_array = []
    monomers_array = p_monomers_file.readlines()
    del monomers_array[0]
    strg = ""
    for s in monomers_array:
        if '>' not in s:
            strg += s[:-1]
        else:
            # strg = str(Seq(strg).reverse_complement())
            monomers_str_array.append(strg)
            strg = ""
    if strg != "":
        # strg = str(Seq(strg).reverse_complement())
        monomers_str_array.append(strg)
    return monomers_str_array


centromere_file = open("cenX.fa")
tsv_file = open("final_decomposition.tsv")
monomers_file = open("monomers.fa")

centromere_str = parse_centromerefa(centromere_file)  # absolutely right
monomers_str_array = parse_monomersfa(monomers_file)  # absolutely right

tsv_array = [s.split() for s in tsv_file.readlines()]  # absolutely right
# blocks_list_char_array = [list(centromere_str[int(s[2]):int(s[3]) + 1]) for s in
#                           tsv_array]  # needs to be List[List[char]]  # absolutely right

blocks_list_char_array = []
for s in tsv_array:
    tmp = centromere_str[int(s[2]):int(s[3]) + 1]
    if "'" in s[1]:
        blocks_list_char_array.append(list(str(Seq(tmp).reverse_complement())))
    else:
        blocks_list_char_array.append(list(tmp))

# blocks_str_array = [str(centromere_str[int(s[2]):int(s[3]) + 1]) for s in tsv_array]  # absolutely right

blocks_str_array = []
for s in tsv_array:
    tmp = str(centromere_str[int(s[2]):int(s[3]) + 1])
    if "'" in s[1]:
        blocks_str_array.append(str(Seq(tmp).reverse_complement()))
    else:
        blocks_str_array.append(tmp)

labels_array = [int(s[1][3:].split("_")[0].replace("'", "")) for s in tsv_array]  # absolutely right
labels_without_repeats_array = []  # absolutely right
for x in labels_array:
    if x not in labels_without_repeats_array:
        labels_without_repeats_array.append(x)

labels_and_monomers_map = dict(
    zip(sorted(labels_without_repeats_array), monomers_str_array))  # map: labels -> centres # absolutely right

monomers_and_labels_map = dict(
    zip(monomers_str_array, sorted(labels_without_repeats_array)))  # map: centres -> labels # absolutely right

centres_and_mblocks_map = {}
for i in range(len(blocks_str_array)):
    tmp = labels_and_monomers_map[labels_array[i]]
    if not centres_and_mblocks_map.keys().__contains__(tmp):
        centres_and_mblocks_map[tmp] = [blocks_str_array[i]]
    else:
        centres_and_mblocks_map[tmp].append(blocks_str_array[i])


def silhouette_score(blocks_char_array, labels, num=500):
    max_size = -1
    for str1 in blocks_char_array:
        max_size = max(max_size, len(str1))

    # reads_lst = []
    # i = 0
    # for x in blocks_char_array:
    #     reads_lst.append(SeqRecord(Seq("".join(x)), id="block_" + str(i), name="block_" + str(i), description=""))
    #     i += 1

    for str1 in blocks_char_array:
        str1 += ['N' for _ in range(max_size - len(str1))]

    blocks_str_np_array = np.array(blocks_char_array)
    blocks_str_np_array_num = blocks_str_np_array.view(np.int32)
    # with open("temp.txt", 'w') as f:
    #     for line in blocks_str_np_array_num:
    #         f.write(" ".join([str(x) for x in list(line)]) + '\n')
    #
    # with open("blocks.fa", "w") as f:
    #     fasta_out = FastaIO.FastaWriter(f, wrap=None)
    #     fasta_out.write_file(reads_lst)

    labels_np_array = np.array(labels)

    return sklearn.metrics.silhouette_score(X=blocks_str_np_array_num[:num], labels=labels_np_array[:num],
                                            metric=edlib_edit_distance)


def dbscan(blocks_char_array):
    max_size = -1
    for str1 in blocks_char_array:
        max_size = max(max_size, len(str1))

    for str1 in blocks_char_array:
        str1 += ['N' for _ in range(max_size - len(str1))]

    blocks_str_np_array = np.array(blocks_char_array)
    blocks_str_np_array_num = blocks_str_np_array.view(np.int32)

    # clustering = OPTICS(max_eps=3, metric=edlib_edit_distance, n_jobs=2, cluster_method='dbscan').fit(
    #     blocks_str_np_array_num)
    clustering = DBSCAN(eps=3, metric=edlib_edit_distance, n_jobs=2).fit(blocks_str_np_array_num)
    return clustering.labels_


# def agglomerative_clustering(blocks_char_array):
#     return


def distortion(blocks, labels,
               labels_and_centres_map):  # blocks array and labels array need to be 'connected' (blocks[i] has labels[i])
    coef = 1 / len(blocks)
    res = 0
    for i in range(len(blocks)):
        res += edlib_edit_distance(blocks[i], labels_and_centres_map[labels[i]]) ** 2
    return coef * res


def DBI(monomers,
        centres_and_blocks_map):  # monomers -- iterable smth of monomers; centres_and_mblocks -- map:str -> list
    coef = 1 / len(monomers)
    res = 0
    for m in monomers:
        maximum = -1
        for m1 in monomers:
            if m != m1:
                dist = edlib_edit_distance(m, m1)
                R1 = R(m, centres_and_blocks_map[m])
                R2 = R(m1, centres_and_blocks_map[m1])
                maximum = max(maximum, (R1 + R2) / dist)
        res += maximum
    return coef * res


def R(m, m_blocks):
    sum_distance = 0
    for mur in m_blocks:
        sum_distance += edlib_edit_distance(m, mur)
    return sum_distance / len(m_blocks)


# print("Silhouette score:", silhouette_score(blocks_list_char_array, labels_array))

# print("DBI:", DBI(monomers_str_array, centres_and_mblocks_map))

# print("Distortion:", distortion(blocks_str_array, labels_array, labels_and_monomers_map))

print("meow")
print("DBSCAN:\n", dbscan(blocks_list_char_array))

# print(silhouette_score(blocks_list_char_array[:500], dbscan(blocks_list_char_array[:500])))
centromere_file.close()
tsv_file.close()
monomers_file.close()
