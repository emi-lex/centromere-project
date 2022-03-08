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

centromere_file = open("centromere.fa")
tsv_file = open("final_decomposition.tsv")

centromere_array = centromere_file.readlines()
del centromere_array[0]
centromere_str = ""
for s in centromere_array:
    centromere_str += s[:-1]

tsv_array = [s.split() for s in tsv_file.readlines()]
str_array = [list(centromere_str[int(s[2]):int(s[3]) + 1]) for s in tsv_array]
labels_array = [int(s[1][3:].split("_")[0].replace("'", "50")) for s in tsv_array]
# str_and_labels_map = dict(zip(str_array, labels_array))

max_size = -1
for str in str_array:
    max_size = max(max_size, len(str))

for str in str_array:
    str += ['N' for _ in range(max_size - len(str))]

str_np_array = np.array(str_array)
str_np_array_num = str_np_array.view(np.int32)

labels_np_array = np.array(labels_array)

num = 500

print(sklearn.metrics.silhouette_score(X=str_np_array_num[:num], labels=labels_np_array[:num], metric=edlib_edit_distance))

centromere_file.close()
tsv_file.close()
