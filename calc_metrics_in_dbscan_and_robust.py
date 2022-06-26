from Bio.Seq import Seq
from task import parse_monomersfa, silhouette_score, distortion, DBI, parse_centromerefa


def parse_files(labels_array):
    tsv_file = open("final_decomposition.tsv")
    monomers_file = open("monomers.fa")
    centromere_file = open("cenX.fa")

    centromere_str = parse_centromerefa(centromere_file)  # absolutely right
    monomers_str_array = parse_monomersfa(monomers_file)  # absolutely right

    tsv_array = [s.split() for s in tsv_file.readlines()[:1000]]  # absolutely right
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

    return blocks_list_char_array, blocks_str_array, labels_and_monomers_map, monomers_str_array, centres_and_mblocks_map


def calc_metrics(labels_filename, text):
    labels_array = []
    with open(labels_filename, "r") as labels_file:
        for line in labels_file:
            labels_array.append(int(line))
    blocks_list_char_array, blocks_str_array, labels_and_monomers_map, monomers_str_array, centres_and_mblocks_map = parse_files(
        labels_array)

    print(text, "silhouette score:", silhouette_score(blocks_list_char_array, labels_array, 1000))
    print(text, "distortion:", distortion(blocks_str_array, labels_array, labels_and_monomers_map))
    print(text, "DBI", DBI(centres_and_mblocks_map.keys(), centres_and_mblocks_map))


calc_metrics("out_hdbscan_1000.txt", "dbscan")
calc_metrics("out_robust_1000.txt", "robust")
