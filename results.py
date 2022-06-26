from Bio.Seq import Seq
from task import parse_centromerefa, parse_monomersfa, silhouette_score, distortion, DBI


def count_metrics(tsv_filename, centromere_filename="cenX.fa", monomere_filename="monomers.fa"):
    tsv_file = open(tsv_filename)
    monomers_file = open(monomere_filename)
    centromere_file = open(centromere_filename)

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

    print("silhouette score:", silhouette_score(blocks_list_char_array, labels_array))
    print("distortion:", distortion(blocks_str_array, labels_array, labels_and_monomers_map))
    print("DBI", DBI(monomers_str_array, centres_and_mblocks_map))


count_metrics("final_decomposition.tsv")
