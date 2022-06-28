from Bio.Seq import Seq
from task import silhouette_score, distortion, DBI, parse_centromerefa
import numpy as np
import Bio.Application

def save_fasta(filename, orfs):
    with open(filename, "w") as output_handle:
        from Bio.SeqIO import FastaIO
        fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(orfs)


def run_clustal(mappings, clustal_dir,
                pair_name):  # mappings -- лист мономерных блоков -- все последовательности данного кластера, pair_name -- номер лейбла
    from Bio.Align.Applications import ClustalwCommandline
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment
    import os

    cluster_seqs_path = os.path.join(clustal_dir, pair_name + "_seq.fasta")
    aln_file = os.path.join(clustal_dir, pair_name + "_seq.clu")
    if not os.path.isdir(clustal_dir):
        os.makedirs(clustal_dir)
    if len(mappings) == 1:
        save_fasta(cluster_seqs_path, mappings + mappings)
    else:
        save_fasta(cluster_seqs_path, mappings)

    cmd = ClustalOmegaCommandline(infile=cluster_seqs_path, outfile=aln_file, force=True, threads=10)
    stdout, stderr = cmd()
    align = AlignIO.read(aln_file, "fasta")

    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0, ambiguous='N')
    consensus = str(consensus).replace('-', '')
    return consensus


def calc_metrics(labels_filename, monomers_filename, text=""):
    with open(labels_filename, "r") as labels_file:
        labels_array = list(map(int, labels_file.readlines()))

    tsv_file = open("final_decomposition.tsv")
    centromere_file = open("cenX.fa")
    centromere_str = parse_centromerefa(centromere_file)  # absolutely right
    tsv_array = [s.split() for s in tsv_file.readlines()[:1000]]  # absolutely right
    blocks_list_char_array = []
    for s in tsv_array:
        tmp = centromere_str[int(s[2]):int(s[3]) + 1]
        if "'" in s[1]:
            blocks_list_char_array.append(list(str(Seq(tmp).reverse_complement())))
        else:
            blocks_list_char_array.append(list(tmp))
    blocks_str_array = []
    for s in tsv_array:
        tmp = str(centromere_str[int(s[2]):int(s[3]) + 1])
        if "'" in s[1]:
            blocks_str_array.append(str(Seq(tmp).reverse_complement()))
        else:
            blocks_str_array.append(tmp)

    labels_and_clasters_map = dict(zip(np.unique(labels_array), [[] for _ in range(len(np.unique(labels_array)))]))
    for label, block in zip(labels_array, blocks_str_array):
        labels_and_clasters_map[label].append(block)

    monomers_str_array = []
    labels_and_monomers_map = {}
    centres_and_mblocks_map = {}
    for label in labels_array:
        blocks = labels_and_clasters_map[label]
        from Bio.SeqRecord import SeqRecord
        import os
        monomer = run_clustal(
            [SeqRecord(Seq(blocks[i]), id=str(i), name=str(i), description="") for i in range(len(blocks))],
            os.getcwd(), str(label))
        labels_and_monomers_map[label] = monomer
        monomers_str_array.append(monomer)
        centres_and_mblocks_map[monomer] = labels_and_clasters_map[label]
    print(text, "silhouette score:", silhouette_score(blocks_list_char_array, labels_array, 1000))
    print(text, "distortion:", distortion(blocks_str_array, labels_array, labels_and_monomers_map))
    print(text, "DBI", DBI(centres_and_mblocks_map.keys(), centres_and_mblocks_map))

calc_metrics("out_hdbscan_1000.txt", "out_monomers_hdbscan_1000", "dbscan")
# calc_metrics("out_robust_1000.txt", "out_monomers_robust_1000", "robust")
