from sys import argv

name, cen_file, num = argv


def parse_centromerefa(p_centromere_file):
    centromere_array = p_centromere_file.readlines()
    del centromere_array[0]
    centromere_str = ""
    for s in centromere_array:
        centromere_str += s[:-1]
    return centromere_str


with open("cutted_cen.fa", 'w') as f:
    cen = open(cen_file)
    cen1 = parse_centromerefa(cen)
    end = 171399 + 57820108
    f.write(">cenX_0:57820108-" + str(end) + '\n')
    f.write(cen1[:int(num)])
