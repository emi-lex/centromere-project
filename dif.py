from sys import argv
import numpy as np
name, first, second = argv

def differences(filename1, filename2):
    cnt = 0
    tsv_file = open("final_decomposition.tsv")
    tsv_array = [s.split() for s in tsv_file.readlines()]
    with open(filename1, "r") as file1:
        table1 = np.array([s.split() for s in file1.readlines()])
    with open(filename2, "r") as file2:
        table2 = np.array([s.split() for s in file2.readlines()])
    for i, a1, a2, b1, b2 in zip(range(len(table1)), table1.transpose()[2], table2.transpose()[2], table1.transpose()[3], table2.transpose()[3]):
        if a1 != a2 or b1 != b2:
            cnt += 1
            print("\n\nline:", i, "\nleft1:", a1, "\nleft2:", a2, "\nright1:", b1, "\nright2:", b2)
    print(cnt)

differences(first, second)
