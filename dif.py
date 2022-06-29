from sys import argv
import numpy as np
import matplotlib.pyplot as plt

name, first, second = argv


def differences(filename1, filename2):
    cnt = 0
    mx = -1
    d = {}
    tsv_file = open("final_decomposition.tsv")
    tsv_array = [s.split() for s in tsv_file.readlines()]
    with open(filename1, "r") as file1:
        table1 = np.array([s.split() for s in file1.readlines()])
    with open(filename2, "r") as file2:
        table2 = np.array([s.split() for s in file2.readlines()])
        for i, a1, a2, b1, b2, id1, id2 in zip(range(len(table1)), table1.transpose()[2], table2.transpose()[2],
                                               table1.transpose()[3], table2.transpose()[3], table1.transpose()[4],
                                               table2.transpose()[4]):
            if id1.astype(np.float) < 85 or id2.astype(np.float) < 85:
                continue
            if a1 != a2 or b1 != b2:
                delta = max(abs(a1.astype(np.float)-a2.astype(np.float)), abs(b1.astype(np.float)-b2.astype(np.float)))
                mx = max(mx, delta)
                if d.__contains__(delta):
                    d[delta] += 1
                else:
                    d[delta] = 1
                #y[cnt] += max(abs(a1.astype(np.float)-a2.astype(np.float)), abs(b1.astype(np.float)-b2.astype(np.float)))
                cnt += 1
                print("\n\nline:", i, "\nleft1:", a1, "\nleft2:", a2, "\nright1:", b1, "\nright2:", b2)
    x = np.array(list(d.values())).astype(float)
    y = np.array(list(d.keys())).astype(float)
    plt.title("Line graph")
    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.plot(x, y)
    plt.show()
    print(cnt)
    print(mx)


differences(first, second)
