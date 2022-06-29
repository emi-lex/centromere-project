from sys import argv
import numpy as np
import matplotlib.pyplot as plt

name, first, second = argv


def diap(diff):
    if diff < 100:
        return str(diff)
    return str(diff // 100 * 100) + "-" + str(diff // 100 * 100 + 100)


def differences(filename1, filename2):
    cnt = 0
    mx = -1
    m = 0
    d = {}
    diap_d = {}
    with open(filename1, "r") as file1:
        table1 = np.array([s.split() for s in file1.readlines()][:16200])
    with open(filename2, "r") as file2:
        table2 = np.array([s.split() for s in file2.readlines()][:16200])
    for i, a1, a2, b1, b2, id1, id2 in zip(range(len(table1)), table1.transpose()[2], table2.transpose()[2],
                                           table1.transpose()[3], table2.transpose()[3], table1.transpose()[4],
                                           table2.transpose()[4]):
        if id1.astype(np.float) < 80 or id2.astype(np.float) < 80:
            continue
        # if a1 != a2 or b1 != b2:
        m = max(m, b1.astype(np.float) - a1.astype(np.float), b2.astype(np.float) - a2.astype(np.float))
        delta = max(abs(a1.astype(np.float) - a2.astype(np.float)), abs(b1.astype(np.float) - b2.astype(np.float)))
        mx = max(mx, delta)
        if d.__contains__(delta):
            d[delta] += 1
        else:
            d[delta] = 1
        dia = diap(delta)
        if diap_d.__contains__(dia):
            diap_d[dia] += 1
        else:
            diap_d[dia] = 1
        cnt += 1
        # print("\n\nline:", i, "\nleft1:", a1, "\nleft2:", a2, "\nright1:", b1, "\nright2:", b2)
    # print(cnt)
    # print(mx)
    y = np.array(list(diap_d.values())).astype(int)
    x = np.array(list(diap_d.keys())).astype(str)
    for a, b in zip(x, y):
        print(a, str(b / cnt * 100) + "%")
    # print(x, "\n", y)
    # fig, ax = plt.subplots()
    plt.title("Сравнение мономерных блоков после 1 итерации и после последней итерации")
    plt.xlabel("различие блоков")
    plt.ylabel("количество таких блоков")
    plt.bar(x, y)
    # plt.xticks(fontsize=5)
    plt.show()


differences(first, second)
