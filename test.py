import numpy as np

filename1, filename2 = "final_decomposition.tsv", "final_decomposition_iter_0.tsv"
with open(filename1, "r") as file1:
    table1 = np.array([s.split() for s in file1.readlines()])
with open(filename2, "r") as file2:
    table2 = np.array([s.split() for s in file2.readlines()])

table1 = np.concatenate((table1[:16220], table1[16230:]))

# for i, a1, a2, id1, id2 in zip(range(len(table1)), table1.transpose()[2], table2.transpose()[2], table1.transpose()[4],
#                                table2.transpose()[4]):
#
#     if abs(int(a1) - int(a2)) > 4 and np.float(id1) > 80 and np.float(id2) > 80:
#         print(i)
#         break
#
mask = (table1.transpose()[4].astype(np.float) > 80) & (table2.transpose()[4].astype(np.float) > 80)
# # print(table2.transpose()[2][16150:16250].astype(int) - table1.transpose()[2][16150:16250].astype(int))
# diffs = zip(np.arange(len(table1))[mask],
#             table2.transpose()[2][mask].astype(int) - table1.transpose()[2][mask].astype(int))
# print(np.arange(len(table1))[mask])
# for i, diff in diffs:
#     if 16100 < i < 16300:
#         print(i, diff)
diffs = table2.transpose()[2].astype(int) - table1.transpose()[2].astype(int)
for i, diff in zip(np.arange(len(table2)), diffs):
    if diff != 0:
        print(i, diff)
# print(table2.transpose()[2][1600:1700].astype(int) - table1.transpose()[2][1600:1700].astype(int))
