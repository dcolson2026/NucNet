import subprocess
import os
import numpy as np

home = os.environ['HOME']
os.chdir(f"{home}")

# DATA=np.loadtxt("mumpower_fission_barriers.dat")
# Z_data = DATA[:, 0]
# A_data = DATA[:, 2]
# FB_MeV_data = DATA[:, 3]
# num_entries = len(Z_data)
# FB_dict = {}
# for i in range(num_entries):
#     FB_dict[(int(Z_data[i]), int(A_data[i]))] = FB_MeV_data[i]

# # Parsing text file
# with open("sf_reactions.txt", "r") as outfile:
#     data1 = outfile.readlines()
# data_list1 = []
# for line in data1:
#     line = line.strip()
#     data_list1.append(line.split("\n"))
# Z_and_A_list = data_list1[3::6]
Z_A_tuples = []
for Z in range(85, 86):
    for A in [272,273,276]:
        if A - Z >= 204:
            continue
        Z_A_tuples.append((Z, A))

other_fw = open(f"{home}/Downloads/ALL_GEF/file.in", "w")
for i in Z_A_tuples:
    Z, A = i
    print(f'"in/{Z}_{A}_SF.in"', file=other_fw)
    fw = open(f"{home}/Downloads/ALL_GEF/in/{Z}_{A}_SF.in", "w")
    print("1", file=fw)
    print("0", file=fw)
    print(f'{Z}, {A}, "GS"', file=fw)
    print("Options(lmd)", file=fw)
    print("END", file=fw)
    fw.close()
other_fw.close()
