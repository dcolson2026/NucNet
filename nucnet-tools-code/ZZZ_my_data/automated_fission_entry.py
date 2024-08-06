"""Automatic fission entry from GEF lmd files"""
import os
import contextlib

from gef_analysis_copy import info_for_network

home = os.environ['HOME']

# Parsing text file
with open(f"{home}/nucnet-tools-code/data_pub/sf_reactions.txt", "r") as outfile:
    data1 = outfile.readlines()

data_list1 = []
for line in data1:
    line = line.strip()
    data_list1.append(line.split("\n"))

fw = open(f"{home}/nucnet-tools-code/ZZZ_my_data/sf_reactions_complete.txt", "w")

for index, element in enumerate(data_list1):
    if (index - 3) % 6 == 0:
        Z = element[0].split("  ")[0]
        A = element[0].split("  ")[1]
        print(Z, A)
    elif (index - 5) % 6 == 0:
        with contextlib.redirect_stdout(fw):
            info_for_network(int(Z), int(A))

        #print(info_for_network(int(Z), int(A)), file=fw)
        #print("poop", file=fw)
    print(*element, file=fw)

fw.close()
