"""Automatic fission entry from GEF lmd files"""
import os
import contextlib
import numpy as np
from gef_analysis_copy import info_for_network

def theoretical_sf_T12(Z, A):
    N = A - Z
    a = -43.25203
    b = 0.49192
    c = 3674.3927
    d = -9360.6
    e = 580.75058
    t_half_sf_yr = 10**(a*(Z**2/A) + b*(Z**2/A)**2 + c*((N-Z)/(N+Z)) + d*((N-Z)/(N+Z))**2 + e)
    t_half_sf = t_half_sf_yr * 365.25 * 24 * 60 * 60
    return t_half_sf

Z_A_tuples = []
for Z in range(100, 107):
    for A in range(260, 316):
        Z_A_tuples.append((Z, A))

home = os.environ['HOME']
os.chdir(f"{os.environ['HOME']}/out_v2")

fw = open(f"{home}/nucnet-tools-code/ZZZ_my_data/poopa.txt", "w")

for i in Z_A_tuples:
    Z, A = i
    print(Z, A)
    if not os.path.isfile(f"Z{Z}_A{A}_sf_E0MeV.lmd"):
        print("skip")
        continue
    print("single_rate", file=fw)
    print("hypothetical SFs", file=fw)
    print("1", file=fw)
    print(str(Z) + "  " + str(A), file=fw)
    print(10000, file=fw)
    with contextlib.redirect_stdout(fw):
        info_for_network(Z, A)
        print("", file=fw)
fw.close()
