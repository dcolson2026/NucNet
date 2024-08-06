import re
import numpy as np
import os
import sys
import subprocess


def read_lmd(filename):
    """Extracts all fission events with (Z1, Z2, A1, A2) post neutron drip. Returns
    the events as a list of tuples, because a dictionary is made later."""

    # Parsing text file
    with open(filename, "r", encoding='ascii') as outfile: # with open(filename, 'r') as file:
        data1 = outfile.readlines()
    event_list = []

    # look at each line
    for line in data1:
        # Strip of whitepsaces and newlines
        stripped = line.strip()

        # else continue and skip junk
        if stripped == "" or stripped[0] == "*":
            continue

        temp = re.findall(r'-?\d+\.?\d*', stripped)
        event_list.append(tuple(temp[3:5] + temp[7:9])) # make immutable for dict later

    # returns all data
    return len(event_list), event_list

def half_life_theory_eq(Z, A, fission_type, t_half_bdf):
    """Half-lives according to theory papers."""

    # 07/10/24 theory paper from Jin. gives lifetime for spontaneous fission
    if fission_type == "spontaneous":
        N = A - Z
        C0 = -195.09227
        C1 = 3.10156
        C2 = -0.04386
        C3 = 1.4030e-6
        C4 = -0.03199
        # yr_to_s = 3.154e7
        t_half_sf = np.exp(2*np.pi*(C0 + C1*A + C2*Z**2 + C3*Z**4 + C4*(N - Z)**2 - (0.13323*(Z**2/A**(1/3)) - 11.64)))
        lifetime = np.log(2) / t_half_sf

    # from theory paper in table
    elif fission_type == "beta delayed":
        lifetime = np.log(2) / t_half_bdf

    # idk
    elif fission_type == "neutron induced":
        lifetime = None
    return lifetime
    


# def pppfind_channels():
#     """takes high yield nuclides and finds reaction(s) associated with them"""
def info_for_network(Z, A):
    """Makes data for fission.txt. We need: rate type, source, number of reactants, Z A, rate (in time), reactions w/ probability.
    If neutron induced, we need: rate type, source, number of reactants, Z1 A1, Z2 A2, # of xs(T9), xs T9, reactions w/ probability."""
    filename = f"Z{Z}_A{A}_sf_E0MeV.lmd"

    num_events = read_lmd(filename)[0]
    event_list = read_lmd(filename)[1]
    reaction_dict = {}
    for i in event_list:
        if reaction_dict.get(i) is None:
            reaction_dict[i] = 1
        else:
            reaction_dict[i] += 1

    pruned_dict = {} # terrible coding here
    threshold = 0.005
    while threshold < 0.030:
        pruned_dict = {}
        for k, v in reaction_dict.items():
            if v > threshold*num_events:
                pruned_dict[k] = v
        if len(pruned_dict) <= 5:
            break
        threshold += 0.001

    # for k, v in reaction_dict.items():
    #     if v > 0.025*num_events:
    #         pruned_dict[k] = v

    # find total for later renormalization
    prune_total = 0
    for v in pruned_dict.values():
        prune_total += v/num_events

    # if fission_type != "neutron induced":
    #     print("single_rate")
    #     print("GEF", fission_type)
    #     print("1")
    #     print(str(Z) + "  " + str(A))
    #     print(half_life_theory_eq(Z, A, fission_type, t_half_bdf))
    # else:
    #     print("rate_table")
    #     print("GEF neutron induced")
    #     print("2")
    #     print("0  1")
    #     print(Z, A)
    #     print("10")
    #     print()

    # take pruned dict that contains most common reactions.
    # print those reactions, include neutrons based on mass number, renormalize
    for k, v in pruned_dict.items():
        num_neutrons = A - int(k[2]) - int(k[3])
        # if fission_type == "neutron induced":
        #     num_neutrons += 1
        neutron_string = "0 1 " * num_neutrons
        #print(*k, neutron_string[:-1], (v/num_events)/prune_total)
        print(k[0], k[2], k[1], k[3], neutron_string[:-1], (v/num_events)/prune_total)


#os.chdir(f"{os.environ['HOME']}/out")
#info_for_network(Z = 106, A = 267)
