"""Gradient descent approach to r-process calculations"""

import numpy as np
# import sys

def mass_num_vs_abundance_lists(filename):
    """Extracts data from last timestep in r-process calculation.
    Returns mass number list, abundance list, and error of r-process peaks.
    In addition, it now makes sure that the integrals of the r-process regions are equivalent.
    Therefore, the text files are processed and unchanged, but the lists here 
    are scaled accordingly."""

    # All Solar data stored in tuple. Gives mass number, abundance, y error down, y error up
    # Data taken from:
    # https://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1999A%26A...342..881G&defaultprint=YES&filetype=.pdf
    actual_mass_abun_down_up = ([69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, 170.0, 171.0, 172.0, 173.0, 174.0, 175.0, 176.0, 177.0, 178.0, 179.0, 180.0, 181.0, 182.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 189.0, 190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 199.0, 200.0, 201.0, 202.0, 203.0, 204.0, 205.0, 206.0, 207.0, 208.0, 209.0], [6.18, 7.74, 1.96, 0.0, 6.31, 19.7, 3.78, 8.78, 3.76, 0.0, 4.81, 28.1, 4.07, 6.2, 4.38, 23.6, 2.87, 0.0, 0.292, 4.09, 1.11, 2.61, 0.21, 0.062, 0.0987, 0.0, 0.14, 0.0, 0.0808, 0.0739, 0.173, 0.226, 0.267, 0.315, 0.297, 0.337, 0.266, 0.171, 0.211, 0.119, 0.172, 0.156, 0.152, 0.176, 0.124, 0.172, 0.111, 0.0955, 0.15, 0.244, 0.184, 0.214, 0.0836, 0.152, 0.113, 0.22, 0.256, 0.492, 0.848, 1.47, 1.08, 1.58, 0.822, 0.653, 0.309, 0.385, 0.248, 0.33, 0.17, 0.214, 0.157, 0.161, 0.11, 0.066, 0.0706, 0.0998, 0.054, 0.0533, 0.0334, 0.0421, 0.0323, 0.049, 0.0452, 0.0571, 0.0495, 0.0595, 0.0468, 0.0579, 0.0471, 0.0614, 0.0601, 0.0741, 0.0741, 0.09, 0.0972, 0.103, 0.0839, 0.0753, 0.0546, 0.0506, 0.034, 0.0369, 0.0297, 0.0381, 0.0316, 0.0391, 0.0305, 0.0292, 0.0238, 0.0192, 0.0138, 0.0145, 0.0106, 0.0136, 0.0065, 0.0106, 0.0151, 0.0245, 0.0318, 0.0708, 0.103, 0.152, 0.229, 0.273, 0.388, 0.421, 0.445, 0.302, 0.191, 0.095, 0.0507, 0.0334, 0.0265, 0.0257, 0.0033, 0.0266, 0.0497, 0.197, 0.142, 0.0003, 0.0501], [0.0, 6.8, 0.0, 0.0, 0.0, 9.94, 3.24, 7.84, 3.48, 0.0, 0.918, 24.8, 3.04, 5.83, 3.05, 14.2, 1.05, 0.0, 0.0, 0.0, 0.0, 1.26, 0.0, 0.0, 0.0, 0.0, 0.0976, 0.0, 0.0496, 0.0, 0.146, 0.21, 0.23, 0.244, 0.209, 0.298, 0.224, 0.113, 0.178, 0.066, 0.131, 0.136, 0.127, 0.0921, 0.0916, 0.0515, 0.0816, 0.0697, 0.103, 0.151, 0.115, 0.0634, 0.0578, 0.0, 0.0925, 0.195, 0.217, 0.381, 0.663, 1.29, 0.851, 1.42, 0.627, 0.389, 0.283, 0.231, 0.0, 0.26, 0.0, 0.0, 0.0183, 0.0, 0.0545, 0.0, 0.0526, 0.0582, 0.0456, 0.0145, 0.0156, 0.0221, 0.0278, 0.0459, 0.0267, 0.0498, 0.046, 0.0505, 0.0364, 0.0501, 0.0429, 0.0497, 0.0517, 0.0655, 0.0684, 0.0795, 0.089, 0.0827, 0.0728, 0.0691, 0.0495, 0.042, 0.025, 0.0283, 0.0107, 0.0323, 0.0266, 0.0229, 0.0156, 0.0177, 0.0186, 0.01, 0.0109, 0.0, 0.0042, 0.0, 0.0, 0.0, 0.011, 0.0073, 0.027, 0.0633, 0.0961, 0.137, 0.221, 0.252, 0.374, 0.362, 0.394, 0.256, 0.179, 0.0805, 0.0357, 0.0061, 0.0111, 0.0, 0.0, 0.0171, 0.0, 0.0364, 0.0, 0.0, 0.01], [9.37, 8.55, 9.61, 9.93, 8.19, 28.9, 4.68, 9.68, 4.65, 10.3, 5.71, 32.2, 4.87, 6.51, 5.68, 34.5, 4.01, 0.587, 1.01, 4.75, 1.81, 3.01, 0.484, 0.437, 0.27, 0.0602, 0.226, 0.025, 0.112, 0.153, 0.2, 0.25, 0.305, 0.426, 0.375, 0.383, 0.303, 0.228, 0.244, 0.193, 0.207, 0.174, 0.179, 0.25, 0.155, 0.291, 0.136, 0.127, 0.193, 0.375, 0.247, 0.412, 0.113, 0.18, 0.131, 0.242, 0.295, 0.601, 1.03, 1.62, 1.31, 1.74, 1.01, 0.926, 0.341, 0.477, 0.272, 0.396, 0.296, 1.0, 0.248, 0.357, 0.136, 0.131, 0.0811, 0.124, 0.0611, 0.0711, 0.0347, 0.0522, 0.0328, 0.0515, 0.0482, 0.0622, 0.0526, 0.0609, 0.05, 0.0634, 0.0508, 0.0694, 0.0672, 0.0787, 0.0745, 0.0917, 0.098, 0.104, 0.0941, 0.0833, 0.0586, 0.057, 0.0391, 0.0407, 0.0326, 0.0432, 0.0353, 0.0515, 0.0374, 0.0334, 0.0263, 0.0236, 0.016, 0.0214, 0.0144, 0.0215, 0.01, 0.0179, 0.0176, 0.0337, 0.0359, 0.0781, 0.109, 0.168, 0.237, 0.289, 0.402, 0.47, 0.493, 0.347, 0.204, 0.105, 0.0682, 0.064, 0.0426, 0.0677, 0.0271, 0.033, 0.115, 0.379, 0.433, 1.78, 0.164])

    # Parsing text file
    with open(filename + ".txt", "r", encoding='ascii') \
        as outfile: # with open(filename, 'r') as file:
        data1 = outfile.readlines()
    data_list1 = []

    for line in data1:
        if line[0] in ["#", "\n"]:
            continue
        line = line.strip()
        data_list1.append(line.split(" "))
    cut = data_list1[2:] # cuts out top two lines with parameter info
    mass_number_list = [int(i[0]) for i in cut]
    abundance_list = [float(i[2]) for i in cut]

    sum_simulation_peaks = sum(abundance_list[69:210])
    sum_actual_peaks = sum(actual_mass_abun_down_up[1])
    scale_factor = sum_actual_peaks/sum_simulation_peaks

    # make integrals of r-process region same for comparison
    corrected_abundance_list = [i*scale_factor for i in abundance_list]

    error = 0
    for i in range(69, 210): # range specific for the actual data of solar abundances
        if actual_mass_abun_down_up[1][i-69] == 0:
            continue
        error += (actual_mass_abun_down_up[1][i-69] - corrected_abundance_list[i])**2

    return mass_number_list, abundance_list, error, corrected_abundance_list

# fix this function. if it is a list of indices, there is no difference between
# index of the element (because the element is the index)
def find_min_error(filename, file_index_list):
    """CHANGE range later"""
    min_error = np.inf
    for list_index, element in enumerate(file_index_list):
        current_error = mass_num_vs_abundance_lists(filename + str(element))[2]
        print("File index and its error: " + str(list_index) + " and " + str(current_error))
        if  current_error < min_error:
            min_error = current_error
            min_index = list_index
    return min_index, min_error


if __name__ == "__main__":
    print(find_min_error("ya_with_fission", [str(i) for i in range(1, 21)]))
    #print(find_min_error("is_it_still_best", [str(""),]))
