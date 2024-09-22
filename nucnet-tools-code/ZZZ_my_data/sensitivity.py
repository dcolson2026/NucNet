"""Runs all sensitivity calculations"""
# sensititivty paper: https://link.aps.org/accepted/10.1103/PhysRevLett.118.072701
# outline:
# what nuclides are we looking at?
# change XML reaction file to divide by 10, no change, multiply by 10. only change XML file if property exists

# run baseline calculation, get all values
# for ONE nuclide, /10 and run calc
# then x10 and run calc
# from what i understand, the indexing is only over the baseline values


import os
import subprocess
import concurrent.futures
import xml.etree.ElementTree as ET

elements_dict = {
    1: "h",
    2: "he",
    3: "li",
    4: "be",
    5: "b",
    6: "c",
    7: "n",
    8: "o",
    9: "f",
    10: "ne",
    11: "na",
    12: "mg",
    13: "al",
    14: "si",
    15: "p",
    16: "s",
    17: "cl",
    18: "ar",
    19: "k",
    20: "ca",
    21: "sc",
    22: "ti",
    23: "v",
    24: "cr",
    25: "mn",
    26: "fe",
    27: "co",
    28: "ni",
    29: "cu",
    30: "zn",
    31: "ga",
    32: "ge",
    33: "as",
    34: "se",
    35: "br",
    36: "kr",
    37: "rb",
    38: "sr",
    39: "y",
    40: "zr",
    41: "nb",
    42: "mo",
    43: "tc",
    44: "ru",
    45: "rh",
    46: "pd",
    47: "ag",
    48: "cd",
    49: "in",
    50: "sn",
    51: "sb",
    52: "te",
    53: "i",
    54: "xe",
    55: "cs",
    56: "ba",
    57: "la",
    58: "ce",
    59: "pr",
    60: "nd",
    61: "pm",
    62: "sm",
    63: "eu",
    64: "gd",
    65: "tb",
    66: "dy",
    67: "ho",
    68: "er",
    69: "tm",
    70: "yb",
    71: "lu",
    72: "hf",
    73: "ta",
    74: "w",
    75: "re",
    76: "os",
    77: "ir",
    78: "pt",
    79: "au",
    80: "hg",
    81: "tl",
    82: "pb",
    83: "bi",
    84: "po",
    85: "at",
    86: "rn",
    87: "fr",
    88: "ra",
    89: "ac",
    90: "th",
    91: "pa",
    92: "u",
    93: "np",
    94: "pu",
    95: "am",
    96: "cm",
    97: "bk",
    98: "cf",
    99: "es",
    100: "fm",
    101: "md",
    102: "no",
    103: "lr",
    104: "rf",
    105: "db",
    106: "sg",
    107: "bh",
    108: "hs",
    109: "mt",
    110: "ds",
    111: "rg",
    112: "cn",
    113: "nh",
    114: "fl",
    115: "mc",
    116: "lv",
    117: "ts",
    118: "og"
}

# Number of CPUs (default is 8 if not specified in the environment)
NUM_CPUS = int(os.environ.get('MAPP_NR_CPUS', 12))

def run_parallel(tasks):
    """Runs commands in parallel. Accepts tasks, which is a list of tuples of form (function, args).
    Args must be iterables"""
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_CPUS) as executor:
        # Submit all tasks to the executor
        # futures = [executor.submit(func) for _ in range(num_instances)]
        futures = [executor.submit(func, *args) for func, args in tasks]

        # Wait for all futures to complete
        for future in concurrent.futures.as_completed(futures):
            pass

def update_single_rate(xml_file, Z, N, val):
    """Updates the XML file rate"""
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    the_reactant = elements_dict[Z] + str(Z+N)

    # Iterate over each <reaction> element
    for reaction in root.findall('reaction'):
        # Find the <reactant> element and check its value
        reactant = reaction.find('reactant')
        if reactant is not None and reactant.text == the_reactant:
            # Update the <single_rate> element
            single_rate = reaction.find('single_rate')
            if single_rate is not None:
                print(single_rate.text)
                single_rate.text = str(float(single_rate.text) * val)
                print(f"Updated <single_rate> for reactant {the_reactant}")

    # Write the updated XML back to the file
    tree.write(xml_file, encoding='UTF-8', xml_declaration=True)

def run_single_zone(net_file_param, zone_file_param, output_file_param, condition_param):
    """Does ./run_single_zone in Python."""
    # Change into proper directory with ./run_single_zone
    os.chdir(f"{os.environ['HOME']}/nucnet-tools-code/my_examples/network")

    # Build the command as a list
    command = [
        "./run_single_zone",
        net_file_param,
        zone_file_param,
        output_file_param,
        condition_param
    ]

    # Run the command and capture the output
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=False)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"Error running command: {result.stderr.decode()}", file=sys.stderr)
    else:
        print("run_single_zone executed successfully.")

def txt_files(output_xml, output_txt):
    """Makes txt file from output xml. Can later be analyzed."""
    # Command for printing mass num and abundances
    command = [
        f"{os.environ['HOME']}/nucnet-tools-code/examples/analysis/print_zone_abundances", # https://sourceforge.net/p/nucnet-tools/discussion/help/thread/79b2fa10/
        output_xml,
        "[last()]"
    ]

    # Redirect output to the output file
    with open(output_txt, 'w', encoding='ascii') as outfile:
        result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, check=False)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"Error making TXT file: {result.stderr.decode()}")
    else:
        print(f"TXT file made successfully, saved to {output_txt}")

def main(Z_N_tuple, val):
    # user edit these
    # variables
    Z = Z_N_tuple[0]
    N = Z_N_tuple[1]

    home = os.environ['HOME']
    edited_xml_file = f"{home}/PycharmProjects/NucNet_data/data_pub/hypo_SFs_nudat3_no_err_plz_copy.xml"  # Your input XML file

    overall_net_file = f"{home}/nucnet-tools-code/data_pub/net_NO_fission_combo.xml"
    zone_xml = f"{home}/nucnet-tools-code/ZZZ_my_data/sensitivity/zone.xml"
    new_output_xml = f"{home}/nucnet-tools-code/ZZZ_my_data/sensitivity/xmy_output_with_fission_{Z}_{N}_{val}.xml"
    x_path = "[z <= 106]" #90 #changed to 106 now for new data!
    output_txt = f"{home}/nucnet-tools-code/ZZZ_my_data/sensitivity/ya_with_fission_{Z}_{N}_{val}.txt"

    update_single_rate(edited_xml_file, Z, N, val)
    # run_single_zone(overall_net_file, zone_xml, new_output_xml, x_path)
    # txt_files(new_output_xml, output_txt)
    update_single_rate(edited_xml_file, Z, N, 1/val) #reverts changes




if __name__ == "__main__":
    main((100, 160), 10)

    # for parallelization
    Z_N_tuple_list = []
    for i in range(Z_range):
        for j in range(N_range):
            Z_N_tuple_list.append((i, j))

    list_funcs = [
    (main, (Z_N_tuple_list[k], 10)) for k in range(Z_range * N_range)
    ]
    run_parallel(list_funcs)

# mkdir sensitivity
# need output txt to have Z, N (or A) OR just read straight from output xml
