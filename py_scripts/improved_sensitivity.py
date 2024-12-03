"""Runs all sensitivity calculations"""
# sensititivty paper: https://link.aps.org/accepted/10.1103/PhysRevLett.118.072701

import os
import fcntl #locks file until unlocked
import subprocess
import concurrent.futures
import sys
import xml.etree.ElementTree as ET
import time


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
#entire_Z_N_tuple_list = [(90, 142), (92, 142), (92, 144), (92, 146), (94, 142), (94, 144), (94, 146), (94, 148), (94, 150), (96, 144), (96, 146), (96, 148), (96, 150), (96, 152), (96, 154), (98, 140), (98, 144), (98, 148), (98, 150), (98, 152), (98, 154), (98, 156), (98, 158), (100, 141), (100, 142), (100, 144), (100, 146), (100, 150), (100, 152), (100, 154), (100, 156), (100, 158), (100, 159), (100, 160), (100, 161), (100, 162), (100, 163), (100, 164), (100, 165), (100, 166), (100, 167), (100, 168), (100, 169), (100, 170), (100, 171), (100, 172), (100, 173), (100, 174), (100, 175), (100, 176), (100, 177), (100, 178), (100, 179), (100, 180), (100, 181), (100, 182), (100, 183), (100, 184), (100, 185), (100, 186), (100, 187), (100, 188), (100, 189), (100, 190), (100, 191), (100, 192), (100, 193), (100, 194), (100, 195), (100, 196), (100, 197), (100, 198), (100, 199), (100, 200), (100, 201), (100, 202), (100, 203), (101, 158), (101, 159), (101, 160), (101, 161), (101, 162), (101, 163), (101, 164), (101, 165), (101, 166), (101, 167), (101, 168), (101, 169), (101, 170), (101, 171), (101, 172), (101, 173), (101, 174), (101, 175), (101, 176), (101, 177), (101, 178), (101, 179), (101, 180), (101, 181), (101, 182), (101, 183), (101, 184), (101, 185), (101, 186), (101, 187), (101, 188), (101, 189), (101, 190), (101, 191), (101, 192), (101, 193), (101, 194), (101, 195), (101, 196), (101, 197), (101, 198), (101, 199), (101, 200), (101, 201), (101, 202), (101, 203), (102, 148), (102, 150), (102, 152), (102, 156), (102, 158), (102, 159), (102, 160), (102, 161), (102, 162), (102, 163), (102, 164), (102, 165), (102, 166), (102, 167), (102, 168), (102, 169), (102, 170), (102, 171), (102, 172), (102, 173), (102, 174), (102, 175), (102, 176), (102, 177), (102, 178), (102, 179), (102, 180), (102, 181), (102, 182), (102, 183), (102, 184), (102, 185), (102, 186), (102, 187), (102, 188), (102, 189), (102, 190), (102, 191), (102, 192), (102, 193), (102, 194), (102, 195), (102, 196), (102, 197), (102, 198), (102, 199), (102, 200), (102, 201), (102, 202), (102, 203), (103, 157), (103, 158), (103, 159), (103, 160), (103, 161), (103, 162), (103, 163), (103, 164), (103, 165), (103, 166), (103, 167), (103, 168), (103, 169), (103, 170), (103, 171), (103, 172), (103, 173), (103, 174), (103, 175), (103, 176), (103, 177), (103, 178), (103, 179), (103, 180), (103, 181), (103, 182), (103, 183), (103, 184), (103, 185), (103, 186), (103, 187), (103, 188), (103, 189), (103, 190), (103, 191), (103, 192), (103, 193), (103, 194), (103, 195), (103, 196), (103, 197), (103, 198), (103, 199), (103, 200), (103, 201), (103, 202), (103, 203), (104, 149), (104, 150), (104, 152), (104, 154), (104, 156), (104, 157), (104, 158), (104, 159), (104, 160), (104, 161), (104, 162), (104, 163), (104, 164), (104, 165), (104, 166), (104, 167), (104, 168), (104, 169), (104, 170), (104, 171), (104, 172), (104, 173), (104, 174), (104, 175), (104, 176), (104, 177), (104, 178), (104, 179), (104, 180), (104, 181), (104, 182), (104, 183), (104, 184), (104, 185), (104, 186), (104, 187), (104, 188), (104, 189), (104, 190), (104, 191), (104, 192), (104, 193), (104, 194), (104, 195), (104, 196), (104, 197), (104, 198), (104, 199), (104, 200), (104, 201), (104, 202), (104, 203), (105, 155), (105, 156), (105, 157), (105, 158), (105, 159), (105, 160), (105, 161), (105, 162), (105, 163), (105, 164), (105, 165), (105, 166), (105, 167), (105, 168), (105, 169), (105, 170), (105, 171), (105, 172), (105, 173), (105, 174), (105, 175), (105, 176), (105, 177), (105, 178), (105, 179), (105, 180), (105, 181), (105, 182), (105, 183), (105, 184), (105, 185), (105, 186), (105, 187), (105, 188), (105, 189), (105, 190), (105, 191), (105, 192), (105, 193), (105, 194), (105, 195), (105, 196), (105, 197), (105, 198), (105, 199), (105, 200), (105, 201), (105, 202), (105, 203), (106, 152), (106, 154), (106, 155), (106, 156), (106, 157), (106, 158), (106, 159), (106, 160), (106, 161), (106, 162), (106, 163), (106, 164), (106, 165), (106, 166), (106, 167), (106, 168), (106, 169), (106, 170), (106, 171), (106, 172), (106, 173), (106, 174), (106, 175), (106, 176), (106, 177), (106, 178), (106, 179), (106, 180), (106, 181), (106, 182), (106, 183), (106, 184), (106, 185), (106, 186), (106, 187), (106, 188), (106, 189), (106, 190), (106, 191), (106, 192), (106, 193), (106, 194), (106, 195), (106, 196), (106, 197), (106, 198), (106, 199), (106, 200), (106, 201), (106, 202), (106, 203)]


# Number of CPUs (default is 8 if not specified in the environment)
NUM_CPUS = int(os.environ.get('MAPP_NR_CPUS', 16))

def run_parallel(tasks):
    """Runs commands in parallel. Accepts tasks, which is a list of tuples of form (function, args).
    Args must be iterables"""
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_CPUS) as executor:
        # # Submit all tasks to the executor
        # # futures = [executor.submit(func) for _ in range(num_instances)]
        # futures = [executor.submit(func, *args) for func, args in tasks]

        # # Wait for all futures to complete
        # for future in concurrent.futures.as_completed(futures):
        #     pass
        for x in range(0, len(tasks), 16):
            # Submit a batch of tasks
            batch = tasks[x:x + 16]
            futures = [executor.submit(func, *args) for func, args in batch]

            # Wait for the current batch to complete
            for future in concurrent.futures.as_completed(futures):
                pass

def file_lock(filename):
    """Locks file and can only be used by process"""
    file = open(filename, 'r+', encoding="UTF-8")
    fcntl.flock(file, fcntl.LOCK_EX)
    #print(f"{filename} is now locked.")

def file_unlock(filename):
    """Unlocks file"""
    file = open(filename, 'r+', encoding="UTF-8")
    fcntl.flock(file, fcntl.LOCK_UN)
    file.close()
    #print(f"{filename} is now UNlocked.")

def update_single_rate(xml_file, Z, N, val):
    """Updates the XML file rate"""
    # Parse the XML file
    found_bool = False
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
                single_rate.text = str(float(single_rate.text) * val)
                print(f"Updated <single_rate> for reactant {the_reactant} to {single_rate.text} in {xml_file}")
                tree.write(xml_file, encoding='UTF-8', xml_declaration=True)
                found_bool = True
                time.sleep(3)
    if found_bool:
        return True

    #print(f"reactant {the_reactant} was not found")
    return False

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

    print(f"Single zone is running for output {output_file_param}")

    # Run the command and capture the output
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=False)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"Error running command: {result.stderr.decode()} for {output_file_param}", \
              file=sys.stderr)
    else:
        print(f"run_single_zone executed successfully for {output_file_param}.")

def txt_files(output_xml, output_txt):
    """Makes txt file from output xml. Can later be analyzed."""
    # Command for printing mass num and abundances
    # https://sourceforge.net/p/nucnet-tools/discussion/help/thread/79b2fa10/
    # The output is Z, A, Y, X,
    command = [
        f"{os.environ['HOME']}/nucnet-tools-code/examples/analysis/print_zone_abundances",
        output_xml,
        "[last()]"
    ]

    # Redirect output to the output file
    with open(output_txt, 'w', encoding='ascii') as outfile:
        result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, check=False)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"Error making TXT file: {result.stderr.decode()} for {output_xml}")
    else:
        print(f"TXT file made successfully, saved to {output_txt}")

def main(Z_N_tuple, inc_or_dec, reac_xml_num, output_folder):
    """main function that does it all"""
    # user edit these
    # variables
    Z = Z_N_tuple[0]
    N = Z_N_tuple[1]
    if inc_or_dec == "inc":
        #val = 10
        val = 2
    else:
        #val = 0.1
        val = 1/2

    home = os.environ['HOME']
    edited_xml_file = f"{home}/nucnet-tools-code/cleaner_data_pub/A_inc_bdf_{reac_xml_num}.xml" # bbb for dec
    overall_net_file = f"{home}/nucnet-tools-code/cleaner_data_pub/A_separated_reactions_inc_{reac_xml_num}.xml"
    zone_xml = f"{home}/nucnet-tools-code/ZZZ_my_data/sensitivity/sens_zone.xml"
    new_output_xml = f"{home}/nucnet-tools-code/ZZZ_my_data/sensitivity/{output_folder}/xmy_output_with_fission_{Z}_{N}_{inc_or_dec}.xml"
    x_path = "[z <= 106]" #90 #changed to 106 now for new data!
    output_txt = f"{home}/nucnet-tools-code/ZZZ_my_data/sensitivity/{output_folder}/ya_with_fission_{Z}_{N}_{inc_or_dec}.txt"

    # Function calls
    file_lock(edited_xml_file) # LOCKS file
    reactant_found = update_single_rate(edited_xml_file, Z, N, val)
    if not reactant_found:
        return
    if reactant_found:
        print("OFFICIALLY LOCKED AND FOUND")
    time.sleep(5)
    run_single_zone(overall_net_file, zone_xml, new_output_xml, x_path)
    txt_files(new_output_xml, output_txt)
    update_single_rate(edited_xml_file, Z, N, 1/val) #reverts changes
    file_unlock(edited_xml_file) # UNLOCKS file
    time.sleep(5)



if __name__ == "__main__":
    #main((100, 160), 10)
    #Z_range = list(range(52,84))
    #N_range = list(range(68, 126))

    # # for parallelization
    # Z_N_tuple_list = []
    # for i in range(90, 107): # 90, 107
    #     for j in range(130, 220): # 146, 204 u238 to sg309
    #         Z_N_tuple_list.append((i, j))

    # why have i ever been using this stupid Z_N_tuple list??? JUST GET ALL OF THEM BEFOREHAND!!!! BY PARSING THE FILE omg...brb omg
    bdf_nuclides_list = [(97, 163), (97, 167), (97, 187), (97, 188), (97, 189), (97, 190), (97, 191), (97, 192), (98, 188), (98, 189), (98, 190), (98, 191), (98, 192), (98, 193), (98, 194), (99, 169), (99, 171), (99, 189), (99, 190), (99, 191), (99, 192), (99, 193), (99, 194), (100, 190), (100, 191), (100, 192), (100, 193), (100, 194), (100, 195), (101, 169), (101, 171), (101, 187), (101, 191), (101, 192), (101, 193), (101, 194), (101, 195), (101, 196), (102, 192), (102, 193), (102, 194), (102, 195), (102, 196), (102, 197), (103, 193), (103, 194), (103, 195), (103, 196), (103, 197), (103, 198), (104, 194), (104, 195), (104, 196), (104, 197), (104, 198), (104, 199), (105, 195), (105, 196), (105, 197), (105, 198), (105, 199), (105, 200), (106, 196), (106, 197), (106, 199), (106, 200), (106, 201), (107, 197), (107, 198), (107, 199), (107, 200), (107, 201), (107, 202)]
    list_funcs = [
    (main, (bdf_nuclides_list[k], "inc", k%16, "12_01_2024_bdf_inc_output")) for k in range(len(bdf_nuclides_list)) # k%16 in stead of 0 17*90
    ]
    run_parallel(list_funcs)

# mkdir sensitivity
# need output txt to have Z, N (or A) OR just read straight from output xml. jk, just got Y(A)
# rm sensitivity/test_output/ya*
# EDIT ALL REACTIONS! Not just one
# export NNT_USE_SPARSKIT2=1

# for bdf: for i in {0..15}; do cp ./new_bdf_reactions.xml ./A_inc_bdf_$i.xml; done
# CHANGE VAL IN FUTURE IF NEEDED, currently 5 to avoid integration errors hopefully
