"""Total script"""

import os
import concurrent.futures
import sys
import subprocess
from iterative_r_process import find_min_error

# helpful commands:
# export NNT_USE_SPARSKIT2=1

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


def edit_xml_with_xmlstarlet(input_file_param, output_file_param, new_value_param):
    """Edits XML file in Python by running xmlstarlet in shell."""
    # Command to update the XML file using xmlstarlet
    # h1 = new_value_param
    # n = str(1-float(new_value_param))

    command = [
        "xmlstarlet", "ed",
        "-u", "/zone_data/zone/optional_properties/property[@name='tend']",
        "-v", new_value_param,
        input_file_param
    ]

    # command = [
    #     "xmlstarlet", "ed",
    #     "-u", "/zone_data/zone/mass_fractions/nuclide[@name='n']/x",
    #     "-v", n,
    #     "-u", "/zone_data/zone/mass_fractions/nuclide[@name='h1']/x", # ADDED FOR YE ONLY
    #     "-v", h1, # ADDED FOR YE ONLY
    #     input_file_param
    # ]

    # Redirect output to the output file
    with open(output_file_param, 'w', encoding='UTF-8') as outfile:
        result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, check=False)

    # Check if the command was successful
    if result.returncode != 0:
        print(f"Error editing XML file: {result.stderr.decode()}")
    else:
        print(f"XML file edited successfully, saved to {output_file_param}")


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
        f"{os.environ['HOME']}/nucnet-tools-code/examples/analysis/print_abundances_vs_nucleon_number",
        output_xml,
        "a",
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


def clean_up(tag):
    """Get rid of extra output XMLs, zone XMLs, and output TXTs."""
    # Types of files
    files_to_clean = [f"xmy_output*{tag}.xml", f"zone*{tag}.xml", f"ya*{tag}.txt"]

    # Change into directory
    os.chdir(f"{os.environ['HOME']}/nucnet-tools-code/ZZZ_my_data")

    for pattern in files_to_clean:
        # Execute shell command to remove files matching the pattern
        command = f"rm {pattern}"
        result = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, \
                                stderr=subprocess.PIPE, check=False)
        if result.returncode != 0:
            print(f"Failed to remove files with pattern '{pattern}': {result.stderr.decode()}")


def single_calc(num, new_zone_value, loc):
    """Performs one calculation on a new zone. Saves it to a new text file."""
    # Define variables and paths
    home = os.environ['HOME']
    net_file = f"{home}/nucnet-tools-code/data_pub/net_NO_fission_combo.xml" # CHANGED TO GET HIGHER T9, was prev my_net
    og_zone_file = f"{home}/nucnet-tools-code/data_pub/zone.xml" # CHANGE THIS BACK LATER
    
    new_zone_xml = f"{home}/nucnet-tools-code/ZZZ_my_data/{loc}/zone{num}.xml"
    new_output_xml = f"{home}/nucnet-tools-code/ZZZ_my_data/{loc}/xmy_output_with_fission{num}.xml"
    x_path = "[z <= 106]" #90 #changed to 106 now for new data!
    output_txt = f"{home}/nucnet-tools-code/ZZZ_my_data/{loc}/ya_with_fission{num}.txt"

    # let user know calculation has begun
    print(f"Calculation {num} has started.")

    # edit and create a new zone xml file
    edit_xml_with_xmlstarlet(og_zone_file, new_zone_xml, new_zone_value)

    # run_single_zone calculation
    run_single_zone(net_file, new_zone_xml, new_output_xml, x_path)

    # make txt files from output xml, which was created by network calc
    txt_files(new_output_xml, output_txt)

    # let user know calculation has ended
    print(f"Calculation {num} has finished.")


# Is called only if the script is run directly
if __name__ == "__main__":
    if sys.argv[1] == "clean":
        clean_up("xxx")

    elif sys.argv[1] == "calc":
        list_new_values = [str(1.e25 * (1000**i)) for i in range(0, 12)]
        list_funcs = [
        (single_calc, (str(j) + "_nf", list_new_values[j], "8_06_2024_tend_v1")) for j in range(12)
        ]
        run_parallel(list_funcs)

    elif sys.argv[1] == "test":
        #print(os.environ.get('MAPP_NR_CPUS'))
        list_new_values = [str(.00003*i) for i in range(1,2)]
        list_funcs = [
        (single_calc, (str(j+5000), list_new_values[j])) for j in range(1)
        ]
        run_parallel(list_funcs)

# starlet checks
    # Check for XML attribute to change?

# single_calc() checks
    # correct net file?
    # correct zone file?
    # correct Z value in xpath?

# main checks
    # made new folder?
    # Are new values correct?
    # file tag correct?
    # ranges correct?
    # Did you export NNT_USE_SPARSKIT2=1
    # Change numbers for filenames?



# clean() checks
    # correct tag?


