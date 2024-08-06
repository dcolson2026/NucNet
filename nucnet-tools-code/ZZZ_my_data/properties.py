"""Get T9(t) and Ye(t)"""

import subprocess
import os
import numpy as np

def props_file(output_xml, output_Ye="Ye_test.txt", output_t_T9="t_T9_test.txt"):
    """Makes txt file from output xml. Can later be analyzed."""
    #xpath1 = '"" z 1 ""'

    # Command for printing mass num and abundances
    command1 = [
        f"{os.environ['HOME']}/nucnet-tools-code/examples/analysis/compute_abundance_moment_in_zones",
        output_xml,
        "",
        "z",
        "1",
        "",
    ]
    command2 = [
        f"{os.environ['HOME']}/nucnet-tools-code/examples/analysis/print_properties",
        output_xml,
        "time",
        "t9",
    ]

    # Redirect output to the output file
    with open(output_Ye, 'w', encoding='ascii') as outfile:
        result1 = subprocess.run(command1, stdout=outfile, stderr=subprocess.PIPE, check=False)
    with open(output_t_T9, 'w', encoding='ascii') as outfile:
        result2 = subprocess.run(command2, stdout=outfile, stderr=subprocess.PIPE, check=False)

    # Check if the command was successful
    if result1.returncode != 0:
        print(f"Error making TXT file: {result1.stderr.decode()}")
    else:
        print(f"TXT file made successfully, saved to {output_Ye}")
    if result2.returncode != 0:
        print(f"Error making TXT file: {result2.stderr.decode()}")
    else:
        print(f"TXT file made successfully, saved to {output_t_T9}")
    
    # load data
    DATA1=np.loadtxt(output_Ye)
    Ye_data = DATA1[:, 1]
    
    DATA2=np.loadtxt(output_t_T9)
    t_data = DATA2[:, 1]
    T9_data = DATA2[:, 2]

    return Ye_data, t_data, T9_data


if __name__ == "__main__":
    props_file("my_output_with_fission0_higher_tend.xml", "Ye_test.txt", "t_T9_test.txt")
