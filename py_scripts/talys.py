import subprocess
import os
import sys

def run_talys(my_input):
    
    # Change into directory of our fission folder
    os.chdir(f"{os.environ['HOME']}/talys/my_fission_folder")

    # Define paths
    talys_path = f"{os.environ['HOME']}/talys/bin/talys"
    output_file = "talys.out"

    # Open the input and output files
    with open(my_input, 'r') as infile, open(output_file, 'w') as outfile:
        # Run the command with input and output redirection
        subprocess.run([talys_path], stdin=infile, stdout=outfile)

    # move the one file we care about to keep directory
    move = [
        "mv",
        "fission.tot",
        f"{os.environ['HOME']}/talys/my_fission_folder/keep"
    ]
    subprocess.run(move)


def clean_up():
    """Get rid of excess."""

    # Types of files
    files_to_clean = ["energies", "incident.*", "*.tot", "*.fis", "*.dat", "*.inp", "ecis*", "massinfo.out"]

    # Change into directory
    os.chdir(f"{os.environ['HOME']}/talys/my_fission_folder")

    for pattern in files_to_clean:
        # Execute shell command to remove files matching the pattern
        command = f"rm {pattern}"
        result = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, \
                                stderr=subprocess.PIPE, check=False)
        if result.returncode != 0:
            print(f"Failed to remove files with pattern '{pattern}': {result.stderr.decode()}")

if __name__ == "__main__":
    #run_talys("fission_input.txt")
    clean_up()