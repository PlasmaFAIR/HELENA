import os
import csv
import numpy as np


def get_directories():
    # Define recognized output file data extensions that will be retained in "Dir"
    file_extensions = [".PDT", ".pdt", "PLT", "plt", ".nam", ".dat", ".out"]

    # Obtain home directory and contents
    dir_home = []  # List of all folders in home.
    dir_home_contents = os.listdir(os.path.abspath("."))
    # Determine folders within home directory and add correct 'grammar'.
    for i in range(len(dir_home_contents)):
        if os.path.isdir(dir_home_contents[i]):
            dir_home.append("./" + dir_home_contents[i] + "/")

    directories = []  # List containing all simulation folder directories relative to HELENA
    files = []  # List containing all output file in each Dirlist folder relative to HELENA

    # Extract directories of each sub-folder within home directory
    numfolders = 0
    for i in range(len(dir_home)):
        previous_numfolders = numfolders
        dir_current = dir_home[i]
        dir_contents = os.listdir(dir_current)

        # For each file contained within the subfolders, determine which are datafiles.
        for j in range(len(dir_contents)):
            file_name = dir_contents[j]

            # Save datafiles (with root) to working directory (Dir) and number of datafolders.
            if any([x in file_name for x in file_extensions]):
                files.append(dir_current + file_name)
                if numfolders == previous_numfolders:
                    directories.append(dir_current)
                    numfolders += 1

    # If no folders detected, end analysis script.
    if numfolders == 0:
        print("--------------------------------------------")
        print("No Output Files Detected, Aborting Analysis.")
        print("--------------------------------------------")
        print("")
        exit()

    # Maintain alphanumerical foldername structure (Dirlist) in-sync with dataname structure (Dir)
    return numfolders, sorted(files), sorted(directories), dir_home


def write_to_csv(data, directory, filename, header=None, mode="w"):
    # Takes 1D or 2D array and writes to a datafile in .csv format
    # Inputs,
    # Data = Data to be written, Array (real) :: 1D or 2D
    # Align data array row-wise, 	i.e. [i,j] = [Radius,Height]
    # Origin as per HPEM mesh, 	i.e. [0,0] = top left corner of image
    # Directory = Folder to write to, String
    # Filename = Data file will be named "Filename.csv", String
    # Header = Array of Header information, String
    # Mode = "w" to write new file or "a" to append to existing file
    # Returns,
    # ()
    ###########

    if header is None:
        header = []

    # Write array length and SI dimension to file
    with open(directory + filename, mode) as file:
        # Open Header
        file.write("*START HEADER*")
        file.write("\n")
        # Write supplied header information
        for i in range(0, len(header)):
            file.write(header[i])
            file.write("\n")

        # Close Header
        file.write("*END HEADER*")
        file.write("\n")

        # Append CSV formatted data after header
        writer = csv.writer(file)
        # Write 1D array as a single row
        if np.ndim(data) == 1:
            writer.writerow(data)
        # Write 2D array [i,j] as: 'i rows of j length each'
        elif np.ndim(data) == 2:
            writer.writerows(data)
        else:
            raise ValueError("Data must be 1D or 2D")

    return ()


def read_from_csv(directory, filename, mode="r", cycle=0):
    # Reads a .csv formatted file and returns a 1D or 2D data array and 1D header array
    # Inputs,
    # Directory = Folder to read from. 						[String]
    # Filename = Data file will be named "Filename.csv". 		[String]
    # Mode = "w" to write new file or "a" to append to existing file.
    # Cycle = For time-resolved data, each 2D data array is sequentially
    # written as a new "Cycle". 						[Integer]
    # Returns,
    # Header = 1D array containing each row of header data	[Strings]
    # Data = 1D or 2D array of data following Header			[Floats]
    # Notes,
    # 1D arrays are returned as nested list, i.e. [[DATA]]
    ###########

    # Write array length and SI dimension to file
    with open(directory + filename, mode) as file:
        # Initiate any required lists
        header = []
        data = []

        # Read Header sequentially and identify end index
        raw_data = file.readlines()
        for i in range(0, len(raw_data)):
            # Strip "new line" character (\n) from each entry
            header.append(raw_data[i].strip("\n"))
            # Stop once "*END HEADER*" is reached
            if "*END HEADER*" in raw_data[i]:
                # Data block starts 1 idx after *END HEADER*
                data_start_idx = i + 1
                break

        # Identify end index of data block
        for i in range(data_start_idx, len(raw_data)):
            # Cyclic data blocks are deliminated by next header
            if "*START HEADER*" in raw_data[i]:
                data_end_idx = i
                break
            else:
                # Non-cyclic data continues to end of file
                data_end_idx = len(raw_data)

        # Advance data indices to requested Cycle
        if cycle > 0:
            # Read length of Cycle header - default 3
            for i in range(data_end_idx, len(raw_data)):
                # Stop once "*END HEADER*" is reached
                if "*END HEADER*" in raw_data[i]:
                    cycle_header_idx = i
                    cycle_header_len = cycle_header_idx - data_end_idx + 1
                    break
                else:
                    cycle_header_len = 3

            # Advance data start and end indices to requested cycle
            data_block_len = data_end_idx - data_start_idx
            data_start_idx += cycle * (data_block_len + cycle_header_len)
            data_end_idx += cycle * (data_block_len + cycle_header_len)

        # Read Data and append to output array in row-wise fashion
        for i in range(data_start_idx, data_end_idx):
            # Split each row into scalars, assuming comma delimination
            split_row = raw_data[i].split(",")
            # Convert each scalar from .csv dtype U8 to dtype U16 float
            for j in range(0, len(split_row)):
                split_row[j] = float(split_row[j])

            # Append Row-wise to output data array
            data.append(split_row)

    return data, header
