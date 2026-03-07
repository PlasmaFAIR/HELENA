import os


def get_directories(dir_home, file_extensions):
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
    return numfolders, sorted(files), sorted(directories)
