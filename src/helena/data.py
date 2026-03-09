from .variables import (
    variable_unit_conversion,
    enumerate_variables,
    variable_interpolator,
)


def extract_raw_data(directory, name_string, list_index):
    # Takes full directory list (directory) and data filename type (e.g. .png, .txt)
    # Returns row-wise list of data and length of datafile.
    # data_raw, data_length = extract_raw_data(directory,'.dat',l)

    try:
        data_file_dir = filter(lambda x: name_string in x, directory)
        data_file_dir = sorted(data_file_dir)
        data_raw = open(data_file_dir[list_index]).readlines()
        nn_data = len(data_raw)
    except:
        print("Unable to extract " + str(name_string))
        exit()

    return data_raw, nn_data


def read_TEC2D(data_raw, header, numvariables, Rmesh, Zmesh, offset=0, dimension="2D"):
    # Takes ASCII data in 2/3D format and converts to HELENA friendly structure.
    # Requires rawdata(2D/3D), header and variable number and mesh dimensions.
    # Allows for an optional offset in 'starting' variable number.
    # Returns 2D array of form [Variables,datapoint(R,Z)]
    # current_folder_data = ReadTEC2D(rawdata_2D[l],header_2D,numvariables_2D)

    # Excluding the header, split each row of data and append items to 1D list.
    current_folder_data, data_array1_d = [], []

    for i in range(header, len(data_raw)):
        # If end of phasecycle reached, break. (Applicable to 3D Datafiles only)
        if "CYCL=" in data_raw[i] or "ITER=" in data_raw[i]:
            break
        else:
            current_row = data_raw[i].split()

        # For all elements in the current row, convert to float and save in list.
        for j in range(0, len(current_row)):
            try:
                data_array1_d.append(float(current_row[j]))
            except:
                String_Conversion_Error = 1

    # If data is 1D, seperate into 1D chunks using Zmesh as the chunk size
    if dimension == "1D":
        # Seperate total 1D array into further 1D sub-arrays with data for each variable.
        for i in range(offset, numvariables):
            numstart = Zmesh * i
            numend = Zmesh * (i + 1)
            current_folder_data.append(list(data_array1_d[numstart:numend]))

        return current_folder_data

    # If data is 2D, return array in 2D chunks using Zmesh*Rmesh as the chunk size
    elif dimension == "2D":
        # Seperate total 1D array into 2D array with data for each variable.
        # Offset data by a certain number of variable 'chunks' if requested.
        for i in range(offset, numvariables):
            numstart = int((Zmesh * Rmesh) * i)
            numend = int((Zmesh * Rmesh) * (i + 1))
            current_folder_data.append(list(data_array1_d[numstart:numend]))

        return current_folder_data


def read_TEC2D_phase(
    folder,
    variables,
    directory,
    movie1,
    Header_movie1,
    MinSharedVariables,
    Globalnumvars,
    R_mesh,
    Z_mesh,
    Units,
    AtomicSpecies,
):
    # Extracts all phase and variable data for the provided folder ID.
    # Initial R and Z for CYCL=1 are skipped over and not saved.
    # Takes current folder, returns Data[phase][variable][datapoints,R/Z]
    # Data,Phaselist = ReadTEC2DPhase(folder=l,variables=PhaseVariables)

    # Load data from movie1 file and unpack into 1D array.
    try:
        rawdata, filelength = extract_raw_data(
            directory, movie1[folder].split("/")[-1], folder
        )
    except:
        rawdata, filelength = extract_raw_data(directory, "movie1.pdt", folder)

    # Read through all variables for each file and stop when list ends.
    # Movie1 has geometry at top, therefore len(header) != len(variables).
    # Only the first encountered geometry is used to define variable zone.
    VariableEndMarker, HeaderEndMarker = "GEOMETRY", "ZONE"
    Mov1VariableStrings, numvar = [], 0
    for i in range(2, filelength):
        if HeaderEndMarker in str(rawdata[i]):
            header = i + 2  # plus 2 to skip to first data line.
            break
        if VariableEndMarker in str(rawdata[i]) and numvar == 0:
            numvar = i - 1 - 2  # minus 1 for overshoot, minus 2 for starting at 2.
        if len(rawdata[i]) > 1 and numvar == 0:
            Mov1VariableStrings.append(str(rawdata[i][:-2].strip(' \t\n\r"')))

    # Enumerate variables in the order they appear in movie1.pdt
    proclist, varlist = enumerate_variables(variables, Header_movie1[folder])

    # Interpolate variable lists for all folders to find minimum shared set
    proclist, varlist = variable_interpolator(
        proclist, varlist, MinSharedVariables, Globalnumvars
    )

    # Rough method of obtaining the movie1.pdt cycle locations for data extraction.
    cycleloc = []
    for i in range(0, len(rawdata)):
        if "CYCL=" in rawdata[i]:
            cycleloc.append(i + 1)

    # Cycle through all phases for current datafile, appending per cycle.
    # Variables R and Z only saved for first iteration, they are skipped if i == 0.
    FolderData, Phaselist = [], []
    for i in range(0, len(cycleloc) - 1):  # len(cycleloc)-1 as python starts at idx 0 #
        # R,Z arrays are saved only for first "Cycle", apply +2 variable index offset to ignore
        if i == 0:
            PhaseData = read_TEC2D(
                rawdata,
                cycleloc[i],
                numvar + 2,
                R_mesh[folder],
                Z_mesh[folder],
                offset=2,
            )
            # Convert data from CGS (HPEM DEFAULT) to user requested unit system
            for j in range(0, len(Mov1VariableStrings)):
                PhaseData[j] = variable_unit_conversion(
                    PhaseData[j], Mov1VariableStrings[j], Units, AtomicSpecies
                )
            # PhaseData[j] = AzimuthalPhaseConversion(PhaseData[j],Mov1VariableStrings[j])

            FolderData.append(PhaseData[0:numvar])

        # Later cycles do not save R,Z arrays so no variable index offset is required.
        elif i > 0:
            PhaseData = read_TEC2D(
                rawdata, cycleloc[i], numvar, R_mesh[folder], Z_mesh[folder]
            )
            for j in range(0, len(Mov1VariableStrings)):
                PhaseData[j] = variable_unit_conversion(
                    PhaseData[j], Mov1VariableStrings[j], Units, AtomicSpecies
                )
            # PhaseData[j] = AzimuthalPhaseConversion(PhaseData[j],Mov1VariableStrings[j])

            FolderData.append(PhaseData)

        # Always update phaselist with cycle index
        Phaselist.append("CYCL = " + str(i + 1))

    return FolderData, Phaselist, proclist, varlist
