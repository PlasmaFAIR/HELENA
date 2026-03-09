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
