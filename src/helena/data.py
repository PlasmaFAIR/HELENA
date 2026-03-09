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
