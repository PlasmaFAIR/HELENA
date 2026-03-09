def enumerate_variables(variables, header):
    # Takes:
    # variables 	- 1D array of strings
    # header 		- 1D Datafile Header Array (for 1 simulation folder)
    # Returns:
    # indices 		- 1D array of indices corresponding to variables in header
    # strings		- 1D array of strings corresponding to variable string in header
    #
    # Example: enumerate_variables(variables, header_TEC2D[l])

    indices, strings = [], []

    for variable in variables:
        if variable in header:
            index = header.index(variable)
            indices.append(index)
            strings.append(header[index])

    return indices, strings


def enumerate_vectors(variable, header, prefixes=None):
    # Takes:
    # variable 		- Variable string, e.g. "AR3S"
    # header 		- 1D Datafile Header Array (for 1 simulation folder)
    # prefixes		- 1D array containing vector prefixes
    # prefixes are limited to ['FR','FZ'] or ['VR','VZ']
    # Returns:
    # radial 		- 1D array, name and index of radial vector for input 'variable'
    # axial 		- 1D array, name and index of axial vector for input 'variable'
    #
    # Example:
    # enumerate_vectors('AR3S', Header_TEC2D[l], prefixes['FR','FZ'])
    # radial = ['FR-AR3S',121]
    # axial = ['FZ-AR3S',120]
    # Where, 120,121 are the TECPLOT2D header indices of the named variables

    if prefixes is None:
        prefixes = ["FR", "FZ"]

    # Create header lookup
    header_index = {name: i for i, name in enumerate(header)}

    # =====#=====# SPECIAL CASES #=====#=====#

    # Many variables do not follow a standard prefix format,
    # this section explicitly defines these special cases.
    #
    # Special variables are directly "matched" to their
    # vector equivelants, here "Wanted" variables are
    # the explicit radial & axial vectors of "Variable"

    if variable in ["E FLUX-R", "E FLUX-Z"]:
        # Prefixes are the explicit vector quantities relating to "Variable"
        vectors = ["E FLUX-R", "E FLUX-Z"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        # Match vector variable names
        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        # Index vector variable names
        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        # Return vector variable names and indices
        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["FR-E", "FZ-E"]:
        # Prefixes are the explicit vector quantities relating to "Variable"
        vectors = ["FR-E", "FZ-E"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        # Match vector variable names
        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        # Index vector variable names
        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        # Return vector variable names and indices
        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["EF-TOT", "EAMB-R", "EAMB-Z"]:
        vectors = ["EAMB-R", "EAMB-Z"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    # NOTE, THIS IS FOR STATIC FIELDS ONLY, DO NOT MULTIPLY BY PHASE
    if variable in ["BT", "BR", "BZ"]:
        vectors = ["BR", "BZ"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["JR-NET", "JZ-NET"]:
        vectors = ["JR-NET", "JZ-NET"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["VR-NEUTRAL", "VZ-NEUTRAL"]:
        vectors = ["VR-NEUTRAL", "VZ-NEUTRAL"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["VR-ION+", "VZ-ION+"]:
        vectors = ["VR-ION+", "VZ-ION+"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["VR-ION-", "VZ-ION-"]:
        vectors = ["VR-ION-", "VZ-ION-"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====#

    if variable in ["I MFLUX-R", "I MFLUX-Z"]:
        vectors = ["I MFLUX-R", "I MFLUX-Z"]
        radial_wanted = [f"{vectors[0]}"]
        axial_wanted = [f"{vectors[1]}"]

        radial_match = next((w for w in radial_wanted if w in header_index), None)
        axial_match = next((w for w in axial_wanted if w in header_index), None)

        radial_index = header_index[radial_match]
        axial_index = header_index[axial_match]

        radial = [radial_match, radial_index]
        axial = [axial_match, axial_index]

        return radial, axial

    # =====#=====# GENERAL CASE #=====#=====#

    # The general case applies to species fluxes and velocities
    # where variables are saved in the following formats:
    # "FR-AR3S", "FZ-AR3S"
    # "VR-AR3S", "VR-AR3S"
    # This general approach applies to any given species,
    # e.g. "FR-O2" or "FZ-NH2"

    # IF PLOTTING VECTORS ONTO DENSITY FIGURES
    # Apply prefix to full name,
    # i.e. "AR3S" becomes "FR-AR3S" 		and IS matched
    #        "FR-AR3S" becomes "FR-FR-AR3S" and IS NOT matched
    # radial_wanted = [f"{Prefixes[0]}-{Variable}"]
    # axial_wanted = [f"{Prefixes[1]}-{Variable}"]

    # IF PLOTTING VECTORS ONTO FLUX FIGURES
    # Apply prefix ignoring first three characters,
    # i.e. "AR3S" becomes "FR-S" 			and IS NOT matched
    #        "FR-AR3S" becomes "FR-AR3S" 	and IS matched
    radial_wanted = [f"{prefixes[0]}-{variable[3::]}"]
    axial_wanted = [f"{prefixes[1]}-{variable[3::]}"]

    # Search .PDT file header and match for "Wanted" string, return None if none exist
    radial_match = next((w for w in radial_wanted if w in header_index), None)
    axial_match = next((w for w in axial_wanted if w in header_index), None)

    # Index variables (in Header order) if a match is found
    if radial_match is None:
        radial_index = None
    else:
        radial_index = header_index[radial_match]
    if axial_match is None:
        axial_index = None
    else:
        axial_index = header_index[axial_match]

    # Return vector variable names and indices
    radial = [radial_match, radial_index]
    axial = [axial_match, axial_index]

    return radial, axial
