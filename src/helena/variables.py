import re

from .utility import string_in_variable


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


def variable_interpolator(
    variable_indices, variable_strings, min_shared_variables, globalnumvars
):
    # Identifies if variable exists in all simulations, rejects if not.
    # Allows for the comparison of datasets with different icp.dat files.
    # Takes variable_indices, variable_strings, globalMinSharedVariables
    # Returns variable_indices and variable_strings with largest commonly shared variables.
    # proclist,varlist = VariableInterpolator(variable_indices,variable_strings,min_shared_variables):

    # No interpolation needed if variable count is the same for all datasets.
    # if all(map(lambda x: x == globalnumvars[0], globalnumvars)):		#Py2.x.x Lambda Method
    if globalnumvars.count(globalnumvars[0]) == len(
        globalnumvars
    ):  # Py3.x.x Count Method
        return variable_indices, variable_strings

    # Identify elements in each variable_strings which are not in the min_shared_variables list
    inter = set(min_shared_variables).symmetric_difference(variable_strings)
    inter = list(inter)

    # If at least one element not present in all folders, remove from lists to be plotted
    if len(inter) != 0:
        for i in range(0, len(inter)):
            j = 0
            while j < len(variable_strings):
                # Check for exact string match, e.g to avoid "AR in AR+".
                if inter[i] == variable_strings[j]:
                    del variable_strings[j]
                    del variable_indices[j]
                else:
                    j += 1

    return variable_indices, variable_strings


def variable_unit_conversion(profile, variable, units="SI", atomic_species=None):
    # Converts units and direction (sign) for input 1D, 2D data arrays.
    # Takes profile and variable name, returns profile in required SI unit.
    # Implicitly calculates for common variables, explicitly for densities.

    if atomic_species is None:
        atomic_species = []

    # For Pressures, convert to mTorr or Pa as requested, or retain as default Torr.
    if string_in_variable(variable, ["PRESSURE"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 133.333333333  # [Pa]
            elif units == "CGS":
                profile[i] = profile[i]  # [Torr]

    # For ionisation rates, convert from [cm3 s-1] to [m3 s-1]
    if string_in_variable(variable, ["S-", "SEB-", "SRCE-"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e6  # [m3 s-1]
            elif units == "CGS":
                profile[i] = profile[i]  # [cm3 s-1]

    # For fluxes, convert from [cm-2] to [m-2]. (also reverse axial flux direction)
    if string_in_variable(variable, ["E FLUX-Z", "E FLUX-R", "FZ-", "FR-"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e4  # [m2 s-1]
            elif units == "CGS":
                profile[i] = profile[i]  # [cm2 s-1]

    if string_in_variable(variable, ["E FLUX-Z", "FZ-"]):
        for i in range(0, len(profile)):
            profile[i] = profile[i] * (-1)

    # For velocities, convert from [cms-1] to [ms-1] or [kms-1]. (also reverse axial velocity)
    if string_in_variable(variable, ["VR-NEUTRAL", "VZ-NEUTRAL"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e-2  # [ms-1]
            elif units == "CGS":
                profile[i] = profile[i]  # [cms-1]

    if string_in_variable(variable, ["VR-ION+", "VZ-ION+", "VR-ION-", "VZ-ION-"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e-5  # [kms-1]
            elif units == "CGS":
                profile[i] = profile[i]  # [cms-1]

    if string_in_variable(variable, ["VZ-NEUTRAL", "VZ-ION+", "VZ-ION-"]):
        for i in range(0, len(profile)):
            profile[i] = profile[i] * (-1)

    # For B-field strengths, convert to Tesla or retain as default Gauss.
    if string_in_variable(variable, ["BR", "BT", "BZ", "BRF", "BTHETA"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i]  # [G]
            elif units == "CGS":
                profile[i] = profile[i]  # [G]

    # ~~~ AXIAL MAGNETIC FIELD IS NOT REVERSED HERE - RM: NEED TO LOOK INTO THIS... ~~~#
    if string_in_variable(variable, ["BZ", "BZS"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * (-1)  # [G]
            elif units == "CGS":
                profile[i] = profile[i]  # [G]

    # For E-field strengths, convert from [V cm-1] to [V m-1]. (also reverse axial field)
    if string_in_variable(
        variable, ["EF-TOT", "EAMB-R", "EAMB-Z", "ERADIAL", "ETHETA", "EAXIAL"]
    ):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 100.0  # [V m-1]
            elif units == "CGS":
                profile[i] = profile[i]  # [V cm-1]

    if string_in_variable(variable, ["EAMB-Z"]):
        for i in range(0, len(profile)):
            profile[i] = profile[i] * (-1)

    # For surface charge, convert from [C cm-3] to [C m-3].
    if string_in_variable(variable, ["RHO"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e6  # [C m-3]
            elif units == "CGS":
                profile[i] = profile[i]  # [C cm-3]

    # For plasma conductivity, convert from [S cm-1] to [S m-1].
    if string_in_variable(variable, ["SIGMA"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e2  # [S m-1]
            elif units == "CGS":
                profile[i] = profile[i]  # [S cm-1]

    # For Current Densities, convert from [A cm-2] to [mA cm-2]. (also reverse axial current)
    if string_in_variable(variable, ["JZ-NET", "JR-NET", "J-THETA", "J-TH(MAG)"]):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] / 1000.0  # [A m-2]
            elif units == "CGS":
                profile[i] = profile[i]  # [A cm-2]

    if string_in_variable(variable, ["JZ-NET"]):
        for i in range(0, len(profile)):
            profile[i] = profile[i] * (-1)

    # For power densities, convert from [Wcm-3] to [Wm-3].
    if string_in_variable(
        variable,
        ["POW-ALL", "POWALL", "POW-TOT", "POW-ICP", "POWICP", "POW-RF", "POW-RF-E"],
    ):
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e6  # [W m-3]
            elif units == "CGS":
                profile[i] = profile[i]  # [W cm-3]

    # For densities, convert from [cm-3] to [m-3]. (atomic_species is defined in icp.nam input)
    if variable in atomic_species or variable in [
        x.replace("^", "+") for x in atomic_species
    ]:
        for i in range(0, len(profile)):
            if units == "SI":
                profile[i] = profile[i] * 1.0e6  # [m-3]
            elif units == "CGS":
                profile[i] = profile[i]  # [cm-3]

    return profile


def variable_label_maker(variable_strings, units, image_logplot, atomic_species):
    # Makeshift way of creating units for each legend entry.

    # Define common lists for implicit legend generation.
    power_list = ["POW-ALL", "POW-TOT", "POW-ICP", "POW-RF", "POW-RF-E"]
    flux_list = ["FZ-", "FR-", "EFLUX-R", "EFLUX-Z"]
    ionisation_list = ["S-", "SEB-"]
    velocity_list = ["VZ-", "VR-"]

    # Define Regular Expression lists for numericised ICP coil set variable names
    reg_ex = re.compile("POWICP.")
    POWICP_vars = ["POWICP"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]

    reg_ex = re.compile("ERADIAL.")
    ERADIAL_vars = ["ERADIAL"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("ETHETA.")
    ETHETA_vars = ["ETHETA"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("EAXIAL.")
    EAXIAL_vars = ["EAXIAL"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]

    reg_ex = re.compile("PHASEER.")
    PHASEER_vars = ["PHASEER"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("PHASE.")
    PHASE_vars = ["PHASE"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("PHASEEZ.")
    PHASEEZ_vars = ["PHASEEZ"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]

    reg_ex = re.compile("BR.")
    BR_vars = ["BR"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("BT.")
    BT_vars = ["BT"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("BZ.")
    BZ_vars = ["BZ"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("BRF.")
    BRF_vars = ["BRF"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]

    reg_ex = re.compile("PHASEBR.")
    PHASEBR_vars = ["PHASEBR"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("PHASEBT.")
    PHASEBT_vars = ["PHASEBT"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("PHASEBZ.")
    PHASEBZ_vars = ["PHASEBZ"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]
    reg_ex = re.compile("J-THETA.")
    JTHETA_vars = ["J-THETA"] + [
        string for string in variable_strings if re.match(reg_ex, string)
    ]

    log_string = "Log$_{10}$" if image_logplot else ""

    variable_legends = []
    for i in range(0, len(variable_strings)):
        # Explicit Pressure and Species Densities.
        if variable_strings[i] in ["PRESSURE"]:
            variable = "Pressure"
            if units == "SI":
                variable_unit = log_string + "[Pa]"
            elif units == "CGS":
                variable_unit = log_string + "[Torr]"
        elif variable_strings[i] in ["AR", "AR3S"]:
            variable = "Neutral Ar Density"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$]"
        elif variable_strings[i] in ["E"]:
            variable = "Electron Density n$_{e}$"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$]"
        elif variable_strings[i] in ["AR+"]:
            variable = "Ar+ Density"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$]"

        # Explicit Ionization Rates.
        elif variable_strings[i] == "S-E":
            variable = "Bulk e$^-$ Source Rate \n"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"
        elif variable_strings[i] == "SEB-E":
            variable = "Secondry e$^-$ Source Rate \n"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"
        elif variable_strings[i] == "EB-ESORC":
            variable = "Secondry e$^-$ Relaxation Rate \n"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"
        elif variable_strings[i] == "S-AR+":
            variable = "Bulk Ar+ Ionization Rate \n"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"
        elif variable_strings[i] == "SEB-AR+":
            variable = "Secondry Ar+ Ionization Rate \n"
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"

        # Explicit Vibrational States.
        elif variable_strings[i] == "GSH2V1":
            variable = "1st Vibrational Excited State \n"
            variable_unit = log_string + "[cm$^{-3}$]"
        elif variable_strings[i] == "GSH2V4":
            variable = "4th Vibrational Excited State \n"
            variable_unit = log_string + "[cm$^{-3}$]"
        elif variable_strings[i] == "GSH2V14":
            variable = "14th Vibrational Excited State \n"
            variable_unit = log_string + "[cm$^{-3}$]"

        # Explicit Species Temperatures.
        elif variable_strings[i] == "TE":
            variable = "Electron Temperature T$_{e}$"
            variable_unit = log_string + "[eV]"
        elif variable_strings[i] == "TG-AVE":
            variable = "Neutral Gas Temperature"
            variable_unit = log_string + "[K]"

        # Explicit Species Velocities.
        elif variable_strings[i] == "VZ-NEUTRAL":
            variable = "Neutral Axial Velocity"
            if units == "SI":
                variable_unit = log_string + "[ms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"
        elif variable_strings[i] == "VR-NEUTRAL":
            variable = "Neutral Radial Velocity"
            if units == "SI":
                variable_unit = log_string + "[ms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"
        elif variable_strings[i] == "VZ-ION+":
            variable = "+Ion Axial Velocity"
            if units == "SI":
                variable_unit = log_string + "[kms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"
        elif variable_strings[i] == "VR-ION+":
            variable = "+Ion Radial Velocity"
            if units == "SI":
                variable_unit = log_string + "[kms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"
        elif variable_strings[i] == "VZ-ION-":
            variable = "-Ion Axial Velocity"
            if units == "SI":
                variable_unit = log_string + "[kms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"
        elif variable_strings[i] == "VR-ION-":
            variable = "-Ion Radial Velocity"
            if units == "SI":
                variable_unit = log_string + "[kms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"

        # Explicit Species Fluxes.
        elif variable_strings[i] == "E FLUX-Z":
            variable = "Electron Axial Flux"
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"
        elif variable_strings[i] == "E FLUX-R":
            variable = "Electron Radial Flux"
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"
        elif variable_strings[i] == "FZ-AR+":
            variable = "Ar+ Axial Flux"
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"
        elif variable_strings[i] == "FR-AR+":
            variable = "Ar+ Radial Flux"
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"
        elif variable_strings[i] == "FZ-AR":
            variable = "Ar Axial Flux"
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"
        elif variable_strings[i] == "FR-AR":
            variable = "Ar Radial Flux"
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"

        # Explicit Electrodynamic Properties
        elif variable_strings[i] in ["RHO"]:
            variable = "Charge Density $\\rho$"
            if units == "SI":
                variable_unit = log_string + "[C m$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[C cm$^{-3}$]"
        elif variable_strings[i] in ["PPOT", "P-POT"]:
            variable = "Plasma Potential V$_{p}$"
            variable_unit = log_string + "[V]"
        elif variable_strings[i] in ["SIGMA"]:
            variable = "Conductivity $\\sigma$"
            if units == "SI":
                variable_unit = log_string + "[S m$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[S cm$^{-1}$]"

        elif variable_strings[i] in ["EF-TOT"]:
            variable = "Absolute E-Field Amplitude"
            if units == "SI":
                variable_unit = log_string + "[Vm$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Vcm$^{-1}$]"
        elif variable_strings[i] in ERADIAL_vars + ["ER", "EAMB-R"]:
            variable = "Radial E-Field Amplitude $E_{R}$"
            if units == "SI":
                variable_unit = log_string + "[Vm$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Vcm$^{-1}$]"
        elif variable_strings[i] in ETHETA_vars + ["ET", "EAMB-T"]:
            variable = "Azimuthal E-Field Amplitude $E_{\\theta}$"
            if units == "SI":
                variable_unit = log_string + "[Vm$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Vcm$^{-1}$]"
        elif variable_strings[i] in EAXIAL_vars + ["EZ", "EAMB-Z"]:
            variable = "Axial E-Field Amplitude $E_{Z}$"
            if units == "SI":
                variable_unit = log_string + "[Vm$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Vcm$^{-1}$]"
        elif variable_strings[i] in PHASEER_vars:
            variable = "Radial E-Field Phase"
            variable_unit = "[Radians]"
        elif variable_strings[i] in PHASE_vars:  # RM SJD, CATCHES ALL PHASE VARIABLES
            variable = "Azimuthal E-Field Phase"
            variable_unit = "[Radians]"
        elif variable_strings[i] in PHASEEZ_vars:  # RM SJD, CATCHES ALL PHASE VARIABLES
            variable = "Axial E-Field Phase"
            variable_unit = "[Radians]"

        elif variable_strings[i] == "BRS":
            variable = "Radial Static \n B-field Magnitude $B_{R}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"
        elif variable_strings[i] == "BTS":
            variable = "Azimuthal Static \n B-field Magnitude $B_{\\theta}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"
        elif variable_strings[i] == "BZS":
            variable = "Axial Static \n B-field Magnitude $B_{Z}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"

        elif variable_strings[i] in BR_vars:
            variable = "Radial Induced \n B-field Magnitude $B_{R}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"
        elif variable_strings[i] in BT_vars + ["BTHETA"]:
            variable = "Azimuthal Induced \n B-field Magnitude $B_{\\theta}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"
        elif variable_strings[i] in BZ_vars:
            variable = "Axial Induced \n B-field Magnitude $B_{Z}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"
        elif variable_strings[i] in BRF_vars:
            variable = "Induced B-field Magnitude $B_{RF}$"
            if units == "SI":
                variable_unit = log_string + "[G]"
            elif units == "CGS":
                variable_unit = log_string + "[G]"
        elif variable_strings[i] in PHASEBR_vars:
            variable = "Radial Induced B-field Phase"
            variable_unit = "[Radians]"
        elif variable_strings[i] in PHASEBT_vars:
            variable = "Azimuthal Induced B-field Phase"
            variable_unit = "[Radians]"
        elif variable_strings[i] in PHASEBZ_vars:
            variable = "Axial Induced B-field Phase"
            variable_unit = "[Radians]"

        elif variable_strings[i] in ["JZ-NET"]:
            variable = "Axial Net Current Density $J_{Z}$"
            if units == "SI":
                variable_unit = log_string + "[A m$^{-2}$]"
            elif units == "CGS":
                variable_unit = log_string + "[A cm$^{-2}$]"
        elif variable_strings[i] in ["JR-NET"]:
            variable = "Radial Net Current Density $J_{R}$"
            if units == "SI":
                variable_unit = log_string + "[A m$^{-2}$]"
            elif units == "CGS":
                variable_unit = log_string + "[A cm$^{-2}$]"
        elif variable_strings[i] in JTHETA_vars:
            variable = "Azimuthal Net Current Density $J_{\\theta}$"
            if units == "SI":
                variable_unit = log_string + "[A m$^{-2}$]"
            elif units == "CGS":
                variable_unit = log_string + "[A cm$^{-2}$]"
        elif variable_strings[i] in ["J-TH(MAG)"]:
            variable = "Azimuthal MCS Electron Current Density $J_{e\\theta}$"
            if units == "SI":
                variable_unit = log_string + "[A m$^{-2}$]"
            elif units == "CGS":
                variable_unit = log_string + "[A cm$^{-2}$]"
        elif variable_strings[i] in ["J-TH(PHA)"]:
            variable = "Azimuthal MCS Electron \n Current Phase $J_{e\\theta}$"
            variable_unit = "[Radians]"

        # Explicit Power Deposition.
        elif variable_strings[i] in ["POW-TOT"]:
            variable = "Total Power Density"
            if units == "SI":
                variable_unit = log_string + "[Wm$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Wcm$^{-3}$]"
        elif variable_strings[i] in POWICP_vars + [
            "POW-ICP"
        ]:  # Note: 	POW-ICP:  total icp power
            variable = "Inductive Power Density"  # POWICP-n: coilset #n icp power
            if units == "SI":
                variable_unit = log_string + "[Wm$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Wcm$^{-3}$]"
        elif variable_strings[i] in ["POW-RF"]:
            variable = "RF Power Density"
            if units == "SI":
                variable_unit = log_string + "[Wm$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Wcm$^{-3}$]"
        elif variable_strings[i] in ["POW-RF-E"]:
            variable = "RF Electron Power Density"
            if units == "SI":
                variable_unit = log_string + "[Wm$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Wcm$^{-3}$]"

        # Explicit Collision Rates.
        elif variable_strings[i] == "COLF":
            variable = "Electron Collision Frequency"
            variable_unit = log_string + "[s$^{-1}$]"

        # Implicit Variables.
        elif string_in_variable(variable_strings[i], ionisation_list):
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"
        elif string_in_variable(variable_strings[i], ["SRCE-"]):
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$s$^{-1}$]"
        elif string_in_variable(variable_strings[i], ["T-"]):
            variable = variable_strings[i]
            variable_unit = "[K]"
        elif string_in_variable(variable_strings[i], velocity_list):
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[kms$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cms$^{-1}$]"
        elif string_in_variable(variable_strings[i], flux_list):
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[m$^{-2}$s$^{-1}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-2}$s$^{-1}$]"
        elif string_in_variable(variable_strings[i], ["POW-"]):
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[Wm$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[Wcm$^{-3}$]"
        elif variable_strings[i] in [x.replace("^", "+") for x in atomic_species]:
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$]"
        elif variable_strings[i] in atomic_species:
            variable = variable_strings[i]
            if units == "SI":
                variable_unit = log_string + "[m$^{-3}$]"
            elif units == "CGS":
                variable_unit = log_string + "[cm$^{-3}$]"

        # Default if no fitting variable found.
        else:
            variable = "Variable"
            variable_unit = log_string + "[Unit]"

        variable_legends.append(variable + " " + variable_unit)

    return variable_legends


def azimuthal_phase_conversion(profile, variable):
    # Applies field phase to sign of field amplitude for azimuthal data arrays.
    # Takes 1D, or 2D data profile and variable string.
    # Returns data array multiplied by sin(phase).
    # Returns non-azimuthal data arrays unchanged.

    return profile

    # FIXME

    # Global toggle to enforce plotting of magnitudes only if requested
    # Only Azimuthally varying fields require phase conversion
    # if ConvAzimuthalPhase == False:
    #    return (profile)
    # elif string_in_variable(variable, ['ETHETA']) == True:
    #    phaseprocess, phasevariable = enumerate_variables(['PHASE'], Header_TEC2D[l])
    # elif string_in_variable(variable, ['J-THETA']) == True:
    #    phaseprocess, phasevariable = enumerate_variables(['PHASE'], Header_TEC2D[l])
    # elif string_in_variable(variable, ['J-TH(MAG)']) == True:
    #    phaseprocess, phasevariable = enumerate_variables(['J-TH(PHA)'], Header_TEC2D[l])
    # else:
    #    return (profile)

    # Extract the appropriate phase data for the supplied variable
    # phasemap = ImageExtractor2D(Data[l][phaseprocess[0]], phasevariable[0])

    # Convert from phase [0 --> 2pi] to relative azimuthal direction (0 --> -1 --> +1)
    # for i in range(0, len(phasemap)):
    #    for j in range(0, len(phasemap[i])):
    #        phasemap[i][j] = np.sin(phasemap[i][j])

    # Multiply azimuthal data amplitude by sin( azimuthal phase )
    # if len(profile.shape) == 2:
    #    profile = profile * phasemap  # 2D data profile
    # elif len(profile.shape) == 1:
    #    profile = profile * phasemap.flatten()  # 1D data profile
    # !!! RM SJD, Not Tested 1D Yet!

    # return profile
