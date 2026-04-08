def string_in_variable(variable, stringarray):
    # Takes array of strings and compares to variable string.
    # Returns true if any element of stringarray is in variable.

    # Check if each element of string is inside variable.
    for i in range(0, len(stringarray)):
        if stringarray[i] in variable:
            return True

    return False


def is_radial_variable(variable):
    # Takes variablenames and checks if variable is radial.
    # Returns boolian if variable is radial, used for symmetry options.

    return string_in_variable(variable, stringarray=["VR-", "JR-", "FR-", "FLUX-R"])
