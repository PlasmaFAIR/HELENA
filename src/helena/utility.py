def string_in_variable(variable, stringarray):
    # Takes array of strings and compares to variable string.
    # Returns true if any element of stringarray is in variable.

    # Check if each element of string is inside variable.
    for i in range(0, len(stringarray)):
        if stringarray[i] in variable:
            return True

    return False
