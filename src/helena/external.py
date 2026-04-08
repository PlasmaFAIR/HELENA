import os
from subprocess import Popen, PIPE

import numpy as np


def auto_conv_prof(
    conv_prof_exe, direc, IMATSTATS, IPCMCSPEC, args=None, dir_additions=None
):
    # Runs requested dataconversion script with pre-defined arguments.
    # Takes name of convert script, any predefined arguments and newly created files.
    # Returns nothing, runs script in each folder, expects to run over all folders.
    # *** TECPLOT FILE IS "NPROFILE_TEC2D.PDT" FOR NEUTRAL PROFILE DATA. ***
    # *** TECPLOT FILE IS "IPROFILE_TEC2D.PDT FOR ION PROFILE DATA. ***

    if args is None:
        args = []

    if dir_additions is None:
        dir_additions = []

    home_dir = os.getcwd()
    os.chdir(direc)

    # ENTER FILE NAME FOR RAW PLOTTING DATA (simulation pcmc output file, default pcmc.prof)
    pcmcprof = "pcmc.prof"
    # ENTER TITLE FOR DATA TO BE USED BY TECPLOT (TITLE MUST BE <= 20 CHARACTERS)
    tec_plot_title = "title"
    # ENTER 1 TO "APPROVE" WRITING VARIABLES 0 NOT TO "APPROVE"
    write_vars = "1"
    # Should angular statistics be normalized to f(theta)/solid-angle? (default is f(theta) x d(omega))
    normalise_angular_statistics = "1"
    # Should profiles be averaged across 0 degrees?
    average_across0_degrees = "1"
    # [pcmc.prof, title, write, normalise, average] are basic requirements
    init_args = [
        pcmcprof,
        tec_plot_title,
        write_vars,
        normalise_angular_statistics,
        average_across0_degrees,
    ]

    # Write data for each CPCMCSPEC species, for each CMATSTATS material
    # Asks for material first, then for each species seperately, total = IPCMCSPEC + IMATSTATS*IPCMCSPEC
    # EXAMPLE TEXT: SHOULD DATA FOR POSITION  0.25 CM, MATERIAL 5 BE WRITTEN? (1 OR 0)
    spec_args = np.ones(IMATSTATS * IPCMCSPEC + IPCMCSPEC, dtype=str).tolist()

    # Should flux(energy) be integrated for angle and flux(angle) be integrated for energy
    integrate_flux_energy = "1"
    integrate_args = [integrate_flux_energy]

    # Concat all pcmc.prof args into a single list of strings
    args = init_args + spec_args + integrate_args
    # args = [] 						# <<< UNCOMMENT FOR MANUAL ENTRY OF VALUES.
    # print(args)						# UPDATE FEED ARGS IN NAM.READIN AND REMOVE THE ABOVE LINE ONCE FIXED.
    # print(os.linesep.join(args))	# PRINTS EACH ITEM IN LIST ON SEPERATE LINE
    # #[pcmc.prof,title,1,1,1]+[pcmcmat+pcmcmat*pcmcspecies]+[1] for manual usage

    # =====#=====#

    # Remove old files to avoid any overwrite errors.
    for i in range(0, len(dir_additions)):
        os.system("rm -f " + dir_additions[i])

    # THIS SECTION IS OUT OF DATE, SUBPROCESS.COMMUNICATE REQUIRES "bytes-like objects" NOT STRINGS
    # TypeError: a bytes-like object is required, not 'str'
    # Use predefined arguments if supplied, suppresses output to devnull.
    if len(args) > 0:
        with open(os.devnull, "w") as fp:
            try:
                subprocess = Popen(
                    conv_prof_exe, stdin=PIPE, stdout=fp, encoding="utf8"
                )  # noshell=True
                subprocess.communicate(
                    os.linesep.join(args), timeout=15
                )  # 15 second timeout
            except:
                print("")
                print("### Timeout while using conv_prof.exe   ###")
                print('### Check "args" length in AutoConvProf ###')
                print("")
                exit()

    # If no arguments are supplied, run script and allow user inputs.
    elif len(args) == 0:
        os.system(conv_prof_exe)

    # Return to Home directory (where HELENA3.py is located)
    os.chdir(home_dir)
