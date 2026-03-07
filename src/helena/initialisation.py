import re

import numpy as np


def get_mesh_and_SI(
    numfolders,
    Dirlist,
    icpnam,
    icpdat,
    TEC2D,
    mesh,
    Magmesh,
    image_plotsymmetry,
    IDEBUG,
):
    R_mesh, Z_mesh = [], []
    ISYMlist, IXZlist = [], []
    Radius, Height, Depth = [], [], []
    dr, dz, dy = [], [], []
    VRFM, VRFM2 = [], []
    FREQM, FREQM2 = [], []
    FREQC, FREQGLOB, FREQALL = [], [], []
    FREQMIN, FREQMAX = [], []
    IRFPOW = []
    PRESOUT = []
    IMOVIE_FRAMES = []
    header_icpdat = []
    AtomicSpecies, NeutSpecies, PosSpecies, NegSpecies = [], [], [], []

    # Loop over all folders and retrieve mesh sizes and SI sizes.
    for index in range(numfolders):
        # ==========##===== INITMESH.OUT READIN =====##==========#
        # ==========##===============================##==========#

        # Attempt automated retrieval of mesh sizes.
        try:
            # Identify mesh size from TECPLOT2D file. (Data Array Always Correct Size)
            meshdata = open(TEC2D[index]).readlines()

            # Zone line holds data, split at comma, R&Z values are given by "I=,J=" respectively.
            # R_py3_object = list(filter(lambda x: x.isdigit(), R))
            # Z_py3_object = list(filter(lambda x: x.isdigit(), Z))
            R = list(
                filter(lambda x: "ZONE" in x, meshdata)
            )  # String: 'ZONE I=xxx, J=xxx, F=BLOCK"
            Z = list(
                filter(lambda x: "ZONE" in x, meshdata)
            )  # String: 'ZONE I=xxx, J=xxx, F=BLOCK"
            R = (
                R[0].split(",")[0].strip(" \t\n\r,=ZONE I")
            )  # Split at commas, [0] gives "I=xxx"
            Z = (
                Z[0].split(",")[1].strip(" \t\n\r,=ZONE J")
            )  # Split at commas, [1] gives "J=xxx"
            R_mesh.append(int(R))  # R_mesh (Cells) [int]
            Z_mesh.append(int(Z))  # Z_mesh (Cells) [int]

        # If extraction from TECPLOT2D file fails, attempt to extract from initmesh.out header
        # This is an old method and causes issues with Q-VT meshes and magnified meshes
        except ValueError:
            # Identify mesh size from initmesh.out header:
            meshdata = open(mesh[index]).readline()
            R_mesh.append([int(i) for i in meshdata.split()][1])
            if Magmesh == 1:
                Z_mesh.append([int(i) + 1 for i in meshdata.split()][3])
            elif Magmesh == 2:
                Z_mesh.append([int(i) + 3 for i in meshdata.split()][3])
            elif Magmesh == 3:
                Z_mesh.append([int(i) + 5 for i in meshdata.split()][3])
        # endif

        # If all else fails, request manual input of mesh resolution
        except:
            # If data for current file exists:
            if index <= len(TEC2D) - 1:
                # If the initmesh.out file cannot be found, manual input is required.
                print("ERR: ICP.NAM GEOMETRY READIN, USING MANUAL MESH CELL INPUT:")
                print(Dirlist[index])
                r_mesh = int(input("DEFINE NUM RADIAL CELLS:"))
                z_mesh = int(input("DEFINE NUM AXIAL CELLS:"))
                print("")

                R_mesh.append(r_mesh)
                Z_mesh.append(z_mesh)
        # endif
        # endtry

        # ==========##===== ICP.NAM READIN =====##==========#
        # ==========##==========================##==========#

        # Attempt automated retrieval of SI conversion units.
        NamelistData = open(icpnam[index]).readlines()

        # Mesh Geometry Namelist Inputs
        try:
            RADIUS = list(filter(lambda x: "RADIUS=" in x, NamelistData))
            RADIUS = RADIUS[0].split("!!!")[0]
            RADIUS = float(RADIUS.strip(" \t\n\r,=RADIUS"))

            RADIUST = list(filter(lambda x: "RADIUST=" in x, NamelistData))
            RADIUST = RADIUST[0].split("!!!")[0]
            RADIUST = float(RADIUST.strip(" \t\n\r,=RADIUST"))

            HEIGHT = list(filter(lambda x: "HEIGHT=" in x, NamelistData))
            HEIGHT = HEIGHT[0].split("!!!")[0]
            HEIGHT = float(HEIGHT.strip(" \t\n\r,=HEIGHT"))

            HEIGHTT = list(filter(lambda x: "HEIGHTT=" in x, NamelistData))
            HEIGHTT = HEIGHTT[0].split("!!!")[0]
            HEIGHTT = float(HEIGHTT.strip(" \t\n\r,=HEIGHTT"))

            DEPTH = list(filter(lambda x: "DEPTH=" in x, NamelistData))
            DEPTH = DEPTH[0].split("!!!")[0]
            DEPTH = float(DEPTH.strip(" \t\n\r,=DEPTH"))

            IXZ = list(filter(lambda x: "IXZ=" in x, NamelistData))
            IXZ = IXZ[0].split("!!!")[0]
            IXZ = int(IXZ.strip(" \t\n\r,=IXZ"))

            ISYM = list(filter(lambda x: "ISYM=" in x, NamelistData))
            ISYM = ISYM[0].split("!!!")[0]
            ISYM = int(ISYM.strip(" \t\n\r,=ISYM"))

            # ISYMlist[l] = 1 if mesh uses radial symmetry, = 0 if not
            if image_plotsymmetry:
                ISYMlist = np.append(ISYMlist, ISYM)
            else:
                ISYMlist.append(0)

            # IXZlist[l] = 1 if mesh uses cartesian coordinates, = 0 if cylindrical
            IXZlist = np.append(IXZlist, IXZ)

            # Determine if mesh RADIUS or RADIUST was used, save 'Radius' used for further calculations
            if RADIUS > 0.0:
                Radius = np.append(Radius, RADIUS)  # [cm]
            elif RADIUST > 0.0:
                Radius = np.append(Radius, RADIUST)  # [cm]
            if HEIGHT > 0.0:
                Height = np.append(Height, HEIGHT)  # [cm]
            elif HEIGHTT > 0.0:
                Height = np.append(Height, HEIGHTT)  # [cm]

            # Determine mesh cell radial (dr), axial (dz), and depth (dy) resolutions
            dr = np.append(dr, Radius[-1] / (R_mesh[-1] - 1))  # [cm/cell]
            dz = np.append(dz, Height[-1] / (Z_mesh[-1] - 1))  # [cm/cell]
            dy = np.append(dy, DEPTH)  # [cm/cell] (DEPTH is always 1 cell)

        except:
            # If the geometry section cannot be found, manual input is required.
            print("ERR: ICP.NAM GEOMETRY READIN, USING MANUAL MESH SI INPUT:")
            print(Dirlist[index])
            radius = float(input("DEFINE RADIUST [cm]:"))
            height = float(input("DEFINE HEIGHTT [cm]:"))
            depth = float(input("DEFINE DEPTH [cm]:"))
            print("")

            Radius.append(radius)
            Height.append(height)
            Depth.append(depth)
            dr.append(Radius[-1] / (R_mesh[-1] - 1))
            dz.append(Height[-1] / (Z_mesh[-1] - 1))
        # endtry

        # =====#=====#

        # Namelist Timescale inputs

        DTPOS = list(filter(lambda x: "DTPOS=" in x, NamelistData))
        DTPOS = DTPOS[0].split("!!!")[0]
        DTPOS = float(DTPOS.strip(" \t\n\r,=DTPOS"))

        # =====#=====#

        # Material Namelist Inputs (frequencies/voltages/powers)
        # NUMMETALS,NUMCOILS					::	Integer, [-]
        # CMETALS,CCOILS						::	1D list, [String]
        # FREQM,FREQM2,FREQC,FREQGLOB 		::	1D list, [Hz]
        # VRFM,VRFM2 							::	1D list, [V]
        # IRFPOW 								::	1D list, [W]
        # PRESOUT 							::	1D list, [Torr]
        try:
            NUMMETALS = list(filter(lambda x: "IMETALS" in x, NamelistData))[0].strip(
                " \t\n\r,=IMETALS"
            )
            NUMMETALS = int(NUMMETALS.split("!!!")[0].strip(" \t\n\r,=IMETALS"))
            CMETALS = (
                list(filter(lambda x: "CMETAL=" in x, NamelistData))[0]
                .strip(" \t\n\r,")
                .split()[1 : NUMMETALS + 1]
            )
            for i in range(len(CMETALS)):
                CMETALS[i] = str(CMETALS[i].strip(","))
            IETRODEM = list(filter(lambda x: "IETRODEM=" in x, NamelistData))[
                0
            ].split()[1 : NUMMETALS + 1]
            for i in range(len(IETRODEM)):
                IETRODEM[i] = int(IETRODEM[i].strip(","))
        except:
            print("ERROR: FAILED TO READ IMETALS, CMETAL, IETRODEM")
            print('SEE "#Material Namelist Inputs (frequencies/voltages/powers)"')
            exit()
        ##
        try:
            NUMCOILS = list(filter(lambda x: "ICOILS" in x, NamelistData))
            for i in range(len(NUMCOILS)):
                if "ICOILSOL" not in NUMCOILS[i]:
                    NUMCOILS = NUMCOILS[i]
                    break
            # endif
            # endfor
            NUMCOILS = int(NUMCOILS.split("!!!")[0].strip(" \t\n\r,=ICOILS"))
            CCOILS = (
                list(filter(lambda x: "CCOIL=" in x, NamelistData))[0]
                .strip(" \t\n\r,")
                .split()[1:NUMCOILS]
            )
            for i in range(len(CCOILS)):
                CCOILS[i] = str(CCOILS[i].strip(","))
        except:
            print("ERROR: FAILED TO READ ICOILS, CCOIL")
            print('SEE "#Material Namelist Inputs (frequencies/voltages/powers)"')
            exit()
        ##
        try:
            VRFM.append(
                list(filter(lambda x: "VRFM=" in x, NamelistData))[0].split()[
                    1:NUMMETALS
                ]
            )
            for i in range(len(VRFM[-1])):
                VRFM[-1][i] = VRFM[-1][i].strip(",")
                VRFM[-1][i] = float(VRFM[-1][i].replace("D", "E"))
            # endfor
            VRFM2.append(
                list(filter(lambda x: "VRFM_2=" in x, NamelistData))[0].split()[
                    1:NUMMETALS
                ]
            )
            for i in range(len(VRFM2[-1])):
                VRFM2[-1][i] = VRFM2[-1][i].strip(",")
                VRFM2[-1][i] = float(VRFM2[-1][i].replace("D", "E"))
        # endfor
        except:
            print("ERROR: FAILED TO READ VRFM, VRFM_2")
            print('SEE "#Material Namelist Inputs (frequencies/voltages/powers)"')
            exit()
        ##
        try:
            FREQM.append(
                list(filter(lambda x: "FREQM=" in x, NamelistData))[0].split()[
                    1:NUMMETALS
                ]
            )
            for i in range(len(FREQM[-1])):
                FREQM[-1][i] = FREQM[-1][i].strip(",")
                FREQM[-1][i] = float(FREQM[-1][i].replace("D", "E"))
            # endfor
            FREQM2.append(
                list(filter(lambda x: "FREQM_2=" in x, NamelistData))[0].split()[
                    1:NUMMETALS
                ]
            )
            for i in range(len(FREQM2[-1])):
                FREQM2[-1][i] = FREQM2[-1][i].strip(",")
                FREQM2[-1][i] = float(FREQM2[-1][i].replace("D", "E"))
        # endfor
        except:
            print("ERROR: FAILED TO READ FREQM, FREQM_2")
            print('SEE "#Material Namelist Inputs (frequencies/voltages/powers)"')
            exit()
        ##
        try:
            FREQC.append(
                list(filter(lambda x: "FREQC=" in x, NamelistData))[0].split(",")[0:-1]
            )
            for i in range(len(FREQC[index])):
                FREQC[-1][i] = FREQC[index][i].strip(" \t\n\r,=FREQC")
                FREQC[-1][i] = float(FREQC[index][i].replace("D", "E"))
        # endfor
        except:
            print("ERROR: FAILED TO READ FREQC")
            print('SEE "#Material Namelist Inputs (frequencies/voltages/powers)"')
            exit()
        try:
            ##
            FREQGLOB.append(
                list(filter(lambda x: "FREQ=" in x, NamelistData))[0].split()
            )
            FREQGLOB[index] = FREQGLOB[index][0].strip(
                " \t\n\r,=FREQ"
            )  # NOTE: Assumes 1 entry for "FREQ="
            FREQGLOB[index] = float(
                FREQGLOB[index][0].replace("D", "E")
            )  # NOTE: Assumes 1 entry for "FREQ="
            ##
            IRFPOW.append(
                list(filter(lambda x: "IRFPOW=" in x, NamelistData))[0].strip(
                    " \t\n\r,=IRFPOW"
                )
            )
            IRFPOW[index] = float(IRFPOW[index].split("!!!")[0].strip(" \t\n\r,"))
            ##
            PRESOUT.append(
                list(filter(lambda x: "PRESOUT=" in x, NamelistData))[0].strip(
                    " \t\n\r,=PRESOUT"
                )
            )
            PRESOUT[index] = float(PRESOUT[index].split("!!!")[0].strip(" \t\n\r,"))
        except:
            print("ERROR: FAILED TO CONVERT FREQ")
            print('SEE "#Material Namelist Inputs (frequencies/voltages/powers)"')
            exit()
        # endtry

        # No ICP coils are on - Ignore ICP frequencies in FREQALL
        if int(NUMCOILS) == 0:
            FREQALL = FREQM[index] + FREQM2[index] + FREQGLOB
            FREQALL = [x for x in FREQALL if x > 0]  # Ignore any DC voltages (FREQ = 0)
        # ICP coils are on - include ICP frequencies in FREQALL
        elif int(NUMCOILS) > 0:
            FREQALL = FREQM[index] + FREQM2[index] + FREQC[index] + FREQGLOB
            FREQALL = [x for x in FREQALL if x > 0]  # Ignore any DC voltages (FREQ = 0)
        # endif
        FREQMIN.append(min(FREQALL))  # Minimum global frequency 	[Hz]
        FREQMAX.append(max(FREQALL))  # Maximum global frequency	[Hz]

        # Debug to check all RF frequencies
        if IDEBUG == 1:
            print("Material Frequencies FREQM", FREQM[index])
            print("Material Frequencies FREQM2", FREQM2[index])
            print("Coil Freqencies FREQC", FREQC[index])
            print("Global Frequency FREQ", FREQGLOB)
        # endif

        # =====#=====#

        # Plasma Chemistry Monte-Carlo (PCMC) Namelist Inputs
        try:
            IMATSTATS = list(filter(lambda x: "IMATSTATS" in x, NamelistData))[0].strip(
                " \t\n\r,=IMATSTATS"
            )
            IMATSTATS = int(IMATSTATS[0].split("!!!")[0])
            CMATSTATS = list(filter(lambda x: "CMATSTATS=" in x, NamelistData))[
                0
            ].strip(" \t\n\r,")
            CMATSTATS = CMATSTATS.split("=")[1].split(",")[0:IMATSTATS]
            ##
            IPCMCSPEC = list(filter(lambda x: "IPCMCSPEC" in x, NamelistData))[0].strip(
                " \t\n\r,=IPCMCSPEC"
            )
            IPCMCSPEC = (
                int(IPCMCSPEC[0].split("!!!")[0]) + 1
            )  # Add 1 to account for default ION-TOT
            CPCMCSPEC = list(filter(lambda x: "CPCMCSPEC=" in x, NamelistData))[
                0
            ].strip(" \t\n\r,")
            CPCMCSPEC = CPCMCSPEC.split("=")[1].split(",")[0:IPCMCSPEC]
            ##
            IEBINSPCMC = list(filter(lambda x: "IEBINSPCMC=" in x, NamelistData))
            IEBINSPCMC = IEBINSPCMC[0].split("!!!")[0]
            IEBINSPCMC = float(IEBINSPCMC.split()[0].strip(" \t\n\r,=IEBINSPCMC"))
            EMAXIPCMC = list(filter(lambda x: "EMAXIPCMC=" in x, NamelistData))
            EMAXIPCMC = EMAXIPCMC[0].split("!!!")[0]
            EMAXIPCMC = float(EMAXIPCMC.split()[0].strip(" \t\n\r,=EMAXIPCMC"))
        except:
            print("ERROR: ICP.NAM PCMC READIN, USING DEFAULT PCMC PROPERTIES")
            IEBINSPCMC = 1000
            EMAXIPCMC = 50
        # endtry

        # Phase-Resolved IMOVIE Namelist Inputs
        try:
            imovie_frames = list(
                filter(lambda x: "IMOVIE_FRAMES=" in x, NamelistData)
            )  # OLD VERSION HAD [0]
            imovie_frames = imovie_frames[0].split("!!!")[0]
            IMOVIE_FRAMES.append(int(imovie_frames.strip(" \t\n\r,=IMOVIE_FRAMES")))
        except:
            print(
                "ERROR: ICP.NAM IMOVIE READIN, USING DEFAULT PHASE RESOLVED PROPERTIES"
            )
            IMOVIE_FRAMES.append(180)
        # endtry

        # ==========##===== ICP.DAT READIN =====##==========#
        # ==========##==========================##==========#

        # Automated retrieval of atomic species
        try:
            ChemistryData = open(icpdat[index], encoding="utf-8").readlines()
        except:
            ChemistryData = open(icpdat[index], encoding="iso-8859-15").readlines()
        # endtry

        # Determine reaction mechanism (.dat) file header length
        for i in range(len(ChemistryData)):
            Line = ChemistryData[i].split()
            if Line[0] == "*":
                DatHeaderLen = i
                break
        # endif
        # endfor

        # Extract reaction mechanism (.dat) header, atomic species definitions
        try:
            for i in range(DatHeaderLen):
                Line = ChemistryData[i]
                SplitLine = ChemistryData[i].split()

                # Ignore lines starting with '>' as these are surface return fractions
                # Ignore lines where the first entry is '!' as this is a comment
                # Ignore lines where the first entry contains '!' as this is a comment
                # Minimum split length for atomic entry is 13 when split
                if (
                    SplitLine[0] not in "> !"
                    and "!" not in SplitLine[0]
                    and len(SplitLine) >= 13
                ):
                    # Split line in order of ":, ;, &, ], [, @, and !", and strip these characters
                    StripLine = re.split(r"\s*[:;&\]\[@!]\s*", Line.strip())
                    header_icpdat.append(StripLine)

                    SpeciesName = str(StripLine[0])  # -
                    Charge = int(StripLine[1])  # [e]
                    MolecularWeight = float(StripLine[2])  # [amu]
                    StickingCoeff = float(StripLine[3])  # -
                    TransportBool = int(StripLine[4])  # -
                    ReturnFrac = float(StripLine[5])  # -
                    ReturnSpecies = str(StripLine[6])  # -
                    Comments = str(StripLine[7::])  # -

                    # Identify atomic species (including electrons)
                    if SpeciesName not in AtomicSpecies:
                        AtomicSpecies.append(SpeciesName)

                    # Identify ground-state neutral species
                    # NOTE :: NEED TO DEVELOP ALGORITHM FOR GROUND STATE
                    FluidSpecies = ["AR", "AR3S", "O2", "O"]

                    # Identify positive and negative species (including electrons)
                    if Charge == 0 and SpeciesName not in NeutSpecies:
                        NeutSpecies.append(SpeciesName)
                    elif Charge >= 1 and SpeciesName not in PosSpecies:
                        PosSpecies.append(SpeciesName)
                    elif Charge <= -1 and SpeciesName not in NegSpecies:
                        NegSpecies.append(SpeciesName)
            # endif
        # endfor

        except:
            print(
                "ERROR: CANNOT READ .DAT HEADER ATOMIC DEFINITIONS, CHECK ICP.DAT HEADER"
            )
    # endtry

    return (
        R_mesh,
        Z_mesh,
        ISYMlist,
        IXZlist,
        Radius,
        Height,
        Depth,
        dr,
        dz,
        dy,
        VRFM,
        VRFM2,
        FREQM,
        FREQM2,
        FREQC,
        FREQGLOB,
        FREQALL,
        FREQMIN,
        FREQMAX,
        IRFPOW,
        PRESOUT,
        IMOVIE_FRAMES,
        header_icpdat,
        AtomicSpecies,
        FluidSpecies,
        NeutSpecies,
        PosSpecies,
        NegSpecies,
    )
