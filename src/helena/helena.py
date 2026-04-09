#!/usr/bin/env python3

#################################
#		Point of Contact		#
#								#
#	    Dr. Scott J. Doyle      #
#	  University of Aberdeen    #
#	  Department of Physics     #
#     Kings College, Aberdeen   #
#           AB21 0HD, UK        #
#	  Scott.Doyle@Physics.org   #
#                               #
#################################
#           'HELENA'            #
# Hpem ELectric ENgine Analysis #
#################################


#====================================================================#
				 #PROGRAM FLAGS AND MODULE IMPORTS#
#====================================================================#

from pylab import *
from .io import get_directories, write_to_csv, write_data_to_file, read_data_from_file, create_new_folder, \
	folder_name_trimmer, make_movie
from .initialisation import get_mesh_and_SI
from .data import extract_raw_data, read_TEC2D, read_TEC2D_phase, read_geometry
from .variables import enumerate_variables, enumerate_vectors, variable_interpolator, variable_unit_conversion, \
	azimuthal_phase_conversion, variable_label_maker
from .utility import string_in_variable, is_radial_variable
from .external import auto_conv_prof

def run(argv=None):
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-f", "--first", action="store_true", dest="install", default=False, help="Install prompt for required python modules")
	import sys
	if argv is None:
		argv = sys.argv[1:]
	(options, args) = parser.parse_args(argv)
	input_file = "input.toml" if len(args) == 0 else args[0]

	if 'True' in str(options):
		import os, sys
		import os.path

		print( '')
		print( 'First time use requires installation of additional python modules')
		print( 'Please type your password when prompted to allow installation:')
		print( '')
		try:
			os.system('sudo apt-get install python-pip')
			os.system('sudo apt-get install python-matplotlib')
			os.system('sudo apt-get install python-numpy')
			os.system('sudo apt-get install python-scipy')
			os.system('sudo apt-get install ffmpeg')
			os.system('pip install tqdm')
		except:
			print( '')
			print( 'Error installing required packages')
			print( 'Please attempt manual installation')
			print( '')
		#endtry
		print( '')
		print( '')
	#endif

	#==============#

	#Import Core Modules
	import matplotlib.cm as cm
	import numpy as np
	import scipy as sp
	import math as m
	import subprocess
	import warnings
	import os, sys
	import os.path
	import csv
	import re
	import gc
	import toml
	from typing import TypeVar, Callable

	#Enforce matplotlib to avoid instancing undisplayed windows
	#matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
	import matplotlib as mpl
	#mpl.use('Agg')			# Uncomment to fix mem leak, disables GUI backend and "plot()"

	#Import Additional Modules
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from scipy.signal import savgol_filter
	from matplotlib import pyplot as plt
	from matplotlib import colors as col
	from matplotlib import ticker
	from scipy import ndimage
	from tqdm import tqdm
	from pydantic import Field
	from pydantic.dataclasses import dataclass
	#from Read_data_functions.py import overlay_GEC_geometry


    #Commonly used variable sets.
	variable_sets = {
		"Phys": ['E','S-E','SEB-E','TE','PPOT','P-POT','POW-RF','POW-RF-E','POW-ICP','POW-ICP1','POW-ICP2','POW-ICP3',
				 'POW-ICP4','POW-ALL','EB-ESORC','COLF','SIGMA','EF-TOT', 'ERADIAL','ETHETA','EAXIAL','PHASEER','PHASE',
				 'PHASEEZ','EAMB-Z','EAMB-R','RHO','BR','BRS','BZ','BZS','BT','BTS','BRF', 'PHASEBR','PHASEBT',
				 'PHASEBZ','VR-ION+','VZ-ION+','E FLUX-R','FR-E','E FLUX-Z','FZ-E','JZ-NET','JR-NET','J-THETA',
				 'J-TH(MAG)','J-TH(PHA)','PRESSURE','TG-AVE','VR-NEUTRAL','VZ-NEUTRAL'],
		"PhysCoilsEF": ['ERADIAL-2','ETHETA-2','EAXIAL-2','PHASEER-2','PHASEEZ-2','ERADIAL-3','ETHETA-3','EAXIAL-3',
						'PHASEER-3','PHASEEZ-3','ERADIAL-4','ETHETA-4','EAXIAL-4','PHASEER-4','PHASEEZ-4','ERADIAL-5',
						'ETHETA-5','EAXIAL-5','PHASEER-5','PHASEEZ-5','ERADIAL-6','ETHETA-6','EAXIAL-6','PHASEER-6',
						'PHASEEZ-6','ERADIAL-7','ETHETA-7','EAXIAL-7','PHASEER-7','PHASEEZ-7','ERADIAL-8','ETHETA-8',
						'EAXIAL-8','PHASEER-8','PHASEEZ-8'],
		"PhysCoilsBF": ['BT-2','BT-3','BT-4','BT-5','BT-6','BT-7','BT-8','BRF-2','BRF-3','BRF-4','BRF-5','BRF-6',
						'BRF-7','BRF-8','PHASEBT-2','PHASEBT-3','PHASEBT-4','PHASEBT-5','PHASEBT-6','PHASEBT-7',
						'PHASEBT-8'],

		"Conv": ['E','TE','PPOT','POW-RF','SIGMA','EF-TOT','TG-AVE'],

		"TEST": ['AR2+'],
		"Ar": ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR','AR+','AR2+','AR2*','S-AR+','S-AR4P',
			   'SEB-AR+','SEB-AR4P','FZ-AR3S','FR-AR3S','FR-AR+','FZ-AR+','FZ-AR3S','FR-AR3S'],
		"O2": ['O3','O2','O2V','O2*','O2*1S','O2+','O2-','O','O1S','O+','O-','O*','S-O3','S-O2+','S-O+','S-O-','SEB-O3',
			   'SEB-O+','SEB-O2+','SEB-O-','FR-O+','FZ-O+','FR-O-','FZ-O-'],
		"H2": ['H2V0','H2V1','H2V2','H2V3','H1','H*','H**','H2+','H+','H-','S-H+','SEB-H+','S-2H+','SEB-2H+','S-H-',
			   'SEB-H-','FZ-H2V0','FR-H2V0','FZ-H1','FR-H1','FZ-H+','FR-H+','FZ-H2+','FR-H2+','FZ-H-','FR-H-'],
		"N2": ['N2','N2V','N2*','N2**','N2+','N','N*','N+'],
		"Cl": ['Cl2','Cl','CL+','CL-','Cl2V','Cl2+','CL*','CL**','CL***'],
		"F": ['F2','F2*','F2+','F','F*','F+','F-','S-F','S-F+','S-F-','SEB-F','SEB-F+','SEB-F-','FZ-F','FR-F','FZ-F+',
			  'FR-F+','FZ-F-','FR-F-','FZ-F+','FR-F+'],
		"H2O": ['H2O','H2O+','OH','OH-','H2OV','H2O2','S-H2O','SEB-H2O','S-H2OV','SEB-H2OV','S-H2O+','SEB-H2O+','S-OH',
				'SEB-OH','S-OH-','SEB-OH-','S-OH+','SEB-OH+'],
		"COx": ['CO2','CO2V','CO+','CO','CO+','C','C+'],
		"CHx": ['CH4','CH3','CH2','CH','C','CH5+','CH4+','CH3+','CH2+','CH+','C+'],
		"NHx": ['NH3','NH2','NH','NH4+','NH3+','NH2+','NH+'],
		"NOx": ['NO2','NO2+','N2O','N2O+','NO','NO+'],
		"NFx": ['NF3A','NF2A','NFA','NF3B','NF2B','NFB','NF3+','NF2+','NF+'],
		"SFx": ['SF6','SF5','SF4','SF3','SF2','SF','S','SF5+','SF4+','SF3+','SF2+','SF+','S+','SF6-','SF5-'],
		"Al": ['AL','AL*','AL**','AL+','S-AL','SEB-AL','S-AL*','SEB-AL*','S-AL**','SEB-AL**','S-AL+','SEB-AL+','FZ-AL+',
			   'FR-AL+'],
		"Be": ['BE','BE1','BE2','BE3','BE4','BE5','BE6','BE7','BE8','BE9','BE+','FR-BE','FZ-BE','FR-BE+','FZ-BE+'],

		"Ar_Phase": ['S-E','S-AR+','S-AR4P','SEB-E','SEB-AR+','SEB-AR4P','SRCE-2437','FR-E','FZ-E','TE','PPOT',
					 'POW-ALL'],
		"O2_Phase": ['S-E','S-O+','S-O-','S-O2+','SEB-O+','SEB-O-','SEB-O2+','S-O3P3P','SEB-O3P3P','TE','PPOT','FR-E',
					 'FZ-E'],

		"PRCCPAr_PCMC": ['AR^0.35','EB-0.35','ION-TOT0.35'],
		"PRCCPO2_PCMC": ['O^0.35','EB-0.35','ION-TOT0.35'],
		"GECCCP2a_BE_PCMC": ['AR^    2.6 L','BE^    2.6 L','EB-    2.6 L'],
	}


	_T = TypeVar("T")


	@dataclass
	class ChemistryInput:
		variables: list[str] = Field(default_factory=lambda: variable_sets["Phys"] + variable_sets["Ar"])
		multivar: list[str] = Field(default_factory=list)
		radial_profiles: list[int] = Field(default_factory=list)
		axial_profiles: list[int] = Field(default_factory=lambda: [0])
		probe_loc: list[int] = Field(default_factory=list)

		# EDF.
		IEDF_variables: list[str] = Field(default_factory=lambda: variable_sets["GECCCP2a_BE_PCMC"].copy())
		NEDF_variables: list[str] = Field(default_factory=list)

		# Movie.
		phase_variables: list[str] = Field(default_factory=lambda: variable_sets["Ar_Phase"].copy())
		electrode_loc: list[int] = Field(default_factory=lambda: [0, 0])
		waveform_locs: list[int] = Field(default_factory=list)

		# Diagnostic.
		sheath_ROI: list[int] = Field(default_factory=list)
		source_width: list[int] = Field(default_factory=list)

		def __post_init__(self):
			self.variables = self.expand_variable_set(self.variables)
			self.multivar = self.expand_variable_set(self.multivar)
			self.IEDF_variables = self.expand_variable_set(self.IEDF_variables)
			self.NEDF_variables = self.expand_variable_set(self.NEDF_variables)
			self.phase_variables = self.expand_variable_set(self.phase_variables)

		@staticmethod
		def expand_variable_set(items: list[str]) -> list[str]:
			out = []

			for item in items:
				if item in variable_sets:
					out.extend(variable_sets[item])
				else:
					out.append(item)

			return out

	@dataclass
	class FiguresInput:
		tecplot2D: bool = False

		# Movie.
		movieicp2D: bool = False
		movieicp1D: bool = False
		timeaxis1D: bool = False
		convergence: bool = False
		iterstep: int = 1

		# Profiles.
		monoprofiles: bool = False
		multiprofiles: bool = False
		compare_profiles: bool = False

		# Trends.
		trend_phase_averaged: bool = False
		trend_phase_resolved: bool = False
		thrust_loc: int = 45

		# Resolve.
		phase_resolve2D: bool = False
		phase_resolve1D: bool = False
		sheath_dynamics: bool = False
		PROES: bool = False
		phase_cycles: float = 1.01
		do_Fwidth: int = 0

		# EDF.
		IEDF_angular: bool = False
		IEDF_trends: bool = False
		EEDF: bool = False

	@dataclass
	class WriteInput:
		ASCII: bool = False
		CSV: bool = True

	@dataclass
	class PrintoutInput:
		general_trends: bool = False
		Knudsen_number: bool = False
		total_power: bool = False
		Reynolds: bool = False
		DC_bias: bool = False
		thrust: bool = False
		sheath: bool = False

	@dataclass
	class ImageInput:
		extension: str = ".png"

		# Style.
		interp: str = "spline36"
		cmap: str = "plasma"

		# Geometry.
		aspect_ratio: list[float] = Field(default_factory=lambda: [10.0, 10.0])
		radial_crop: list[float] = Field(default_factory=list)
		axial_crop: list[float] = Field(default_factory=list)
		rotate: bool = False

		# Display.
		plot_symmetry: bool = False
		plot_mesh: bool = True
		plot_grid: bool = False

		# Contours.
		plot_colourfill: bool = True
		plot_contours: bool = True
		contour_lvls: int = 10

		# Axes.
		axis_ticks: bool = True
		axis_labels: bool = True
		legend_loc: str = 'best'

		# Colorbar.
		cbar_ticks: bool = True
		cbar_bins: int = 5
		cbar_limit: list[float] = Field(default_factory=list)

		# Vector.
		plot_vector: bool = True
		vector_density: float = 1.5
		vector_lw: float = 1.0

		# Processing.
		normalise: bool = False
		log_plot: bool = False

		# Overlays.
		plot_sheath: bool = False
		plot_phase_waveform: bool = False
		plot_overlay: bool = False


	@dataclass
	class OverridesInput:
		title: list[str] = Field(default_factory=list)
		legend: list[str] = Field(default_factory=list)
		xaxis: list[str] = Field(default_factory=list)
		xlabel: list[str] = Field(default_factory=list)
		ylabel: list[str] = Field(default_factory=list)

	@dataclass
	class ExpertInput:
		# Debug.
		magmesh: int = 1
		ffmpeg_movies: bool = False
		IDEBUG: bool = False

		# Warnings.
		ignore_div_zero: bool = True
		ignore_empty_contours: bool = True

		# Numerical.
		sheath_method: str = 'AbsDensity'
		thrust_method: str = 'AxialMomentum'
		mean_calculation: str = 'MeanFraction'

		# Overrides.
		DC_bias_axis: str = 'Auto'
		sheath_ion_species: list[str] = Field(default_factory=lambda: ['AR+'])

		# Filtering.
		kinetic_filtering: bool = True
		plot_kinetic_filtering: bool = False
		sav_window: int = 25
		sav_poly_order: int = 3

		# Misc.
		sheath_ion_ratio_threshold: float = 1.03
		conv_azimuthal_phase: bool = True
		EDF_threshold: float = 0.01
		units: str = 'SI'


	@dataclass
	class Config:
		chemistry: ChemistryInput = Field(default_factory=ChemistryInput)
		figures: FiguresInput = Field(default_factory=FiguresInput)
		write: WriteInput = Field(default_factory=WriteInput)
		printout: PrintoutInput = Field(default_factory=PrintoutInput)
		image: ImageInput = Field(default_factory=ImageInput)
		overrides: OverridesInput = Field(default_factory=OverridesInput)
		expert: ExpertInput = Field(default_factory=ExpertInput)


	def create_input(data: dict, category: str, output_obj: type[_T], sections: list[str] | None = None) -> _T:
		# data: the whole toml data
		# category: the specific category we are creating an input for
		# output_obj: the category's output object
		# sections: the additional sections of this category

		if sections is None:
			sections = []

		category_dict = data.get(category, {})

		output_dict = {k: v for k, v in category_dict.items() if k not in sections}

		for section in sections:
			output_dict.update(category_dict.get(section, {}))

		return output_obj(**output_dict)


	def load_config(toml_dict):
		chemistry_input = create_input(toml_dict, "chemistry", ChemistryInput, ["EDF", "movie", "diagnostic"])
		figures_input = create_input(toml_dict, "figures", FiguresInput, ["movie", "profiles", "trends", "resolve", "EDF"])
		write_input = create_input(toml_dict, "write", WriteInput)
		printout_input = create_input(toml_dict, "print", PrintoutInput)
		image_input = create_input(toml_dict, "image", ImageInput, ["style", "geometry", "display", "contours", "axes", "colorbar", "vector", "processing", "overlays"])
		overrides_input = create_input(toml_dict, "overrides", OverridesInput)
		expert_input = create_input(toml_dict, "expert", ExpertInput, ["debug", "warnings", "numerical", "overrides", "filtering", "misc"])

		return Config(
			chemistry=chemistry_input,
			figures=figures_input,
			write=write_input,
			printout=printout_input,
			image=image_input,
			overrides=overrides_input,
			expert=expert_input,
		)

	with open(input_file, "r") as toml_file:
		toml_dict = toml.load(toml_file)

	config = load_config(toml_dict)

	#====================================================================#
							 #LOW LEVEL INPUTS#
	#====================================================================#

	#Various debug and streamlining options.
	Magmesh = config.expert.magmesh							#initmesh.exe magnification factor. (Obsolete - legacy)
	ffmpegMovies = config.expert.ffmpeg_movies				#If False: Suppresses ffmpeg routines, saves RAM.
	IDEBUG = config.expert.IDEBUG							#Produces debug outputs for most diagnostics.

	#Warning suppressions
	if config.expert.ignore_div_zero:
		np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
		#Fix: "can't invoke "event" command: application has been destroyed" error with PROES images
		#Fix: "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

	if config.expert.ignore_empty_contours:
		warnings.filterwarnings("ignore", message="No contour levels were found within the data range.")
		#Fix: Suppress above warning on plotting of empty contour plots

	#Numerical Calculation Methods:
	GlobSheathMethod = config.expert.sheath_method			#Set Global Sheath Calculation Method.
	GlobThrustMethod = config.expert.thrust_method			#Set Global Thrust Calculation Method.
	GlobMeanCalculation = config.expert.mean_calculation	#Definition of 'mean' EDF value

	#Overrides or 'fudge factors' for diagnostics
	DCbiasaxis = config.expert.DC_bias_axis					#Force Direction Over Which DCBias is Calculated
	SheathIonSpecies = config.expert.sheath_ion_species		#Force Sheath Ion Species (blank for auto)

	#Ratio of electrons to ions that determines the edge of the sheath
	#Ideally this should be 1.00, but in practice it's slightly over 1 at lower resolutions
	Sheath_IonRatio_Threshold = config.expert.sheath_ion_ratio_threshold

	#Data Filtering and Smoothing Methods:
	KineticFiltering = config.expert.kinetic_filtering				#Pre-fit kinetic data employing a SavGol filter
	PlotKineticFiltering = config.expert.plot_kinetic_filtering		#Plot Filtered Profiles, or employ only in trends.
	Glob_SavWindow = config.expert.sav_window						#Window > FeatureSize
	Glob_SavPolyOrder = config.expert.sav_poly_order				#Polyorder ~= Smoothness

	#Apply azimuthal direction (phase) to relevant variables if true, else plot magnitude only
	ConvAzimuthalPhase = config.expert.conv_azimuthal_phase

	#Minimum plotted EDF energy fraction, cuts x-axis at index where :: f(e) = f(e)*EDF_Threshold
	#Note: IEDF/EEDF trends are only taken within range :: EDF_threshold < f(e) < 1.0
	EDF_Threshold = config.expert.EDF_threshold				# i.e. = 0.0 to plot all

	#Define units for particular variables
	Units = config.expert.units					#'SI','CGS'
												# FUNCTION MATHS ASSUMES SI, CGS WILL GIVE INCORRECT RESULTS

	#====================================================================#
						#SWITCHBOARD AND DIAGNOSTICS#
	#====================================================================#

	# Requested IEDF/NEDF Variables.
	IEDFVariables = config.chemistry.IEDF_variables			# Requested Variables from iprofile_2d.pdt
	NEDFVariables = config.chemistry.NEDF_variables			# Requested Variables from nprofile_2d.pdt

	# Requested movie1/movie_icp Variables.
	PhaseVariables = config.chemistry.phase_variables		# Requested Movie1 (phase) Variables.
	electrodeloc = config.chemistry.electrode_loc			# Cell location of powered electrode [R,Z].
	waveformlocs = config.chemistry.waveform_locs			# Cell locations of additional waveforms [R,Z].

	# Requested variables and plotting locations.
	Variables = config.chemistry.variables					# Requested Variables from Tecplot2D.pdt, tecplot_kin.pdt, and movie_icp.pdt
	multivar = config.chemistry.multivar					# Additional variables plotted ontop of [Variables]
	radialprofiles = config.chemistry.radial_profiles		# Radial 1D-Profiles to be plotted (fixed Z-mesh) --
	axialprofiles = config.chemistry.axial_profiles			# Axial 1D-Profiles to be plotted (fixed R-mesh) |
	probeloc = config.chemistry.probe_loc					# Cell location For Trend Analysis [R,Z], (leave empty for global min/max)

	# Various Diagnostic Settings			>>> OUTDATED, TO BE RETIRED <<<
	sheathROI = config.chemistry.sheath_ROI					# Sheath Region of Interest, (Start,End) [cells]
	sourcewidth = config.chemistry.source_width				# Source Dimension at ROI, leave empty for auto. [cells]


	# Requested diagnostics and plotting routines.
	savefig_tecplot2D = config.figures.tecplot2D			# 2D Single-Variables: TECPLOT2D.PDT				< .csv File Save

	savefig_movieicp2D = config.figures.movieicp2D			# 2D Variables against space-axis:	movie_icp.pdt	< MAXITER SHOULD BE AN ARRAY
	savefig_movieicp1D = config.figures.movieicp1D			# 1D Variables against space-axis:	movie_icp.pdt	< MAXITER SHOULD BE AN ARRAY
	savefig_timeaxis1D = config.figures.timeaxis1D			# 1D Variables against time-axis:	movie_icp.pdt
	savefig_convergence = config.figures.convergence		# 1D variables against ITER-axis:	movie_icp.pdt
	iterstep = config.figures.iterstep						# movie_icp.pdt iteration step size

	savefig_monoprofiles = config.figures.monoprofiles			# 1D Variables against space-axis:		TECPLOT2D	< .csv File Save
	savefig_multiprofiles = config.figures.multiprofiles		# 1D Variables Compared Same Sims:		TECPLOT2D
	savefig_compareprofiles = config.figures.compare_profiles	# 1D Variables Compared Between Sims:	TECPLOT2D
	# ^^^^
	# NOTE: IMAGEPLOTTER1D RETURNS ARRAY ORDERED AS [0,height]  <<< REVERSED RELATIVE TO HPEM
	# 		IMAGEPLOTTER2D RETURNS ARRAY ORDERED AS [height,0] 	<<< SAME ORIENTATION AS HPEM
	#		1D PROFILES ARE MANUALLY REVERSED ([::-1]) IN THE 1D DIAGNOSTICS TO ACCOUNT FOR THIS
	#		IT WOULD BE BETTER TO USE "DataExtent" TO PROVIDE THE CORRECT ORIENTATION

	savefig_trendphaseaveraged = config.figures.trend_phase_averaged	# Phase averaged trends at axial/radial cells		# CHANGE TO 'ProbeLoc' cell
	savefig_trendphaseresolved = config.figures.trend_phase_resolved	# Phase resolved trends at axial/radial cells		# CHANGE TO 'ProbeLoc' cell
	thrustloc = config.figures.thrust_loc								# Z-axis cell for thrust calculation  [Cells]

	savefig_phaseresolve2D = config.figures.phase_resolve2D		# 2D Phase Resolved Images							< .csv File Save
	savefig_phaseresolve1D = config.figures.phase_resolve1D		# 1D Phase Resolved Images							< .csv File Save
	savefig_sheathdynamics = config.figures.sheath_dynamics		# 1D and 2D sheath dynamics images
	savefig_PROES =	config.figures.PROES						# Simulated PROES Diagnostic
	phasecycles = config.figures.phase_cycles					# Vaveform phase cycles to be plotted. 				 [Float]
	DoFwidth = config.figures.do_Fwidth							# PROES Depth of Field (symmetric about image plane) [Cells]

	savefig_IEDFangular = config.figures.IEDF_angular			# 2D images of angular IEDF; single folders			< .csv File Save
	savefig_IEDFtrends = config.figures.IEDF_trends				# 1D IEDF trends; all folders
	savefig_EEDF = config.figures.EEDF							# 1D EEDF trends; all folders						< No Routine


	# Write processed data to ASCII files.
	write_ASCII = config.write.ASCII				# Data underpinning figs written in ASCII format	< Outdated
	Write_CSV = config.write.CSV					# Data underpinning figs written in .csv format
	# ^^^^
	# NOTE: SHEATH EXTENT SAVES WITH WRONG NAME, IN WRONG FILE FORMAT, IN ROOT DIRECTORY


	# Steady-State diagnostics terminal output toggles.
	print_generaltrends = config.printout.general_trends		# Verbose Min/Max Trend Outputs.
	print_Knudsennumber = config.printout.Knudsen_number		# Print cell averaged Knudsen Number
	print_totalpower = config.printout.total_power				# Print all requested total powers
	print_Reynolds = config.printout.Reynolds					# Print cell averaged sound speed
	print_DCbias = config.printout.DC_bias						# Print DC bias at electrodeloc
	print_thrust = config.printout.thrust						# Print neutral, ion and total thrust
	print_sheath = config.printout.sheath						# Print sheath width at electrodeloc


	# Image plotting options.
	image_extension = config.image.extension					# Define image extension  ('.png', '.jpg', '.eps')
	image_interp = config.image.interp							# Define image smoothing  ('none', 'bilinear','quadric','spline36')
	image_cmap = config.image.cmap								# Define global colourmap ('jet','plasma','inferno','gnuplot','tecmodern')

	image_aspectratio = config.image.aspect_ratio				# Real Size of [X,Y] in cm [Doesn't Rotate - X is always horizontal]
	image_radialcrop = config.image.radial_crop					# Crops 2D images to [R1,R2] in cm
	image_axialcrop = config.image.axial_crop					# Crops 2D images to [Z1,Z2] in cm

	image_plotsymmetry = config.image.plot_symmetry				# Plot radial symmetry - mirrors across the ISYM axis
	image_plotmesh = config.image.plot_mesh						# Plot material mesh outlines
	image_plotgrid = config.image.plot_grid						# Plot major/minor gridlines on 1D profiles
	image_rotate = config.image.rotate							# Rotate 2D images 90 degrees to the right.

	image_plotcolourfill = config.image.plot_colourfill			# Plot 2D image colour fill
	image_plotcontours = config.image.plot_contours				# Plot 2D image contour lines
	image_contourlvls = config.image.contour_lvls				# Number of contour levels

	image_axisticks = config.image.axis_ticks					# Toggle to show axis ticks and associated values or not
	image_axislabels = config.image.axis_labels					# Toggle to show axis labels or not
	image_legendloc = config.image.legend_loc					# Set Legend Location, "1-9" or 'best' for automatic
	image_cbarticks = config.image.cbar_ticks					# Toggle to show cbar ticks and associated values or not
	image_cbarbins = config.image.cbar_bins						# Set number of colourbar bins
	image_cbarlimit = config.image.cbar_limit					# Set arbitrary [min,max] colourbar limits

	image_plotvector = config.image.plot_vector					# Plot vector arrows onto 2D images (uses FR-XX, FZ-XX if they exist)
	image_vectordensity = config.image.vector_density			# Vector line density, higher means more dense streamlines
	image_vectorlw = config.image.vector_lw						# Vector line width, higher means thicker streamlines

	image_normalise = config.image.normalise					# Plot Data normlised to maximum value (Applies to all outputs)
	image_logplot = config.image.log_plot						# Plot log10(Data) (Applies to all outputs)

	image_plotsheath = config.image.plot_sheath					# Plot sheath extent onto 2D images 'Axial','Radial'
	image_plotphasewaveform = config.image.plot_phase_waveform	# Plot waveform sub-figure on phaseresolve2D images
	image_plotoverlay = config.image.plot_overlay				# Plot location(s) of 1D radial/axial profiles onto 2D images


	# Image Overrides (Applies to all images)
	titleoverride = config.overrides.title
	legendoverride = config.overrides.legend
	xaxisoverride = config.overrides.xaxis
	xlabeloverride = config.overrides.xlabel
	ylabeloverride = config.overrides.ylabel

	#============================#




	#####TODO#####

	### Immediate ###

	#
	#	NEED PERFORM ALL VARIABLE CONVERSION AND ROTATION AT THE READIN STAGE
	#		- SAVE "ORIENTATION" AND "ROTATION" AS GLOBAL VARIABLES THAT CAN BE USED IN DIAGNOSTICS
	#
	# 	NOTE: 	IMAGEPLOTTER1D RETURNS ARRAY ORDERED AS [0,height]  <<< !!! REVERSED RELATIVE TO HPEM !!!
	# 			IMAGEPLOTTER2D RETURNS ARRAY ORDERED AS [height,0] 	<<< SAME ORIENTATION AS HPEM
	#			1D PROFILES ARE MANUALLY REVERSED ([::-1]) IN THE 1D DIAGNOSTICS TO ACCOUNT FOR THIS
	#			IT WOULD BE BETTER TO USE "DataExtent" TO PROVIDE THE CORRECT ORIENTATION
	#
	# 	STILL HAVE ISSUE WITH READING movie FILES OF DIFFERENT LENGTHS...
	#		NEED TO SIZE THE INITIAL DATA ARRAY TO THE SHAPE OF THE LARGEST MOVIE FILE, NOT THE FIRST.
	#		ValueError: could not broadcast input array from shape (100,41,8804) into shape (24,41,8804)
	#
	#	INT SHEATH METHOD NOT COMPLETED FOR AXIAL SHEATH CALCULATION
	#		NOTE, RADIAL SHEATH CALC SCANS FROM R=0 INCREASING, I.E ASSUMES BULK IS ON-AXIS AT R=0 (mostly true)
	#		NOTE, AXIAL SHEATH CALC SCANS FROM Z=0 INCREASING< I.E. ASSUMES BULK IS NOT AT Z=0 (always true)
	#		image_plotsheath DOESN'T WORK VERY WELL WITH PROES IMAGES... NEEDS SPECIAL CASES!!!
	#		SHEATHCALC FUNCTION IS A LITTLE DODGY WHEN image_plotsheath IS NOT IN 'radial' or 'axial'
	#			NEED TO ENSURE THAT IT RETURNS NaN ARRAYS OF PROPER SIZE WHEN NOT CALLED EXPLICITLY
	#
	#	ANY VARIABLE NAMES CONTAINING "PHASE" ARE ALL READ IN TOGETHER
	#		NEED TO DIFFERENTIATE BETWEEN E, B PHASES AND ALSO R,Z,THETA PHASES
	#		CTRL+F FOR RM SJD AND SCROLL DOWN TO SEE NOTES IN VARIABLE NAMING FUNCTION
	#
	#	NEED TO FLIP THE SOUND SPEED DIAGNOSTIC, AND CHECK ORIENTATION OF THE KNUDSEN DIAGNOSTIC
	#		- ADD KNUDSEN, REYNOLDS, ETC... TO SEPARATE FUNCTION CALLS
	#
	#	SAVEFIG_SHEATHDYNAMICS USES THE ABSOLUTE SHEATH LOCATION, RATHER THAN THE SHEATH EXTENT
	#		NEED TO DEFINE THE ELECTRODE EDGE, REMOVE AND UPDATE 'sourcewidth' INPUT
	#		PERHAPS: "SheathZero" OR SOMETHING SIMILAR. DEFINE CELL TO USE AS SHEATH ZERO POINT
	#		THIS WILL THEN NEED TO MULTIPLY BY dr[l] or dz[l] DEPDNDING ON image_plotsheath ORIENTATION
	#
	#
	#
	#	NEED TO ADD ICOILP-n TOT OPTION WITH ALL COILSETS OVERLAYED
	#
	#	NEED to examine all maths in all trends to check for CGS or SI units
	#		- Each "trend" that employs a calculation will need an initial branch to
	#			check for global "Units" string, and convert accordingly.
	#		- Maths will always be performed in SI, so always convert back to SI for
	#			each current and future "Units" option.
	#

	#Variable Interpolator needs to work with phasedata - Take variables from batch?

	#Refactor PROES into 3D array (R,Z,Phase) and perform slice/integration
	#Movie1 data is flipped in the axial direction (for some reason).
	#Typically data is saved with Z=0 at top left (as is HPEM internal definition)
	#But Movie1 saves data with Z=0 at the bottom left
	#The following ad-hoc fixes have been applied
		#1DPhaseMovie applies [::-1] reversal on ExtractAxialProfile
		#2PhaseMovie applies [::-1] reversal on ImageExtractor2D
		#PROES applies [::-1] reversal on ExtractAxialProfile
	#IT WOULD BE GOOD IF THESE COULD BE APPLIED AT THE PHASE READIN STAGE
	#KEEP THE MESH DEFINITION AS Z=0 AT THE TOP LEFT, TO BE HOMOGENIOUS WITH HPEM.

	#Sheathwidth function integrates axially or radially depending on mesh geometry.
	#Sheathwidth function has 1D (Scott) and 2D (Greg) capabilities
	#SheathWidth function can deal with image rotations.

	#IEDF diagnostic capable of comparing between different material surfaces in single image
	#IEDF diagnostic saves different material surfaces in different folders

	#Thrust diagnostic split into functions performing the same task as before.
	#Thrust diagnostic enforces image symmetry, correcting the half-thrust error.

	#Add 'probeloc' option to savefig_temporal, where the meshavg is the default
	#but if there are any supplied probeloc cells then those get saved into their own seperate folder.

	#Fix Rynolds diagnostic - Current version has issue with NaNs and incorrect averaging
	#implement Rynolds number properly (was converted from sound speed diagnostic)
	#implement magnetic Rynolds number as well if possible






	### For Future ###

	#implement tests for key I/O functions
	#Implement tests for key mathematical/diagnostic functions (dc self-bias/thrust/soundspeed/knudsen/etc...)
	#Implement tests for each diagnostic
		#AND/OR include a tiny parallel plate folder in github with output figures for comparison

	#Add EEDF section and functionalise.

	#Fix issue with ffmpeg "convert-im6.q16: DistributedPixelCache"
	#	- Address the ignorance of os.system.

	#include 'garbage ./collection' at the end of each diagnostic.
	#Clean up unused functions and ensure homogeneity.
	#Split into modules: HELENAIO, HELENADIAGNOSTIC, HELENAPLOTTING
	#Update README, include all diagnostics and examples.
	#Create developer handbook describing functions.

	#Include Andor, OceanOptics, and LeCroy readin functions.














	#====================================================================#
							#INITIATE GLOBAL LISTS#
	#====================================================================#

	#Variables and lists for basic processing
	#Dir = list()					#List of all datafile directories inside folders
	#Dirlist = list()				#List of all folder directories (not including filenames)
	#IEDFVariableStrings = list()	#List of all variable names in pcmc.prof in header order
	#Geometrylist = list()			#List containing commonly used geometries [LEGACY: NOT USED]

	Globalvarlist = list()			#List of all commonly shared variable names between all folders
	Globalnumvars = list()			#Number of commonly shared variables between all folders.

	#Variables for mesh_size lists and SI conversion
	#ISYMlist = list()				#list of ISYM values in folder order in Dirlist
	#IXZlist = list()				#List of IXZ values in folder order in Dirlist
	#R_mesh = list()					#List of radial mesh cells for initmesh.out in folder order in Dirlist
	#Z_mesh = list()					#List of axial mesh cells for initmesh.out in folder order in Dirlist
	#Raxis = list()					#Radial SI [cm] axis for plotting
	#Zaxis = list()					#Axial SI [cm] axis for plotting

	#Depth = list()					#icp.nam Depth input [cm] in folder order in Dirlist
	#Radius = list()					#icp.nam Radius input [cm] in folder order in Dirlist
	#Height = list()					#icp.nam Height input [cm] in folder order in Dirlist
	#dr = list()						#Radial mesh resolution [cm/cell] in folder order in Dirlist
	#dz = list()						#Axial mesh resolution [cm/cell] in folder order in Dirlist
	#dy = list()						#Depth mesh resolution [cm/cell] in folder order in Dirlist

	#Variables and lists for icp.nam parameters
	#VRFM,VRFM2   = list(),list()							# Array of reals (In Material Order)
	#FREQM,FREQM2 = list(),list()							# Array of reals (In Material Order)
	#FREQC        = list()									# Array of reals (In Material Order)
	#FREQMAX,FREQMIN  = list(),list()						# Array of reals (In Material Order)
	#FREQGLOB,FREQALL = list(),list()						# real
	#IRFPOW = list()										# real
	#PRESOUT = list()										# real
	#IMOVIE_FRAMES = list()									# real
	#NUMMETALS=0; CMETALS,IETRODEM = list(),list()			# int;	Array of strings and ints, respectively
	#NUMCOILS=0; CCOILS = list()								# int;	Array of strings
	IMATSTATS=0#; CMATSTATS = list()							# int;	Array of strings
	IPCMCSPEC=0#; CPCMCSPEC = list()							# int;	Array of strings
	IEBINSPCMC=0; EMAXIPCMC=0								# int; 	int

	#Lists for icp.dat variables
	#header_icpdat = list()			#[SpeciesName, Charge, MolecularWeight, StickingCoeff,	- Array of strings
									# Transport, ReturnFrac, ReturnName]
	#AtomicSpecies = list()			#All species contained within chemistry set				- Array of strings
	#FluidSpecies  = list() 			#All neutral fluid species (for fluid dynamics) 		- Array of strings
	#NeutSpecies	= list()			#All neutral and metastable species						- Array of strings
	#PosSpecies = list()				#All positive ion species								- Array of strings
	#NegSpecies = list()				#All negative ion species								- Array of strings

	#Lists to store raw data
	rawdata_2D = list()				#ASCII format TECPLOT2D.pdt data string list			- Variable,Radius,Axis
	rawdata_kin = list()			#ASCII format kin.pdt data string list					- Variable,Radius,Axis
	rawdata_phasemovie = list()		#ASCII format movie1.pdt data string list				- Variable,Radius,Axis
	rawdata_itermovie = list()		#ASCII format movie_icp.pdt data string list 	  		- Variable,Radius,Axis
	rawdata_IEDF = list()			#ASCII format iprofile_tec2d.pdt data string list 		- Variable,Radius,Axis
	rawdata_mcs = list()			#ASCII format mcs.pdt data string list 			  		- Variable,Radius,Axis

	Data = list()					#Data[folder][Variable][Data]						-Data = 2D (R,Z) of reals
	DataKin = list()				#Data[folder][Variable][Data]						-Data = 1D (Avg) of reals
	DataIEDF = list()				#Data[folder][Variable][Data]						-Data = 2D (R,Z) of reals
	DataEEDF = list()				#Data[folder][Variable][Data]						-Data = 2D (R,Z) of Reals
	IterMovieData = list()			#ITERMovieData[folder][Timestep][Variable][Data]	-Data = 2D (R,Z) of Reals
	PhaseMovieData = list()			#PhaseMovieData[folder][Timestep][Variable][Data]	-Data = 2D (R,Z) of Reals

	Moviephaselist = list()			#'CYCL = n'												-int
	MovieIterlist = list()			#'ITER = n'												-int
	EEDF_TDlist = list()			#'???'													-???

	header_itermovie = list()		#Header num rows for movie_icp.pdt						-1D array of ints
	header_phasemovie = list()		#Header num rows for movie1.pdt							-1D array of ints
	header_IEDFlist = list()		#Header num rows for iprofile_tec2d.pdt					-1D array of ints
	header_kinlist = list()			#Header num rows for kin.pdt							-1D array of ints
	header_2Dlist = list()			#Header num rows for TECPLOT2D.pdt						-1D array of ints

	Header_TEC2D = list()			# TECPLOT2D Header Strings								-[FolderIdx,[Header Array]]
	Header_KIN = list()				# TECPLOT2D Header Strings								-[FolderIdx,[Header Array]]
	Header_movieicp = list()		# movie_icp Header Strings								-[FolderIdx,[Header Array]]
	Header_movie1 = list()			# movie1 Header Strings									-[FolderIdx,[Header Array]]
	Header_IEDF = list()			# IEDF Header Strings									-[FolderIdx,[Header Array]]





	#====================================================================#
						#WELCOME SPLASH AND INFORMATION#
	#====================================================================#

	print( '' )
	print( '--------------------------------------------------------------------')
	print( '    __    __   _______  __       _______  __   __      ___          ')
	print( '   |  |  |  | |   ____||  |     |   ____||  \ |  |    /   \         ')
	print( '   |  |__|  | |  |__   |  |     |  |__   |   \|  |   /  ^  \        ')
	print( '   |   __   | |   __|  |  |     |   __|  |  . `  |  /  /_\  \       ')
	print( '   |  |  |  | |  |____ |  `----.|  |____ |  |\   | /  _____  \      ')
	print( '   |__|  |__| |_______||_______||_______||__| \__|/__/     \__\     ')
	print( '                                                              v3.1.4')
	print( '--------------------------------------------------------------------')
	print( '' )
	print( 'The following diagnostics were requested:')
	print( '-----------------------------------------')
	if True in [savefig_tecplot2D]:
		print('# Plot 2D TECPLOT2D.PDT images')
	if True in [savefig_monoprofiles, savefig_multiprofiles, savefig_compareprofiles]:
		print('# Plot 1D TECPLOT2D.PDT images')
	if True in [savefig_movieicp2D]:
		print('# Plot 2D movieicp.pdt images')
	if True in [savefig_movieicp1D, savefig_convergence, savefig_timeaxis1D]:
		print('# Plot 1D movieicp.pdt images')
	if True in [savefig_phaseresolve2D, savefig_PROES]:
		print('# Plot 2D movie1.pdt images')
	if True in [savefig_phaseresolve1D, savefig_sheathdynamics]:
		print('# Plot 1D movie1.pdt images')
	if True in [savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]:
		print('# Plot PCMC distribution functions')
	if True in [savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_Reynolds, print_totalpower, print_DCbias, print_thrust]:
		print('# Plot Trends Between Folders')
	#endif

	print( '-----------------------------------------')
	print( '')





	#====================================================================#
						#OBTAINING FILE DIRECTORIES#
	#====================================================================#

	#Obtain system RAM. (and rename enviroment variable)
	mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
	mem_gib = mem_bytes/(1024.**3)
	ext = image_extension

	#Determine number of folders containing accepted file extensions (i.e. simulation folders)
	numfolders, Dir, Dirlist, HomeDir = get_directories()

	#Extract directories for all required data I/O files
	#These directories are relative to HELENA.py directory
	icpnam = [x for x in Dir if "icp.nam" in x]
	icpdat = [x for x in Dir if "icp.dat" in x]
	icpout = [x for x in Dir if "icp.out" in x]
	mesh = [x for x in Dir if "initmesh.out" in x]
	TEC2D = [x for x in Dir if "TECPLOT2D." in x]
	movieicp = [x for x in Dir if "movie_icp." in x]
	movie1 = [x for x in Dir if "movie1." in x]
	iprofiletec2d = [x for x in Dir if "iprofile_tec2d" in x]
	boltztec = [x for x in Dir if "boltz_tec" in x]

	R_mesh, Z_mesh, \
		ISYMlist, IXZlist, \
		Radius, Height, Depth, \
		dr, dz, dy, \
		VRFM, VRFM2, \
		FREQM, FREQM2, \
		FREQC, FREQGLOB, FREQALL, \
		FREQMIN, FREQMAX, \
		IRFPOW, \
		PRESOUT, \
		IMOVIE_FRAMES, \
		header_icpdat, \
		AtomicSpecies, FluidSpecies, NeutSpecies, PosSpecies, NegSpecies = get_mesh_and_SI(numfolders,
																						   Dirlist,
																						   icpnam,
																						   icpdat,
																						   TEC2D,
																						   mesh,
																						   Magmesh,
																						   image_plotsymmetry,
																						   IDEBUG,
																						   )

	#==========##========================##==========#
	#==========##========================##==========#




	#====================================================================#
					#UNPACKING AND ORGANIZATION OF DATA#
	#====================================================================#

	def PlotGeometry(ax,MeshCoordinates,MeshConnections,image_plotsymmetry,LabelNodes=False):
	#	Plots mesh geometry onto a 2D plot
	#	MeshCoordinates: 2D array of [[radius,height],[radius,height],...]		:: [cm]
	#					 List of mesh node (vertice) coordinates as floats
	#					 Expects origin to be top-left corner
	#	MeshConnections: 2D array of [[i, j],[i, j],...] 						:: [-]
	#					 list of indices that identify connections between nodes
	#					 indices are sequential to MeshCoordinates.
	###############

		# Convert strings in scientific notation to floats
		Coords = [(float(Radius), float(Height)) for Radius, Height in MeshCoordinates]

		# Plot nodes if requested (diagnostic tool)
		if LabelNodes == True:
			for idx, (Radius, Height) in enumerate(coords, start=1):
				ax.scatter(Radius, Height, color="red", s=20)
				ax.text(Radius, Height, f"{idx}", fontsize=8, ha="right", va="bottom")
			#endfor
		#endif

		# Draw connections
		for n in range(0,len(MeshConnections)):

			# Sequentially read mesh connection indices
			Conn = MeshConnections[n]

			# Adjust for zero-based index
			i = int(Conn[0]) - 1
			j =  int(Conn[1]) - 1

			# Extract coordinates at relevant connection indices
			Radius1, Height1 = Coords[i]
			Radius2, Height2 = Coords[j]

			# Determine if line to be plotted is on symmetry axis
			# If so, 'break' and skip to next 'n' to avoid plotting line at symmetry axis
			if image_plotsymmetry == True:

				# Get max height from current directory [l]
				Height_Max = round( dz[l]*(Z_mesh[l]-1),2)

				# Check if coordinate one is at the origin
				if Radius1 == 0.0 and round(Height1,1) == 0.0:
					# Check if either corresponding height coordinates are at maximum height
					if abs(Height_Max - Height1) <= dz[l] or abs(Height_Max - Height2) <= dz[l]:
						# If both are true, line is on symmetry axis
						# Skip current line plotting
						break
					#endif
				#endif

				# Check if coordinate two is at the origin
				if Radius2 == 0.0 and round(Height2,1) == 0.0:
					# Check if either corresponding height coordinates are at maximum height
					if abs(Height_Max - Height1) <= dz[l] or abs(Height_Max - Height2) <= dz[l]:
						# If both are true, line is on symmetry axis
						# Skip current line plotting
						break
					#endif
				#endif
			#endif

			#=====#=====#

			# Only works if rotating 90 to the right
			if image_rotate == True and image_plotsymmetry == True:
				Radius1,Height1 = Height1,Radius1
				Radius2,Height2 = Height2,Radius2
				ax.plot([Radius1, Radius2], [Height1, Height2], color='dimgrey', lw=2)
				ax.plot([Radius2, Radius1], [-Height1, -Height2], color='dimgrey', lw=2)

			# Only works if rotating 90 to the right
			elif image_rotate == True and image_plotsymmetry == False:
				Radius1,Height1 = Height1,Radius1
				Radius2,Height2 = Height2,Radius2
				Height1,Height2 = -Height2+R_mesh[l]*dr[l], -Height1+R_mesh[l]*dr[l]
				ax.plot([Radius1, Radius2], [Height1, Height2], color='dimgrey', lw=2)

			# Works for all cases
			elif image_rotate == False and image_plotsymmetry == True:
				ax.plot([Radius1, Radius2], [Height1, Height2], color='dimgrey', lw=2)
				ax.plot([-Radius2, -Radius1], [Height1, Height2], color='dimgrey', lw=2)

			# if no rotation or symmetry
			else:
				ax.plot([Radius1, Radius2], [Height1, Height2], color='dimgrey', lw=2)
			#endif
		#endfor

		return()
	#enddef

	#=============#

	def ManualPRCCPMesh(Ax=plt.gca()):
		#Plot pocket rocket material dimensions.
		Ax.plot((27*dz[l],27*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((75*dz[l],75*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((27*dz[l],27*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((75*dz[l],75*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((27*dz[l],75*dz[l]),   ( 0.21, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((27*dz[l],75*dz[l]),   (-0.21,-0.21), '-', color='dimgrey', linewidth=2)

		#Alumina Dielectric
		Ax.plot((34*dz[l],74*dz[l]),  ( 0.21,  0.21), 'c-', linewidth=2)
		Ax.plot((34*dz[l],74*dz[l]),  (-0.21, -0.21), 'c-', linewidth=2)
		Ax.plot((34*dz[l],74*dz[l]),  ( 0.31,  0.31), 'c-', linewidth=2)
		Ax.plot((34*dz[l],74*dz[l]),  (-0.31, -0.31), 'c-', linewidth=2)

		#Macor Dielectric
		Ax.plot((34*dz[l],34*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
		Ax.plot((34*dz[l],34*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
		Ax.plot((58*dz[l],58*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
		Ax.plot((58*dz[l],58*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
		Ax.plot((34*dz[l],58*dz[l]),  ( 47*dr[l], 47*dr[l]), 'b-', linewidth=2)
		Ax.plot((34*dz[l],58*dz[l]),  (-47*dr[l],-47*dr[l]), 'b-', linewidth=2)

		#Powered Electrode
		Ax.plot((39*dz[l],49*dz[l]),  ( 25*dr[l], 25*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],39*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],39*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
		Ax.plot((49*dz[l],49*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
		Ax.plot((49*dz[l],49*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  ( 43*dr[l], 43*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-43*dr[l],-43*dr[l]), 'r-', linewidth=2)

		#Grounded electrodes
		Ax.plot((34*dz[l],34*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((34*dz[l],34*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((74*dz[l],74*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((74*dz[l],74*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	#enddef

	#=============#

	def ManualPRCCPMMesh(Ax=plt.gca()):
		#Plot pocket rocket material dimensions.
		Ax.plot((27*dz[l],27*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((75*dz[l],75*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((27*dz[l],27*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((75*dz[l],75*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((27*dz[l],75*dz[l]),   ( 0.21, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((27*dz[l],75*dz[l]),   (-0.21,-0.21), '-', color='dimgrey', linewidth=2)

		#Alumina Dielectric
		Ax.plot((34*dz[l],74*dz[l]),  ( 0.21,  0.21), 'c-', linewidth=2)
		Ax.plot((34*dz[l],74*dz[l]),  (-0.21, -0.21), 'c-', linewidth=2)
		Ax.plot((34*dz[l],74*dz[l]),  ( 0.31,  0.31), 'c-', linewidth=2)
		Ax.plot((34*dz[l],74*dz[l]),  (-0.31, -0.31), 'c-', linewidth=2)

		#Macor Dielectric
		Ax.plot((34*dz[l],34*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
		Ax.plot((34*dz[l],34*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
		Ax.plot((58*dz[l],58*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
		Ax.plot((58*dz[l],58*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
		Ax.plot((34*dz[l],58*dz[l]),  ( 47*dr[l], 47*dr[l]), 'b-', linewidth=2)
		Ax.plot((34*dz[l],58*dz[l]),  (-47*dr[l],-47*dr[l]), 'b-', linewidth=2)

		#Powered Electrode
		Ax.plot((39*dz[l],49*dz[l]),  ( 25*dr[l], 25*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],39*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],39*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
		Ax.plot((49*dz[l],49*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
		Ax.plot((49*dz[l],49*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  ( 43*dr[l], 43*dr[l]), 'r-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-43*dr[l],-43*dr[l]), 'r-', linewidth=2)

		#Co-axial magnetic ring
		Ax.plot((39*dz[l],49*dz[l]),  ( 48*dr[l], 48*dr[l]), 'g-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-48*dr[l]), 'g-', linewidth=2)
		Ax.plot((39*dz[l],39*dz[l]),  ( 48*dr[l], 62*dr[l]), 'g-', linewidth=2)
		Ax.plot((39*dz[l],39*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)
		Ax.plot((49*dz[l],49*dz[l]),  ( 48*dr[l], 62*dr[l]), 'g-', linewidth=2)
		Ax.plot((49*dz[l],49*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-48*dr[l]), 'g-', linewidth=2)
		Ax.plot((39*dz[l],49*dz[l]),  (-62*dr[l],-62*dr[l]), 'g-', linewidth=2)

		#Grounded electrodes
		Ax.plot((34*dz[l],34*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((34*dz[l],34*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((74*dz[l],74*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
		Ax.plot((74*dz[l],74*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	#enddef

	#=============#

	def ManualHyperionIMesh(Ax=plt.gca()):
		#Plot upstream ICP material dimensions.
		Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68+40)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68-39)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((5*dz[l],51*dz[l]),    ((68+40)*dr[l],(68+40)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((5*dz[l],51*dz[l]),    ((68-39)*dr[l],(68-39)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((51*dz[l],51*dz[l]),   ((68+40)*dr[l],(68+65)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((51*dz[l],51*dz[l]),   ((68-39)*dr[l],(68-64)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((51*dz[l],111*dz[l]),  ((68+65)*dr[l],(68+65)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((51*dz[l],111*dz[l]),  ((68-64)*dr[l],(68-64)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((111*dz[l],111*dz[l]), (68*dr[l],(68+65)*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((111*dz[l],111*dz[l]), (68*dr[l],(68-64)*dr[l]), '-', color='lightgrey', linewidth=4)
		#4cm = 74*dz[l], 12cm=111*dz[l]

		#Macor Dielectric
		Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68+40)*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68-39)*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((5*dz[l],51*dz[l]),    ((68+40)*dr[l],(68+40)*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((5*dz[l],51*dz[l]),    ((68-39)*dr[l],(68-39)*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((51*dz[l],51*dz[l]),   ((68+40)*dr[l],(68+65)*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((51*dz[l],51*dz[l]),   ((68-39)*dr[l],(68-64)*dr[l]), '-', color='c', linewidth=4)

		#Powered Electrode - 'Metal' Anode
	#	Ax.plot((65*dz[l],65*dz[l]),   (2*dr[l],70*dr[l]), '-', color='red', linewidth=4)
	#	Ax.plot((65*dz[l],65*dz[l]),   (-2*dr[l],-70*dr[l]), '-', color='red', linewidth=4)

		#Powered ICP Coils - 'Metal'
		Ax.plot((12*dz[l],16*dz[l]),   (114*dr[l],114*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((12*dz[l],16*dz[l]),   (122*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((12*dz[l],12*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((16*dz[l],16*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Downstream
		Ax.plot((7*dz[l],11*dz[l]),    (23*dr[l],23*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((7*dz[l],11*dz[l]),    (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((7*dz[l],7*dz[l]),     (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((11*dz[l],11*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Downstream

		Ax.plot((22*dz[l],26*dz[l]),   (114*dr[l],114*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((22*dz[l],26*dz[l]),   (122*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((22*dz[l],22*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((26*dz[l],26*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Downtream
		Ax.plot((17*dz[l],21*dz[l]),   (23*dr[l],23*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((17*dz[l],21*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((17*dz[l],17*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((21*dz[l],21*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Downstream

		Ax.plot((32*dz[l],36*dz[l]),   (114*dr[l],114*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((32*dz[l],36*dz[l]),   (122*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((32*dz[l],32*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((36*dz[l],36*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Downstream
		Ax.plot((27*dz[l],31*dz[l]),   (23*dr[l],23*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((27*dz[l],31*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((31*dz[l],31*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((27*dz[l],27*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Downstream

		#Horseshoe Solenoid
		Ax.plot((44*dz[l],48*dz[l]),   (124*dr[l],124*dr[l]), '-', color='lightgreen', linewidth=5)	#Inboard
		Ax.plot((44*dz[l],48*dz[l]),   (13*dr[l],13*dr[l]), '-', color='lightgreen', linewidth=5)	#Inboard
		Ax.plot((44*dz[l],48*dz[l]),   (132*dr[l],132*dr[l]), '-', color='lightgreen', linewidth=5)	#Outboard
		Ax.plot((44*dz[l],48*dz[l]),   (5*dr[l],5*dr[l]), '-', color='lightgreen', linewidth=5)		#Outboard
		Ax.plot((44*dz[l],44*dz[l]),   (124*dr[l],132*dr[l]), '-', color='lightgreen', linewidth=5)	#Upstream
		Ax.plot((48*dz[l],48*dz[l]),   (124*dr[l],132*dr[l]), '-', color='lightgreen', linewidth=5)	#Downstream
		Ax.plot((44*dz[l],44*dz[l]),   (13*dr[l],5*dr[l]), '-', color='lightgreen', linewidth=5)	#Upstream
		Ax.plot((48*dz[l],48*dz[l]),   (13*dr[l],5*dr[l]), '-', color='lightgreen', linewidth=5)	#Downstream
	#enddef

	def ManualHyperionIIMesh(Ax=plt.gca()):
		#Plot upstream ICP material dimensions.
		Ax.plot((1*dz[l],1*dz[l]),       (3*dr[l],109*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((1*dz[l],84*dz[l]),      (3*dr[l],3*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((1*dz[l],84*dz[l]),      (109*dr[l],109*dr[l]), '-', color='lightgrey', linewidth=4)
		Ax.plot((84*dz[l],84*dz[l]),     (3*dr[l],109*dr[l]), '-', color='lightgrey', linewidth=4)

		#Alumina Dielectric
		Ax.plot((2*dz[l],2*dz[l]),     (25*dr[l],87*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((2*dz[l],84*dz[l]),    (25*dr[l],25*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((2*dz[l],84*dz[l]),    (87*dr[l],87*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((14*dz[l],14*dz[l]),   (37*dr[l],75*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((14*dz[l],84*dz[l]),   (37*dr[l],37*dr[l]), '-', color='c', linewidth=4)
		Ax.plot((14*dz[l],84*dz[l]),   (75*dr[l],75*dr[l]), '-', color='c', linewidth=4)

		#Moly Extraction Aperture
		Ax.plot((71*dz[l],84*dz[l]),   (38*dr[l],38*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((71*dz[l],84*dz[l]),   (41*dr[l],41*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((71*dz[l],71*dz[l]),   (38*dr[l],41*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((73*dz[l],73*dz[l]),   (42*dr[l],53*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((74*dz[l],74*dz[l]),   (42*dr[l],53*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((73*dz[l],74*dz[l]),   (53*dr[l],53*dr[l]), '-', color='m', linewidth=4)

		Ax.plot((71*dz[l],84*dz[l]),   (74*dr[l],74*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((71*dz[l],84*dz[l]),   (70*dr[l],70*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((71*dz[l],71*dz[l]),   (70*dr[l],74*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((73*dz[l],73*dz[l]),   (58*dr[l],70*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((74*dz[l],74*dz[l]),   (58*dr[l],70*dr[l]), '-', color='m', linewidth=4)
		Ax.plot((73*dz[l],74*dz[l]),   (58*dr[l],58*dr[l]), '-', color='m', linewidth=4)

		#Powered Electrode - 'Metal' Anode
	#	Ax.plot((65*dz[l],65*dz[l]),   (2*dr[l],70*dr[l]), '-', color='red', linewidth=4)
	#	Ax.plot((65*dz[l],65*dz[l]),   (-2*dr[l],-70*dr[l]), '-', color='red', linewidth=4)

		#Powered ICP Coils - 'Metal'
		Ax.plot((27*dz[l],27*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((29*dz[l],29*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((27*dz[l],29*dz[l]),   (9*dr[l],9*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((27*dz[l],29*dz[l]),   (14*dr[l],14*dr[l]), '-', color='red', linewidth=5)

		Ax.plot((33*dz[l],33*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((35*dz[l],35*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((33*dz[l],35*dz[l]),   (9*dr[l],9*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((33*dz[l],35*dz[l]),   (14*dr[l],14*dr[l]), '-', color='red', linewidth=5)

		Ax.plot((39*dz[l],39*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((41*dz[l],41*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((39*dz[l],41*dz[l]),   (9*dr[l],9*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((39*dz[l],41*dz[l]),   (14*dr[l],14*dr[l]), '-', color='red', linewidth=5)

		Ax.plot((30*dz[l],30*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((32*dz[l],32*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((30*dz[l],32*dz[l]),   (98*dr[l],98*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((30*dz[l],32*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)

		Ax.plot((36*dz[l],36*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((38*dz[l],38*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((36*dz[l],38*dz[l]),   (98*dr[l],98*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((36*dz[l],38*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)

		Ax.plot((42*dz[l],42*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((44*dz[l],44*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((42*dz[l],44*dz[l]),   (98*dr[l],98*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((42*dz[l],44*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((42*dz[l],44*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)

		#Horseshoe Solenoid
		Ax.plot((68*dz[l],68*dz[l]),   (4*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((70*dz[l],70*dz[l]),   (4*dr[l],15*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((68*dz[l],70*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((70*dz[l],72*dz[l]),   (15*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((72*dz[l],72*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((70*dz[l],72*dz[l]),   (22*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)

		Ax.plot((68*dz[l],68*dz[l]),   (4*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((70*dz[l],70*dz[l]),   (4*dr[l],15*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((68*dz[l],70*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((70*dz[l],72*dz[l]),   (15*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((72*dz[l],72*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((70*dz[l],72*dz[l]),   (22*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)

		#Iron Magnetic Focus
		Ax.plot((71*dz[l],71*dz[l]),   (43*dr[l],46*dr[l]), '-', color='g', linewidth=6.5)
		Ax.plot((72*dz[l],72*dz[l]),   (43*dr[l],53*dr[l]), '-', color='g', linewidth=6.5)

		Ax.plot((71*dz[l],71*dz[l]),   (65*dr[l],68*dr[l]), '-', color='g', linewidth=6.5)
		Ax.plot((72*dz[l],72*dz[l]),   (58*dr[l],68*dr[l]), '-', color='g', linewidth=6.5)
	#enddef

	#=============#

	def ManualEVgenyMesh(Ax=plt.gca()):
		#Plot upstream ICP material dimensions.
		Ax.plot((2.5*dz[l],2.5*dz[l]), (0*dr[l],20*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((2.5*dz[l],2.5*dz[l]), (-20*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((2.5*dz[l],41*dz[l]),  (20*dr[l],20*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((2.5*dz[l],41*dz[l]),  (-20*dr[l],-20*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((41*dz[l],41*dz[l]),   (20*dr[l],90*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((41*dz[l],41*dz[l]),   (-20*dr[l],-90*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((41*dz[l],87*dz[l]),   (90*dr[l],90*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((41*dz[l],87*dz[l]),   (-90*dr[l],-90*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((87*dz[l],87*dz[l]),   (0*dr[l],90*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((87*dz[l],87*dz[l]),   (-90*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)

		#Macor Dielectric
		Ax.plot((2.5*dz[l],2.5*dz[l]), (0*dr[l],20*dr[l]), 'c-', linewidth=4)
		Ax.plot((2.5*dz[l],2.5*dz[l]), (-20*dr[l],0*dr[l]), 'c-', linewidth=4)
		Ax.plot((3*dz[l],6*dz[l]),    (20*dr[l],20*dr[l]), 'c-', linewidth=4)
		Ax.plot((3*dz[l],6*dz[l]),    (-20*dr[l],-20*dr[l]), 'c-', linewidth=4)
		Ax.plot((23*dz[l],41*dz[l]),  (20*dr[l],20*dr[l]), 'c-', linewidth=4)
		Ax.plot((23*dz[l],41*dz[l]),  (-20*dr[l],-20*dr[l]), 'c-', linewidth=4)
		Ax.plot((41*dz[l],41*dz[l]),   (20*dr[l],90*dr[l]), 'c-', linewidth=4)
		Ax.plot((41*dz[l],41*dz[l]),   (-20*dr[l],-90*dr[l]), 'c-', linewidth=4)

		#Powered Electrode - LaB6 Cathode
		Ax.plot((6*dz[l],23*dz[l]),   (20*dr[l],20*dr[l]), '-', color='orange', linewidth=5)
		Ax.plot((6*dz[l],23*dz[l]),   (-20*dr[l],-20*dr[l]), '-', color='orange', linewidth=5)

		#Powered Electrode - 'Metal' Anode
		Ax.plot((42*dz[l],86*dz[l]),   (91*dr[l],91*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((42*dz[l],86*dz[l]),   (-91*dr[l],-91*dr[l]), '-', color='red', linewidth=5)

		#Powered ICP Coils - 'Metal'
		Ax.plot((6*dz[l],8*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((6*dz[l],8*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((6*dz[l],8*dz[l]),   (35*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((6*dz[l],8*dz[l]),   (-35*dr[l],-35*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((6*dz[l],6*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((6*dz[l],6*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((8*dz[l],8*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Downstream
		Ax.plot((8*dz[l],8*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)		#Downstream

		Ax.plot((13.5*dz[l],15.5*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((13.5*dz[l],15.5*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Inboard
		Ax.plot((13.5*dz[l],15.5*dz[l]),   (35*dr[l],35*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((13.5*dz[l],15.5*dz[l]),   (-35*dr[l],-35*dr[l]), '-', color='red', linewidth=5)#Outboard
		Ax.plot((13.5*dz[l],13.5*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((15.5*dz[l],15.5*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)	#Downtream
		Ax.plot((13.5*dz[l],13.5*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)#Upstream
		Ax.plot((15.5*dz[l],15.5*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)#Downstream

		Ax.plot((21*dz[l],23*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((21*dz[l],23*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((21*dz[l],23*dz[l]),   (35*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((21*dz[l],23*dz[l]),   (-35*dr[l],-35*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((21*dz[l],21*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((23*dz[l],23*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Downstream
		Ax.plot((21*dz[l],21*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((23*dz[l],23*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)	#Downstream

		#Solenoid
		Ax.plot((42*dz[l],86*dz[l]),   (93*dr[l],93*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((42*dz[l],86*dz[l]),   (-93*dr[l],-93*dr[l]), '-', color='lightgreen', linewidth=5)
	#enddef

	def ManualEVgenyMeshOLD(Ax=plt.gca()):
		#Plot upstream ICP material dimensions.
		Ax.plot((0.5*dz[l],0.5*dz[l]), (0*dr[l],10*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((0.5*dz[l],0.5*dz[l]), (-10*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((0.5*dz[l],39*dz[l]),  (10*dr[l],10*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((0.5*dz[l],39*dz[l]),  (-10*dr[l],-10*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((39*dz[l],39*dz[l]),   (10*dr[l], 80*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((39*dz[l],39*dz[l]),   (-10*dr[l], -80*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((39*dz[l],85*dz[l]),   (80*dr[l],80*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((39*dz[l],85*dz[l]),   (-80*dr[l],-80*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((85*dz[l],85*dz[l]),   (0*dr[l],80*dr[l]), '-', color='dimgrey', linewidth=4)
		Ax.plot((85*dz[l],85*dz[l]),   (-80*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)

		#Macor Dielectric
		Ax.plot((1*dz[l],4*dz[l]),    (10*dr[l],10*dr[l]), 'c-', linewidth=4)
		Ax.plot((1*dz[l],4*dz[l]),    (-10*dr[l],-10*dr[l]), 'c-', linewidth=4)
		Ax.plot((21*dz[l],39*dz[l]),  (10*dr[l],10*dr[l]), 'c-', linewidth=4)
		Ax.plot((21*dz[l],39*dz[l]),  (-10*dr[l],-10*dr[l]), 'c-', linewidth=4)

		#Powered Electrode - LaB6 Cathode
		Ax.plot((4*dz[l],21*dz[l]),   (10*dr[l],10*dr[l]), '-', color='orange', linewidth=5)
		Ax.plot((4*dz[l],21*dz[l]),   (-10*dr[l],-10*dr[l]), '-', color='orange', linewidth=5)

		#Powered Electrode - 'Metal' Anode
		Ax.plot((40*dz[l],84*dz[l]),   (81*dr[l],81*dr[l]), '-', color='red', linewidth=5)
		Ax.plot((40*dz[l],84*dz[l]),   (-81*dr[l],-81*dr[l]), '-', color='red', linewidth=5)

		#Powered ICP Coils - 'Metal'
		Ax.plot((4*dz[l],6*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((4*dz[l],6*dz[l]),   (-15*dr[l],-15*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((4*dz[l],6*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((4*dz[l],6*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((4*dz[l],4*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((4*dz[l],4*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((6*dz[l],6*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Downstream
		Ax.plot((6*dz[l],6*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Downstream

		Ax.plot((11.5*dz[l],13.5*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((11.5*dz[l],13.5*dz[l]),   (-15*dr[l],-15*dr[l]), '-', color='red', linewidth=5)#Inboard
		Ax.plot((11.5*dz[l],13.5*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((11.5*dz[l],13.5*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Outboard
		Ax.plot((11.5*dz[l],11.5*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((13.5*dz[l],13.5*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Downtream
		Ax.plot((11.5*dz[l],11.5*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Upstream
		Ax.plot((13.5*dz[l],13.5*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Downstream

		Ax.plot((19*dz[l],21*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Inboard
		Ax.plot((19*dz[l],21*dz[l]),   (-15*dr[l],-15*dr[l]), '-', color='red', linewidth=5)	#Inboard
		Ax.plot((19*dz[l],21*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Outboard
		Ax.plot((19*dz[l],21*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Outboard
		Ax.plot((19*dz[l],19*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Upstream
		Ax.plot((21*dz[l],21*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Downstream
		Ax.plot((19*dz[l],19*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Upstream
		Ax.plot((21*dz[l],21*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Downstream

		#Solenoid
		Ax.plot((40*dz[l],84*dz[l]),   (83*dr[l],83*dr[l]), '-', color='lightgreen', linewidth=5)
		Ax.plot((40*dz[l],84*dz[l]),   (-83*dr[l],-83*dr[l]), '-', color='lightgreen', linewidth=5)
	#enddef

	#=============#

	def ManualGECMesh(Ax,Sym=True,ICP_Powered=True):

	#	Notes, 	these apply to the 0.1R 0.2Z resolution mesh
	#			These are arranged for zero rotation, the other PR meshes assume rotation=True

		Thick = 4

		# Vacuum Chamber Walls (outer boundary)
		Ax.plot((0.0*dr,60.0*dr), (0.0*dz,0.0*dz), '-', color='dimgrey', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-60.0*dr), (0.0*dz,0.0*dz), '-', color='dimgrey', linewidth=Thick)
		Ax.plot((0.0*dr,60.0*dr), (119.0*dz,119.0*dz),  '-', color='dimgrey', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-60.0*dr), (119.0*dz,119.0*dz),  '-', color='dimgrey', linewidth=Thick)
		Ax.plot((60.0*dr,60.0*dr), (0.0*dz,119.0*dz),  '-', color='dimgrey', linewidth=Thick)
		if Sym == True: Ax.plot((-60.0*dr,-60.0*dr), (0.0*dz,119.0*dz),  '-', color='dimgrey', linewidth=Thick)

		# Vacuum Chamber Walls (inner boundary)
		Ax.plot((0.0*dr,58.0*dr), (4.0*dz,4.0*dz), '-', color='dimgrey', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-58.0*dr), (4.0*dz,4.0*dz), '-', color='dimgrey', linewidth=Thick)
		Ax.plot((0.0*dr,58.0*dr), (115.0*dz,115.0*dz),  '-', color='dimgrey', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-58.0*dr), (115.0*dz,115.0*dz),  '-', color='dimgrey', linewidth=Thick)
		Ax.plot((58.0*dr,58.0*dr), (4.0*dz,115.0*dz),  '-', color='dimgrey', linewidth=Thick)
		if Sym == True: Ax.plot((-58.0*dr,-58.0*dr), (4.0*dz,115.0*dz),  '-', color='dimgrey', linewidth=Thick)

		# Upper ICP Antenna Housing
		Ax.plot((35.0*dr,35.0*dr), (4.0*dz,32.0*dz), '-', color='dimgrey', linewidth=Thick)		#Outboard Vertical
		if Sym == True: Ax.plot((-35.0*dr,-35.0*dr), (4.0*dz,32.0*dz), '-', color='dimgrey', linewidth=Thick)	#Outboard Vertical
		Ax.plot((35.0*dr,41.0*dr), (32.0*dz,32.0*dz), '-', color='dimgrey', linewidth=Thick)	#Top of lip
		if Sym == True: Ax.plot((-35.0*dr,-41.0*dr), (32.0*dz,32.0*dz), '-', color='dimgrey', linewidth=Thick)	#Top of lip
		Ax.plot((41.0*dr,41.0*dr), (32.0*dz,43.0*dz), '-', color='dimgrey', linewidth=Thick)	#Outboard lip
		if Sym == True: Ax.plot((-41.0*dr,-41.0*dr), (32.0*dz,43.0*dz), '-', color='dimgrey', linewidth=Thick)	#Outboard lip
		Ax.plot((28.0*dr,41.0*dr), (43.0*dz,43.0*dz), '-', color='dimgrey', linewidth=Thick)	#Bottom of Lip
		if Sym == True: Ax.plot((-28.0*dr,-41.0*dr), (43.0*dz,43.0*dz), '-', color='dimgrey', linewidth=Thick)	#Bottom of lip
		Ax.plot((28.0*dr,28.0*dr), (43.0*dz,37.0*dz), '-', color='dimgrey', linewidth=Thick)	#Inboard Vertical
		if Sym == True: Ax.plot((-28.0*dr,-28.0*dr), (43.0*dz,37.0*dz), '-', color='dimgrey', linewidth=Thick)	#Inboard Vertical

		# Upper ICP Dielectric Window
		Ax.plot((0.0*dr,32.0*dr), (26.0*dz,26.0*dz), '-', color='orange', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-32.0*dr), (26.0*dz,26.0*dz), '-', color='orange', linewidth=Thick)
		Ax.plot((0.0*dr,32.0*dr), (37.0*dz,37.0*dz), '-', color='orange', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-32.0*dr), (37.0*dz,37.0*dz), '-', color='orange', linewidth=Thick)
		Ax.plot((32.0*dr,32.0*dr), (26.0*dz,37.0*dz), '-', color='orange', linewidth=Thick)
		if Sym == True: Ax.plot((-32.0*dr,-32.0*dr), (26.0*dz,37.0*dz), '-', color='orange', linewidth=Thick)

		# Upper ICP Antennae
		if ICP_Powered == True:
			Antenna1 = [1.5,3.0,22.0,25.0]			# R1, R2, Z1, Z2 [Cells]
			Antenna2 = [7.5,9.0,22.0,25.0]			# R1, R2, Z1, Z2 [Cells]
			Antenna3 = [13.5,15.0,22.0,25.0]		# R1, R2, Z1, Z2 [Cells]
			Antenna4 = [19.5,21.0,22.0,25.0]		# R1, R2, Z1, Z2 [Cells]
			Ants = [Antenna1,Antenna2,Antenna3,Antenna4]
			for i in range(0,4):
				Ax.plot((Ants[i][0]*dr,Ants[i][1]*dr), (22.0*dz,22.0*dz), '-', color='red', linewidth=Thick)
				if Sym == True: Ax.plot((-Ants[i][0]*dr,-Ants[i][1]*dr), (22.0*dz,22.0*dz), '-', color='red', linewidth=Thick)
				Ax.plot((Ants[i][0]*dr,Ants[i][1]*dr), (25.0*dz,25.0*dz), '-', color='red', linewidth=Thick)
				if Sym == True: Ax.plot((-Ants[i][0]*dr,-Ants[i][1]*dr), (25.0*dz,25.0*dz), '-', color='red', linewidth=Thick)
				Ax.plot((Ants[i][0]*dr,Ants[i][0]*dr), (22.0*dz,25.0*dz), '-', color='red', linewidth=Thick)
				if Sym == True: Ax.plot((-Ants[i][0]*dr,-Ants[i][0]*dr), (22.0*dz,25.0*dz), '-', color='red', linewidth=Thick)
				Ax.plot((Ants[i][1]*dr,Ants[i][1]*dr), (22.0*dz,25.0*dz), '-', color='red', linewidth=Thick)
				if Sym == True: Ax.plot((-Ants[i][1]*dr,-Ants[i][1]*dr), (22.0*dz,25.0*dz), '-', color='red', linewidth=Thick)
			#endfor
		#endif

		# Lower CCP Electrode Housing
		Ax.plot((28.0*dr,28.0*dr), (82.0*dz,115.0*dz), '-', color='dimgrey', linewidth=Thick)	#Outboard Vertical
		if Sym == True: Ax.plot((-28.0*dr,-28.0*dr), (82.0*dz,115.0*dz), '-', color='dimgrey', linewidth=Thick)	#Outboard Vertical
		Ax.plot((0.0*dr,28.0*dr), (82.0*dz,82.0*dz), '-', color='dimgrey', linewidth=Thick)		#Top of platform
		if Sym == True: Ax.plot((-0.0*dr,-28.0*dr), (82.0*dz,82.0*dz), '-', color='dimgrey', linewidth=Thick)	#Top of platform
		Ax.plot((25.0*dr,25.0*dr), (82.0*dz,115.0*dz), '-', color='dimgrey', linewidth=Thick)	#Inboard Vertical
		if Sym == True: Ax.plot((-25.0*dr,-25.0*dr), (82.0*dz,115.0*dz), '-', color='dimgrey', linewidth=Thick)	#Inboard Vertical

		# Lower CCP Electrode
		Ax.plot((0.0*dr,25.0*dr), (82.0*dz,82.0*dz), '-', color='red', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-25.0*dr), (82.0*dz,82.0*dz), '-', color='red', linewidth=Thick)
		Ax.plot((0.0*dr,25.0*dr), (97.0*dz,97.0*dz), '-', color='red', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-25.0*dr), (97.0*dz,97.0*dz), '-', color='red', linewidth=Thick)
		Ax.plot((25.0*dr,25.0*dr), (82.0*dz,97.0*dz), '-', color='red', linewidth=Thick)
		if Sym == True: Ax.plot((-25.0*dr,-25.0*dr), (82.0*dz,97.0*dz), '-', color='red', linewidth=Thick)

		# Wafer
		Ax.plot((0.0*dr,25.5*dr), (80.0*dz,80.0*dz), '-', color='cyan', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-25.5*dr), (80.0*dz,80.0*dz), '-', color='cyan', linewidth=Thick)
		Ax.plot((0.0*dr,25.5*dr), (81.0*dz,81.0*dz), '-', color='cyan', linewidth=Thick)
		if Sym == True: Ax.plot((-0.0*dr,-25.5*dr), (81.0*dz,81.0*dz), '-', color='cyan', linewidth=Thick)
	#enddef

	#=============#

	def ManualASTRONmk1Mesh(Ax=plt.gca()):

		#Line Plotting Format
		#	Ax.plot((Yend,Ystart),(Xend,Xstart))

		dZ = dz[l]*2
		dR = dr[l]*2

		#Plot Reactor Wall material dimensions.
		Ax.plot((1.5*dZ,1.5*dZ), (1.5*dR,17.5*dR), '-', color='grey', linewidth=4)
		Ax.plot((17.5*dZ,17.5*dZ), (1.5*dR,17.5*dR), '-', color='grey', linewidth=4)
		Ax.plot((1.5*dZ,17.5*dZ), (1.5*dR,1.5*dR), '-', color='grey', linewidth=4)
		Ax.plot((1.5*dZ,17.5*dZ), (17.5*dR,17.5*dR), '-', color='grey', linewidth=4)

		#Plot Magnetic Ferrous Core material dimensions.
		Ax.plot((8.0*dZ,8.0*dZ), (8.0*dR,11.0*dR), '-', color='green', linewidth=4)
		Ax.plot((11.0*dZ,11.0*dZ), (8.0*dR,11.0*dR), '-', color='green', linewidth=4)
		Ax.plot((8.0*dZ,11.0*dZ), (8.0*dR,8.0*dR), '-', color='green', linewidth=4)
		Ax.plot((8.0*dZ,11.0*dZ), (11.0*dR,11.0*dR), '-', color='green', linewidth=4)
		#Magnetic Material Cross-Hatch
	#	Ax.plot((11.0*dZ,8.0*dZ), (11.0*dR,8.0*dR), '-', color='green', linewidth=2)
	#	Ax.plot((8.0*dZ,11.0*dZ), (11.0*dR,8.0*dR), '-', color='green', linewidth=2)

		#Plot ICP Antenna Material Outline
		Ax.plot((6.0*dZ,6.0*dZ), (6.0*dR,13.0*dR), '-', color='red', linewidth=4)
		Ax.plot((13.0*dZ,13.0*dZ), (6.0*dR,13.0*dR), '-', color='red', linewidth=4)
		Ax.plot((6.0*dZ,13.0*dZ), (6.0*dR,6.0*dR), '-', color='red', linewidth=4)
		Ax.plot((6.0*dZ,13.0*dZ), (13.0*dR,13.0*dR), '-', color='red', linewidth=4)

		#Plot ICP Antenna Material Segmentation
		Ax.plot((6.0*dZ,8.0*dZ), (8.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((6.0*dZ,8.0*dZ), (8.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((6.0*dZ,8.0*dZ), (11.0*dR,11.0*dR), '-', color='red', linewidth=2)
		Ax.plot((6.0*dZ,8.0*dZ), (11.0*dR,11.0*dR), '-', color='red', linewidth=2)

		Ax.plot((11.0*dZ,13.0*dZ), (8.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,13.0*dZ), (8.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,13.0*dZ), (11.0*dR,11.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,13.0*dZ), (11.0*dR,11.0*dR), '-', color='red', linewidth=2)

		Ax.plot((8.0*dZ,8.0*dZ), (11.0*dR,13.0*dR), '-', color='red', linewidth=2)
		Ax.plot((8.0*dZ,8.0*dZ), (11.0*dR,13.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,11.0*dZ), (11.0*dR,13.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,11.0*dZ), (11.0*dR,13.0*dR), '-', color='red', linewidth=2)

		Ax.plot((8.0*dZ,8.0*dZ), (6.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((8.0*dZ,8.0*dZ), (6.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,11.0*dZ), (6.0*dR,8.0*dR), '-', color='red', linewidth=2)
		Ax.plot((11.0*dZ,11.0*dZ), (6.0*dR,8.0*dR), '-', color='red', linewidth=2)

	#enddef

	#===================##===================#
	#===================##===================#

































	#====================================================================#
						 #READING DATA INTO MEMORY#
	#====================================================================#

	print('----------------------')
	print('Beginning Data Readin:')
	print('----------------------')

	#Retrieve mesh geometry for plotting if requested.
	if image_plotmesh == True:

		# Create mesh data arrays
		MeshCoordinates = list()
		MeshConnections = list()

		# Geometry Name String :: Assume TECPLOT2D.PDT always exists
		NameString = "TECPLOT2D"

		# Extract geometry from supplied namestring
		DataFileDir = filter(lambda x: NameString in x, Dir)
		DataFileDir = sorted(DataFileDir)

		# Save mesh coordinates for each folder sequentially
		# MeshCoordinates[Folder][Nodes [R1,Z2], [R2,Z2], ...]
		for l in range (0,numfolders):
			GeomFileDir = DataFileDir[l].rsplit('/',1)[0] + "/meshnodes.dat"
			MeshCoord,MeshConn = read_geometry(DataFileDir[l], GeomFileDir, Z_mesh[l], dz[l])

			MeshCoordinates.append(MeshCoord)
			MeshConnections.append(MeshConn)
		#endfor
	#endif

	#===================##===================#
	#===================##===================#

	#Extraction and organization of data from .PDT files.
	for l in tqdm(range(0,numfolders)):

		#Load data from TECPLOT2D file and unpack into 1D array.
		try: rawdata, nn_2D = extract_raw_data(Dir, TEC2D[l].split('/')[-1], l)
		except: rawdata, nn_2D = extract_raw_data(Dir, 'TECPLOT2D.PDT', l)
		rawdata_2D.append(rawdata)

		#Read through all variables for each file and stop when list ends.
		TEC2DVariableStrings,HeaderEndMarker = ['Radius','Height'],'ZONE'
		for i in range(2,nn_2D):
			if HeaderEndMarker in str(rawdata_2D[l][i]): break
			else: TEC2DVariableStrings.append(str(rawdata_2D[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		Header_TEC2D.append(TEC2DVariableStrings)											# SJD USE THIS HEADER !!!
		numvariables_2D,header_2D = len(TEC2DVariableStrings),len(TEC2DVariableStrings)+2	# REMOVE THIS
		header_2Dlist.append(header_2D)														# REMOVE THIS

		#Seperate total 1D data array into sets of data for each variable.
		CurrentFolderData = read_TEC2D(rawdata_2D[l], header_2D, numvariables_2D, R_mesh[l], Z_mesh[l])

		#Convert data from CGS (HPEM DEFAULT) to user requested unit system
		for i in range(0,len(TEC2DVariableStrings)):
			CurrentFolderData[i] = variable_unit_conversion(CurrentFolderData[i], TEC2DVariableStrings[i], Units, AtomicSpecies)
	#		CurrentFolderData[i] = AzimuthalPhaseConversion(CurrentFolderData[i],TEC2DVariableStrings[i])
		#endfor

		#Save all variables for folder[l] to Data.
		#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
		Data.append(CurrentFolderData)

	#===================##===================#
	#===================##===================#

	#	#Kinetics data readin - NOT CURRENTLY EMPLOYED IN ANY DIAGNOSTICS
	#	if True in [True]:
	#
	#		#Load data from TECPLOT_KIN file and unpack into 1D array.
	#		rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
	#		rawdata_kin.append(rawdata)
	#
	#		#Read through all variables for each file and stop when list ends.
	#		KinVariableStrings,KinHeaderEndMarker = ['T (S)'],'ZONE'
	#		for i in range(2,nn_2D):
	#			if KinHeaderEndMarker in str(rawdata_kin[l][i]):
	#				I = int(filter(lambda x: x.isdigit(), rawdata_kin[l][i].split(',')[0]))
	#				break
	#			else: KinVariableStrings.append(str(rawdata_kin[l][i][:-2].strip(' \t\n\r\"')))
	#			#endif
	#		#endfor
	#		Header_KIN.append(KinVariableStrings)												# SJD USE THIS HEADER !!!
	#		numvariables_kin,header_kin = len(KinVariableStrings),len(KinVariableStrings)+2		# REMOVE THIS
	#		header_kinlist.append(header_kin)													# REMOVE THIS
	#
	#		#Seperate total 1D data array into sets of data for each variable.
	#		CurrentFolderData = ReadTEC2D(rawdata_kin[l],header_kin,numvariables_kin, Zmesh=I,Dimension='1D')
	#
	#		#Convert data from CGS (HPEM DEFAULT) to user requested unit system						!!! Not tested
	#		for i in range(0,len(KinVariableStrings)):
	#			CurrentFolderData[i] = VariableUnitConversion(CurrentFolderData[i],KinVariableStrings[i], Units, AtomicSpecies)
	#			CurrentFolderData[i] = AzimuthalPhaseConversion(CurrentFolderData[i],KinVariableStrings[i])
	#		#endfor
	#
	#		#Save all variables for folder[l] to Data.
	#		#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
	#		DataKin.append(CurrentFolderData)
	#	#endif


	#===================##===================#
	#===================##===================#

		#Movie1 phase-resolved data readin
		if True in [savefig_sheathdynamics,savefig_phaseresolve1D,savefig_phaseresolve2D,savefig_PROES]:

			#Load data from movie_icp file and unpack into 1D array.
			try: rawdata, nn_movie1 = extract_raw_data(Dir, movie1[l].split('/')[-1], l)
			except: rawdata, nn_movie1 = extract_raw_data(Dir, 'movie1.pdt', l)
			rawdata_phasemovie.append(rawdata)

			#Read through all variables for each file and stop when list ends.
			#movie_icp has geometry at top, therefore len(header) != len(variables).
			#Only the first encountered geometry is used to define variable zone.
			VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ZONE'
			PhaseVariableStrings,numvar = list(),0
			for i in range(2,nn_movie1):
				if HeaderEndMarker in str(rawdata[i]):
					header_phase = i+1		# +1 to skip empty row below "ITER".
					break
				if VariableEndMarker in str(rawdata[i]) and numvar == 0:
					numvar = (i-3)			# -3 to ignore "R,Z" and remove overflow.
				if len(rawdata[i]) > 1 and numvar == 0:
					PhaseVariableStrings.append(str(rawdata_phasemovie[l][i][:-2].strip(' \t\n\r\"')))
				#endif
			#endfor
			Header_movie1.append(PhaseVariableStrings)			# SJD USE THIS HEADER !!!
			header_phasemovie.append(header_phase)				# REMOVE THIS

			#Rough method of obtaining the movie_icp iter locations for data extraction.
			Iterloc = list()
			Moviephaselist.append(list())
			for i in range(0,len(rawdata)):
				if "CYCL=" in rawdata[i]:
					Iterloc.append(i+1)

					IterStart=rawdata[i].find('CYCL')
					Moviephaselist[l].append(rawdata[i][IterStart:IterStart+9])
				#endif
			#endfor
		#endif


	#===================##===================#
	#===================##===================#

		if True in [savefig_movieicp2D,savefig_timeaxis1D,savefig_movieicp1D,savefig_convergence]:

			#Load data from movie_icp file and unpack into 1D array.
			try: rawdata, nn_itermovie = extract_raw_data(Dir, movieicp[l].split('/')[-1], l)
			except: rawdata, nn_itermovie = extract_raw_data(Dir, 'movie_icp.pdt', l)
			rawdata_itermovie.append(rawdata)

			#Read through all variables for each file and stop when list ends.
			#movie_icp has geometry at top, therefore len(header) != len(variables).
			#Only the first encountered geometry is used to define variable zone.
			VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ITER'
			MovIcpVariableStrings,numvar = list(),0
			for i in range(2,nn_itermovie):
				if HeaderEndMarker in str(rawdata[i]):
					header_iter = i+1		# +1 to skip empty row below "ITER".
					break
				if VariableEndMarker in str(rawdata[i]) and numvar == 0:
					numvar = (i-3)			# -3 to ignore "R,Z" and remove overflow.
				if len(rawdata[i]) > 1 and numvar == 0:
					MovIcpVariableStrings.append(str(rawdata_itermovie[l][i][:-2].strip(' \t\n\r\"')))
				#endif
			#endfor
			Header_movieicp.append(MovIcpVariableStrings)			# SJD USE THIS HEADER !!!
			header_itermovie.append(header_iter)					# REMOVE THIS

			#Rough method of obtaining the movie_icp iter locations for data extraction.
			Iterloc = list()
			MovieIterlist.append(list())
			for i in range(0,len(rawdata)):
				if "ITER=" in rawdata[i]:
					Iterloc.append(i+1)

					IterStart=rawdata[i].find('ITER')
					MovieIterlist[l].append(rawdata[i][IterStart:IterStart+9])
				#endif
			#endfor

			#Initiate movie_icp.pdt 4D data array
			if l == 0:
				A = numfolders
				B = len(Iterloc)
				C = len(MovIcpVariableStrings)
				D = R_mesh[l]*Z_mesh[l]

				#IterMovieData is saved in form: [folder,iteration,variable,datapoint(R,Z)]
				IterMovieData = np.zeros([A,B,C,D], dtype=float)
			#endif

			#Cycle through all iterations for current datafile, appending per cycle.
			for i in range(0,len(Iterloc)):

				#R,Z arrays are saved only for first "Cycle", apply +2 variable index offset to ignore
				#CurrentFolderData is saved in form: [iteration,variable,datapoint(R,Z)]
				#CurrentIterData is saved in form: [variable,datapoint(R,Z)]
				if i == 0:
					CurrentIterData = read_TEC2D(rawdata, Iterloc[i], numvar + 2, R_mesh[l], Z_mesh[l], offset=2)
	#				print(i, shape(CurrentIterData))
	#				print(shape(CurrentFolderData))

					#Convert data from CGS (HPEM DEFAULT) to user requested unit system
					for j in range(0,len(MovIcpVariableStrings)):
						CurrentIterData[j] = variable_unit_conversion(CurrentIterData[j], MovIcpVariableStrings[j], Units, AtomicSpecies)
	#					CurrentIterData[j] = AzimuthalPhaseConversion(CurrentIterData[j],MovIcpVariableStrings[j])
					#endfor

					#Check if any variables exist in header that do not exist in data
					NumLoadedVariables = shape(CurrentIterData)[0]
					for j in range(1,NumLoadedVariables):
						ExpectedDataLength = R_mesh[l]*Z_mesh[l]
						ReadinDataLength = size(CurrentIterData[:][j])

						if ReadinDataLength != ExpectedDataLength:
							print('')
							print('WARNING:',MovIcpVariableStrings[j],'ARRAY SIZE IS INCORRECT')
							print('ZEROING',MovIcpVariableStrings[j],'DATA AND CONTINUING READIN')
							print('')
							CurrentIterData[j] = np.zeros([ExpectedDataLength])
						#endif
					#endfor
					CurrentFolderData = np.array([CurrentIterData[0:numvar]])

				#Later cycles do not save R,Z arrays so no variable index offset is required.
				#CurrentFolderData is saved in form: [iteration,variable,datapoint(R,Z)]
				#CurrentIterData is saved in form: [variable,datapoint(R,Z)]
				elif i > 0:
					CurrentIterData = read_TEC2D(rawdata, Iterloc[i], numvar, R_mesh[l], Z_mesh[l])
	#				print(i, shape(CurrentIterData))
	#				print(shape(CurrentFolderData))

					#Convert data from CGS (HPEM DEFAULT) to user requested unit system
					for j in range(0,len(MovIcpVariableStrings)):
						CurrentIterData[j] = variable_unit_conversion(CurrentIterData[j], MovIcpVariableStrings[j], Units, AtomicSpecies)
	#					CurrentIterData[j] = AzimuthalPhaseConversion(CurrentIterData[j],MovIcpVariableStrings[j])
					#endfor

					#Check if any variables exist in header that do not exist in data
					NumLoadedVariables = shape(CurrentIterData)[0]
					for j in range(1,NumLoadedVariables):
						ExpectedDataLength = R_mesh[l]*Z_mesh[l]
						ReadinDataLength = size(CurrentIterData[:][j])

						if ReadinDataLength != ExpectedDataLength:
							CurrentIterData[j] = np.zeros([ExpectedDataLength])
						#endif
					#endfor
					CurrentFolderData = np.concatenate((CurrentFolderData, np.array([CurrentIterData])), axis=0)
				#endif
			#endfor

			#Save all variables for folder[l] to IterMovieData.
			MaxIter = len(CurrentFolderData)						#Maximum Iter saved is always in order
			IterMovieData[l,0:MaxIter] = CurrentFolderData
		#endif


	#===================##===================#
	#===================##===================#

		#IEDF/NEDF file readin.
		if True in [savefig_IEDFangular,savefig_IEDFtrends]:

			#Define arguments and autorun conv_prof.exe if possible.
			#### THIS IS HACKY, WON'T ALWAYS WORK, ARGS LIST NEEDS AUTOMATING ####
			IEDFVarArgs = ['1','1','1','1','1'] 						#Works for 2 species 1 surface.
			ExtraArgs = ['1','1','1','1','1','1','1','1','1','1']#[]	#Hack For Additional Species
			Args = ['pcmc.prof','title','1','1','1'] + IEDFVarArgs + ExtraArgs + ['0','0']
			DirAdditions = ['iprofile_tec2d.pdt','nprofile_tec2d.pdt','iprofile_tec1d.pdt', 'nprofile_tec1d.pdt','iprofile_zones_tec1d.pdt','nprofile_zones_tec1d.pdt']
			#try: AutoConvProf('./conv_prof.exe',Args,DirAdditions)
			#except: print('ConvProf Failure:'+Dirlist[l])
			auto_conv_prof('./conv_prof.exe', Dirlist[l],
			               IMATSTATS, IPCMCSPEC,
			               Args, DirAdditions)

			#Load data from IEDFprofile file and unpack into 1D array.
			try: rawdata, nn_IEDF = extract_raw_data(Dir, iprofiletec2d[l].split('/')[-1], l)
			except: rawdata, nn_IEDF = extract_raw_data(Dir, 'iprofile_tec2d.pdt', l)
			rawdata_IEDF.append(rawdata)

			#Read through all variables for each file and stop when list ends.
			IEDFVariableStrings,HeaderEndMarker = ['Theta [deg]','Energy [eV]'],'ZONE'
			for i in range(2,nn_IEDF):
				#Grab EDFangle(I),EDFbins(J) values from the ZONE line, these outline the datablock size.
				if HeaderEndMarker in str(rawdata_IEDF[l][i]):
					I = list(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[0]))	#discrete digits
					I = int( ''.join(I) ); EDFangle = I					#Number of EDF angle bins [Integer]
					J = list(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[1]))	#discrete digits
					J = int( ''.join(J) ); EDFbins = J					#Number of EDF energy bins [Integer]
					break
				else: IEDFVariableStrings.append(str(rawdata_IEDF[l][i][:-2].strip(' \t\n\r\"')))
				#endif
			#endfor
			Header_IEDF.append(IEDFVariableStrings)
			numvariables_IEDF,header_IEDF = len(IEDFVariableStrings),len(IEDFVariableStrings)+2
			header_IEDFlist.append(header_IEDF)

			#Seperate total 1D data array into sets of data for each variable.
			#Data is stored in 2D array of shape: [EDFangle,EDFbins] or [I,J]
			CurrentFolderData = read_TEC2D(rawdata_IEDF[l], header_IEDF, numvariables_IEDF, Rmesh=I, Zmesh=J, offset=0)

			#Convert data from CGS (HPEM DEFAULT) to user requested unit system							!!! Not tested
			for i in range(0,len(IEDFVariableStrings)):
				CurrentFolderData[i] = variable_unit_conversion(CurrentFolderData[i], IEDFVariableStrings[i], Units, AtomicSpecies)
	#			CurrentFolderData[i] = AzimuthalPhaseConversion(CurrentFolderData[i],IEDFVariableStrings[i])
			#endfor

			#Save all variables for folder[l] to Data.
			#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
			DataIEDF.append(CurrentFolderData)
		#endif


	#===================##===================#
	#===================##===================#

		#EEDF data readin.											!!! UNDER DEVELOPMENT
		if True in [savefig_EEDF]:

			#Load data from MCS.PDT file and unpack into 1D array.
			try: rawdata, nn_mcs = extract_raw_data(Dir, boltztec[l].split('/')[-1], l)
			except: rawdata, nn_mcs = extract_raw_data(Dir, 'boltz_tec.pdt', l)
			rawdata_mcs.append(rawdata)

			#Unpack each row of data points into single array of floats.
			#Removing 'spacing' between the floats and ignoring variables above data.
			Energy,Fe = list(),list()
			for i in range(3,len(rawdata_mcs[l])):
				if 'ZONE' in rawdata_mcs[l][i]:
					EEDF_TDlist.append( rawdata_mcs[l][i].split('"')[-2].strip(' ') )
					DataEEDF.append([Energy,Fe])
					Energy,Fe = list(),list()
				#endif
				try:
					Energy.append( float(rawdata_mcs[l][i].split()[0]) )
					Fe.append( float(rawdata_mcs[l][i].split()[1]) )
				except:
					NaN_Value = 1
				#endtry
			#endfor


			a,b = 0,5
			for i in range(a,b):
				plt.plot(DataEEDF[i][0],DataEEDF[i][1], lw=2)
			plt.legend(EEDF_TDlist[a:b])
			plt.xlabel('Energy [eV]')
			plt.ylabel('f(e) [eV-3/2]')
			plt.show()
		#endif														!!! UNDER DEVELOPMENT


	#===================##===================#
	#===================##===================#
	#===================##===================#

	#Create global list of all variable names and find shortest list.
	for l in range(0,numfolders):
		#Alphabetize the VariableStrings and keep global alphabetized list.
		tempvarlist = enumerate_variables(Variables, Header_TEC2D[l])[1]
		tempvarlist = sort(tempvarlist)
		numvars = len(tempvarlist)

		Globalvarlist.append(tempvarlist)
		Globalnumvars.append(numvars)
	#endfor

	#Find the folder with the fewest variables:
	val, idx = min((val, idx) for (idx, val) in enumerate(Globalnumvars))
	MinSharedVariables = Globalvarlist[idx]

	#===================##===================#
	#===================##===================#


	#Empty and delete any non-global data lists.
	tempdata,tempdata2 = list(),list()
	data_array,templineout = list(),list()
	Energy,Fe,rawdata_mcs = list(),list(),list()
	HomeDir,DirContents = list(),list()
	#del RADIUS,RADIUST,HEIGHT,HEIGHTT,DEPTH,ISYM,IXZ
	#del data_array,tempdata,tempdata2,templineout
	#del Energy,Fe,rawdata_mcs
	#del HomeDir,DirContents


	#Alert user that readin process has ended and continue with selected diagnostics.
	if any([savefig_tecplot2D, savefig_phaseresolve2D, savefig_movieicp2D, savefig_convergence, savefig_monoprofiles, savefig_multiprofiles, savefig_compareprofiles, savefig_timeaxis1D, savefig_movieicp1D, savefig_sheathdynamics, savefig_phaseresolve1D, savefig_PROES, savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust, savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]) == True:
		print( '----------------------------------------')
		print( 'Data Readin Complete, Starting Analysis:')
		print( '----------------------------------------')
	else:
		print( '------------------')
		print( 'Analysis Complete.')
		print( '------------------')
	#endif


	#=====================================================================#
	#=====================================================================#


























	#====================================================================#
					  #COMMONLY USED PLOTTING FUNCTIONS#
	#====================================================================#

	#=========================#
	#=========================#

	def tecplot_cmap(NumLevels=256):
	#	Creates a colourmap closely approximating the Tecplot "modern" map
	#
	#		Python tutor:	https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html
	#		RGB Hex Codes:	https://www.rapidtables.com/web/color/RGB_Color.html
	#		Colour Picker:	sudo apt-get install gpick
	#
	# 	Colour dictionary details the vertices of the linear colour scales for rgb
	#		x1		::	Fraction of cbar range (0=min, 1=max)
	#		yleft	::  Percent of colour at start of X1 range
	#		yright	::	Percent of colour at end of X1 range
	#
	#	Colour Order (low to high):
	#		FFFFFF 100.0% red, 100.0% green, 100.0% blue
	#		FF69B4 100.0% red, 41.01% green, 70.31% blue								!!! NOT SWATCHED
	#		43006F 27.17% red, 00.00% green, 43.35% blue
	#		5BC4CB 35.54% red, 76.56% green, 79.29% blue
	#		00C952 00.00% red, 78.51% green, 32.03% blue
	#		CDF100 80.07% red, 94.14% green, 00.00% blue
	#		AC6313	67.18% red, 38.67% green, 07.42% blue
	#		F70000	96.86% red, 00.00% green, 00.00% blue
	#
		from matplotlib.colors import ListedColormap, LinearSegmentedColormap
	#
	#						Cbar%		FadeFrom	FadeTo
							#x1			#yleft		#yright
		cdict0 = {'red':  [[0.00,		0.00,		0.26],
						   [0.17,		0.26,		0.36],
						   [0.33,		0.36,		0.00],
						   [0.50,		0.00,		0.80],
						   [0.66,		0.80,		0.67],
						   [0.83,		0.67,		0.97],
						   [1.00,		0.97,		1.00]],

				 'green': [[0.00,		0.00,		0.00],
						   [0.17,		0.00,		0.77],
						   [0.33,		0.77,		0.79],
						   [0.50,		0.79,		0.94],
						   [0.66,		0.94,		0.38],
						   [0.83,		0.38,		0.00],
						   [1.00,		0.00,		0.00]],

				 'blue':  [[0.00,		0.00,		0.43],
						   [0.17,		0.43,		0.79],
						   [0.33,		0.79,		0.32],
						   [0.50,		0.32,		0.00],
						   [0.66,		0.00,		0.07],
						   [0.83,		0.07,		0.00],
						   [1.00,		0.00,		0.00]]}

							#x1			#yleft		#yright
		cdict1 = {'red':   [[0.00,		1.00,		1.00],
						   [0.14,		1.00,		0.26],
						   [0.28,		0.26,		0.36],
						   [0.42,		0.36,		0.00],
						   [0.56,		0.00,		0.80],
						   [0.70,		0.80,		0.67],
						   [0.84,		0.67,		0.97],
						   [1.00,		0.97,		1.00]],

				 'green': [[0.00,		1.00,		1.00],
						   [0.14,		0.41,		0.00],
						   [0.28,		0.00,		0.77],
						   [0.42,		0.77,		0.79],
						   [0.56,		0.79,		0.94],
						   [0.70,		0.94,		0.38],
						   [0.84,		0.38,		0.00],
						   [1.00,		0.00,		0.00]],

				 'blue':  [[0.00,		1.00,		1.00],
						   [0.14,		0.70,		0.43],
						   [0.28,		0.43,		0.79],
						   [0.42,		0.79,		0.32],
						   [0.56,		0.32,		0.00],
						   [0.70,		0.00,		0.07],
						   [0.84,		0.07,		0.00],
						   [1.00,		0.00,		0.00]]}

		cdict2 = {'red': [[0, 1.0, 1.0],
							[1 / 7, 1.0, 35.0 / 255],
							[2 / 7, 0.0, 35.0 / 255],
							[3 / 7, 0.0, 35.0 / 255],
							[4 / 7, 0.0, 100.0 / 255],
							[5 / 7, 1.0, 85.0 / 255],
							[6 / 7, 217.0 / 255, 100.0 / 255],
							[7 / 7, 1.0, 1.0]],
						'green': [[0, 35.0 / 255, 1.0],
							[1 / 7, 0.0, 35.0 / 255],
							[2 / 7, 0.0, 100.0 / 255],
							[3 / 7, 1.0, 100.0 / 255],
							[4 / 7, 1.0, 100.0 / 255],
							[5 / 7, 1.0, 46.0 / 255],
							[6 / 7, 117.0 / 255, 35.0 / 255],
							[7 / 7, 0.0, 0.0]],
						'blue': [[0, 100.0 / 255, 1.0],
							[1 / 7, 1.0, 100.0 / 255],
							[2 / 7, 1.0, 100.0 / 255],
							[3 / 7, 1.0, 35.0 / 255],
							[4 / 7, 0.0, 35.0 / 255],
							[5 / 7, 0.0, 10.0 / 255],
							[6 / 7, 26.0 / 255, 35.0 / 255],
							[7 / 7, 0.0, 0.0]]}

		cdict = cdict2
		cmap = LinearSegmentedColormap('TecplotModern', segmentdata=cdict, N=NumLevels)
		rgba = cmap(np.linspace(0, 1, 256))

		return(cmap,cdict)
	#enddef

	#=========================#
	#=========================#

	#Load any custom colourmaps

	#Load Tecplot Modern colourmap
	if image_cmap == 'tecmodern':
		tecplotcmap,tecplotcdict = tecplot_cmap(NumLevels=256)
		image_cmap = tecplotcmap
	#endif

	#Load IDL colourmap (std_gamma_II.txt in modules dir)
	Filename = 'Modules/std_gamma_II.txt'
	try:
		map = np.loadtxt(Filename, delimiter=',')
		IDL_Gamma_II = col.ListedColormap(map.T, name='IDL_Gamma_II')
	except:
		print('Warning: gamma_II colourmap not found')
	#endtry

	#=========================#
	#=========================#


	def Matplotlib_GlobalOptions():
	#Takes global inputs from switchboard, returns nothing
	#Alters global image options, run before any diagnostics
	#Attempts to revert matplotlib changes made in 2.0 onwards.
	#See: https://matplotlib.org/users/dflt_style_changes.html

	#	mpl.style.use('classic')								#Resets to classic 1.x.x format

		#Image options
		mpl.rcParams['figure.figsize'] = [10.0,10.0]			#Sets default figure size
		mpl.rcParams['figure.dpi'] = 100						#Sets viewing dpi
		mpl.rcParams['savefig.dpi'] = 100						#Sets saved dpi
		mpl.rcParams['image.interpolation'] = image_interp		#Applies bilinear image 'smoothing'
		mpl.rcParams['image.resample'] = True					#Resamples data before colourmapping
		mpl.rcParams['image.cmap'] = image_cmap 				#Define global colourmap

		#Axis options
		mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'	#View limits coencide with axis ticks
		mpl.rcParams['axes.xmargin'] = 0						#Set default x-axis padding
		mpl.rcParams['axes.ymargin'] = 0						#Set default y-axis padding
		mpl.rcParams['errorbar.capsize'] = 3					#Set error bar end cap width
		mpl.rcParams['font.size'] = 20							#Set global fontsize
		mpl.rcParams['legend.fontsize'] = 'large'				#Set legend fontsize
		mpl.rcParams['figure.titlesize'] = 'medium'				#Set title fontsize

		#Line and Colour options
	#	from cycler import cycler								#See below
	#	mpl.rcParams['axes.prop_cycle']=cycler(color='bgrcmyk')	#Set default colour names
		mpl.rcParams['lines.linewidth'] = 1.0					#Set Default linewidth

		#Maths and Font options
		mpl.rcParams['mathtext.fontset'] = 'cm'					#Sets 'Latex-like' maths font
		mpl.rcParams['mathtext.rm'] = 'serif'					#Sets default string font

		return()
	#enddef
	Matplotlib_GlobalOptions()	#MUST BE RUN BEFORE ANY DIAGNOSTICS!!!!


	#=========================#
	#=========================#


	def plot_linearmap(cdict):
	#	Shows linear rgb colour fractions for cmap colour dictionary
	#
	#	USAGE:
	#		cmap,cdict = tecplot_cmap()
	#		plot_linearmap(cdict)

		from matplotlib.colors import ListedColormap, LinearSegmentedColormap

		newcmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
		rgba = newcmp(np.linspace(0, 1, 256))

		fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)

		col = ['r', 'g', 'b']
		for xx in [0.25, 0.5, 0.75]:
			ax.axvline(xx, color='0.7', linestyle='--')
		#endfor
		for i in range(3):
			ax.plot(np.arange(256)/256, rgba[:, i], color=col[i])
		#endfor

		ax.set_xlabel('index')
		ax.set_ylabel('RGB')
		plt.show()
	#enddef


	#=========================#
	#=========================#


	def ImageExtractor2D(Data,Variable=[],Rmesh=0,Zmesh=0):
	#Returns a 2D array of inputted data with size [R_mesh] x [Z_mesh]
	#Can optionally perform variable unit conversion if required.
	#Image = ImageExtractor2D(Data,Variable=[]):

		#If no mesh sizes supplied, collect sizes for current global folder.
		if Rmesh == 0 or Zmesh == 0:
			Rmesh,Zmesh = int(R_mesh[l]),int(Z_mesh[l])
		#endif

		#Create empty 2D image of required size.
		numrows = int(len(Data)/Rmesh)
		Image = np.zeros([int(numrows),int(Rmesh)])

		#Reshape data into 2D array for further processing.
		for j in range(0,numrows):
			for i in range(0,Rmesh):
				Start = Rmesh*j
				Row = Zmesh-1-j
				Image[Row,i] = Data[Start+i]
			#endfor
		#endfor

		#Convert units if required.								!!! RM SJD, CONVERSION PERFORMED AT READ-IN
	#	Image = VariableUnitConversion(Image,Variable, Units, AtomicSpecies)

		#Convert Azimuthal phase if required.					!!! RM SJD, MOVE THIS CONVERSION TO READ-IN
		Image = azimuthal_phase_conversion(Image, Variable)

		return(Image)
	#enddef


	#=========================#
	#=========================#


	def SymmetryConverter(Image,Radial=False):
	#Takes 2D image array and produces symmetric image about central axis.
	#Returns symmetric 2D image array, allows radially negative values.
	#Image = SymmetryConverter(Image,Radial=False)

		#Create new image by reversing and adding itself on the LHS.
		if image_plotsymmetry == True and int(ISYMlist[l]) == 1:
			SymImage = np.zeros([len(Image),2*len(Image[0])])
			if Radial == False:
				for m in range(0,len(Image)):
					SymImage[m] = np.concatenate([Image[m][::-1],Image[m]])
				#endfor
			elif Radial == True:
				for m in range(0,len(Image)):
					SymImage[m] = np.concatenate([-Image[m][::-1],Image[m]])
				#endfor
			#endif
			Image = SymImage
		#endif

		return(Image)
	#enddef


	#=========================#
	#=========================#


	def ScaleArray(Array,ScaleFactors):
	#Scales a 2D array by the requested X and Y scale factors, can scale asymmetrically.
	#Takes a rectilinear 2D array of floats and returns a scaled rectilinear 2D array.
	#Employs a linear interpolation scheme when mapping new datapoints.
	#ScaledArray = ScaleArray(2DArray,ScaleFactors=[3,3])

		#Naughty module import
		from scipy.ndimage.interpolation import map_coordinates

		#Define old and new array scales based upon supplied factors
		OldScale = [len(Array),len(Array[0])]
		NewScale = [int(OldScale[0]*ScaleFactors[0]),int(OldScale[1]*ScaleFactors[1])]

		#Create new 2D array with required NewScale
		Array,NewDims = np.asarray(Array),list()
		for OriginalLength, NewLength in zip(Array.shape, (NewScale[0],NewScale[1])):
			NewDims.append(np.linspace(0, OriginalLength-1, NewLength))
		#Map old data onto new array employing a linear interpolation
		Coords = np.meshgrid(*NewDims, indexing='ij')
		ScaledArray = map_coordinates(Array, Coords)

		return(ScaledArray)
	#enddef


	#=========================#
	#=========================#

	def figure(aspectratio=[],subplots=1,shareX=False):
	#Create figure of desired size and with variable axes.
	#Returns figure and axes seperately.
	#fig,ax = figure(image_aspectratio,1,shareX=False)

		if len(aspectratio) == 2:
			fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]),sharex=shareX)
		else:
			fig, ax = plt.subplots(subplots, figsize=(10,10), sharex=shareX)
		#endif
		return(fig,ax)
	#enddef

	#=========================#
	#=========================#

	def clearfigures(fig):

		# Clear the current axes.
		plt.cla()
		# Clear the current figure.
		plt.clf()
		# Closes all the figure windows.
		plt.close('all')
		plt.close(fig)
		gc.collect()

		return()
	#enddef

	#=========================#
	#=========================#


	def CropImage(ax,Extent=[],Apply=True,Rotate=True):
	#Crops 2D image taking account of image rotation options.
	#Takes image axis (assumes default axis), use figure()
	#Input Extent format: [ [Rmin,Rmax], [Zmin,Zmax] ] in cm.
	#Returns cropping limits in format: [ [R1,R2],[Z1,Z2] ] in cm.
	#CropImage(ax[0],Extent=[[R1,R2],[Z1,Z2]],Apply=True,Rotate=True),

		#Obtain default limits and rotate if needed, doesn't crash if no crop applied.
		#R1,R2 are radial limits of image, Z1,Z2 are axial limits. (non-rotated)
		R1,R2 = ax.get_xlim()[0],ax.get_xlim()[1]
		Z1,Z2 = ax.get_ylim()[0],ax.get_ylim()[1]
		if Rotate == True and image_rotate == True:
			R1,Z1, R2,Z2 = Z1,R1, Z2,R2
		#endif

		#Set requested cropping limits from function-call or from global.
		if len(Extent) == 2:
			radialcrop,axialcrop = Extent[0],Extent[1]
		else:
			radialcrop = image_radialcrop
			axialcrop = image_axialcrop
		#endif

		#Extract cropping dimentions from image_<input>.
		if len(radialcrop) == 1:
			R1,R2 = -(radialcrop[0]),radialcrop[0]
		elif len(radialcrop) == 2:
			R1,R2 = radialcrop[0],radialcrop[1]
		#endif
		if len(axialcrop) == 1:
			Z1,Z2 = 0,axialcrop[0]
		elif len(axialcrop) == 2:
			Z1,Z2 = axialcrop[0],axialcrop[1]
		#endif

		#Rotate cropping dimentions to match image rotation.
		if Rotate == True:
			if image_rotate == True:
				R1,Z1 = Z1,R1
				R2,Z2 = Z2,R2
			elif image_rotate == False:
				Z1,Z2 = Z1,Z2
			#endif
		#endif

		#Apply cropping dimensions to image.
		if Apply == True:
			ax.set_xlim(R1,R2)
			ax.set_ylim(Z1,Z2)
		#endif

		#Return cropped dimensions in SI units.
		return([[R1,R2],[Z1,Z2]])
	#enddef


	#=========================#
	#=========================#

	###
	#NB THIS NEEDS A COMPLETE OVERHAUL AND SPLIT INTO SIMPLER FUNCTIONS!
	###

	def CbarMinMax(ax,Image,Symmetry=image_plotsymmetry,PROES=False):
	#Determines min/max colourbar scale within the cropped frame of an image array.
	#Takes a 2D image, and returns the min/max value within the cropped region.
	#Assumes image symmetry for best results, otherwise R1 is set to zero.
	#Works for PROES images too, requires PROES='Axial' or 'Radial'.
	#[Minimum,Maximum] = CbarMinMax(Image,Symmetry=False,PROES=False)

		#Return user defined limits if specified.
		if len(image_cbarlimit) == 2:
			cropmin = image_cbarlimit[0]
			cropmax = image_cbarlimit[1]
			return([cropmin,cropmax])
		#endif

		#Extract min/max from cropped region if a region is supplied.
		if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:

			#Import global cell sizes with correct rotation.
			if image_rotate == False: dR,dZ = dr[l],dz[l]
			if image_rotate == True:  dR,dZ = dz[l],dr[l]
			#endif

			#Convert cropped SI region (CropExtent) to cell region (R1,R2,Z1,Z2).
			#CropExtent is extracted using CropImage function, which applies a rotation.
			CropExtent = CropImage(ax,Apply=False)		#CropExtent applies rotation internally

			R1 = int(CropExtent[0][0]/dR)
			R2 = int(CropExtent[0][1]/dR)
			Z1 = int(CropExtent[1][0]/dZ)
			Z2 = int(CropExtent[1][1]/dZ)
			#endif

			#Re-rotate back so that radial and axial crops align correctly.
			#R1,R2 are radial limits of image, Z1,Z2 are axial limits.
			if image_rotate == True:
				R1,Z1 = Z1,R1
				R2,Z2 = Z2,R2
			#Cell Images have no negative R values, unlike SI (top left cell = 0,0)
			#To align properly such that R=0 on axis, add R_mesh to both cropping limits.
			if Symmetry == True:
				R1 = int(R_mesh[l] + R1)
				R2 = int(R_mesh[l] + R2)
			#If image has no symmetry applied, replace negative R1 with zero to align axis.
			if Symmetry == False and R1 < 0:
				R1 = 0
			#endif
			#Reverse axial ROI if required to avoid zero size images, does not affect plotting.
			#HPEM axial zero in top left corner, array must be read with zero in 'bottom left'.
			if Z1 > Z2:
				Z1,Z2 = Z2,Z1
			#endif

			#Crop the cell region to the desired region
			if PROES == False:
				#Default Orientation, image orientated radially then axially Image[R][Z]
				if Symmetry == False:
					Image = Image[Z1:Z2]
					Image = np.asarray(Image).transpose()
					Image = Image[R1:R2]
					Image = np.asarray(Image).transpose()
				#90 Deg Rotated Orientation, image orientated axially then radially Image[Z][R]
				elif Symmetry == True:
					Image = Image[Z1:Z2]
					Image = np.asarray(Image).transpose()
					Image = Image[R1:R2]
					Image = np.asarray(Image).transpose()
				#endif
			#Crop the cell region for PROES images, they only require y-axis cropping (R,Z).
			elif PROES == 'Axial':
				Image = Image[::-1][Z1:Z2]			#Image[::-1] Accounts for reversed PhaseData order
			elif PROES == 'Radial':
				Image = Image[R1:R2]
			#endif
		#endif

		#=====#=====#

		#Flatten image and obtain min/max in region, defaults to full image if no cropping.
		#1D array of 2D image, concats each radial profile in axial order
		flatimage = [item for sublist in Image for item in sublist]
		flatimage = [float(x) for x in flatimage]

		#Filtered 1D array, removing any nan's or inf's
		try:
			flatimage = [x for x in flatimage if not (m.isinf(x) or m.isnan(x))]
		except:
			print( 'Image Filtering Warning - Cbar may be scaled incorrectly' )
		#endtry

		#Take min and max value from the filtered array to use as cbar limits
		if len(flatimage) > 0:
			cropmin, cropmax = min(flatimage), max(flatimage)
		else:
			print('Cbar Limits Unspecified; Setting to 0,1')
			cropmin, cropmax = 0,1
		#endif

		#Return cropped values in list [min,max], as required by colourbar.
		return([cropmin,cropmax])
	#enddef


	#=========================#
	#=========================#

	def ImageGeometry(fig,ax,image_plotsymmetry):

		#Plot mesh outline if requested.			!!! RM SJD, SHOULD EXPLICITLY PASS FOLDER [l]
		if image_plotmesh == True:
			PlotGeometry(ax,MeshCoordinates[l],MeshConnections[l],image_plotsymmetry)
		elif image_plotmesh == 'PRCCP':
			ManualPRCCPMesh(ax)
		elif image_plotmesh == 'PRCCPM':
			ManualPRCCPMMesh(ax)
		elif image_plotmesh == 'GEC':
			ManualGECMesh(ax,image_plotsymmetry)
		elif image_plotmesh == 'EVgeny':
			ManualEVgenyMesh(ax)
		elif image_plotmesh == 'HyperionI':
			ManualHyperionIMesh(ax)
		elif image_plotmesh == 'HyperionII':
			ManualHyperionIIMesh(ax)
		elif image_plotmesh == 'ASTRONmk1':
			ManualASTRONmk1Mesh(ax)
		#endif

		return()
	#enddef

	#=========================#
	#=========================#

	def ImageOptions(fig,ax,Xlabel='',Ylabel='',Title='',Legend=[],Crop=True,Rotate=True):
	#Applies plt.options to current figure based on user input.
	#Returns nothing, open figure is required, use figure().
	#For best results call immediately before saving/displaying figure.
	#ImageOptions(plt.gca(),Xlabel,Ylabel,Title,Legend,Crop=False)

		#Apply user overrides to plots.
		if len(titleoverride) > 0:
			Title = titleoverride
		if len(legendoverride) > 0:
			Legend = legendoverride
		if len(xlabeloverride) > 0:
			Xlabel = xlabeloverride[0]
		if len(ylabeloverride) > 0:
			Ylabel = ylabeloverride[0]
		#endif

		#Set title and legend if one is supplied.
		if len(Title) > 0:
			ax.set_title(Title, fontsize=14, y=1.03)
		if len(Legend) > 0:
			ax.legend(Legend, loc=image_legendloc, fontsize=16, frameon=False)
		#endif

		#Set axis labels
		if image_axislabels == True:
			ax.set_xlabel(Xlabel, fontsize=24)
			ax.set_ylabel(Ylabel, fontsize=24)
		elif image_axislabels == False:
			ax.set_xlabel('')
			ax.set_ylabel('')
		#endif

		#Set axis ticks
		if image_axisticks == True:
			ax.tick_params(axis='x', labelsize=18)
			ax.tick_params(axis='y', labelsize=18)
		elif image_axisticks == False:
			ax.tick_params(axis='x', labelbottom=False, labelleft=False, length=0)
			ax.tick_params(axis='y', labelbottom=False, labelleft=False, length=0)
		#endif

		#Force scientific notation for all axes, accounting for non-scalar x-axes.
		#	RM: THIS SEEMS TO ALSO BE DEPRECIATED - NEED TO REPLACE WITH NEW Matplotlib COMMAND
	#	try: ax.xaxis.get_major_locator().set_params(style='sci',scilimits=(-2,3),axis='both')
	#	except: Axes_Contain_Strings = True
		#
	#	try: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='both')		#Matplotlib v2.x.x TICKFORMAT
	#	except: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='y')		#Matplotlib v2.x.x TICKFORMAT
		#endtry

		#Set grid, default is off.
		if image_plotgrid == True: ax.grid(True)
		#endif

		#Plot mesh outline if requested.
		if image_plotmesh and Crop == True:
			PlotGeometry(ax,MeshCoordinates,MeshConnections,image_plotsymmetry)
		elif image_plotmesh == 'PRCCP' and Crop == True:
			ManualPRCCPMesh(ax)
		elif image_plotmesh == 'PRCCPM' and Crop == True:
			ManualPRCCPMMesh(ax)
		elif image_plotmesh == 'GEC' and Crop == True:
			ManualGECMesh(ax,image_plotsymmetry)
		elif image_plotmesh == 'EVgeny' and Crop == True:
			ManualEVgenyMesh(ax)
		elif image_plotmesh == 'HyperionI' and Crop == True:
			ManualHyperionIMesh(ax)
		elif image_plotmesh == 'HyperionII' and Crop == True:
			ManualHyperionIIMesh(ax)
		elif image_plotmesh == 'ASTRONmk1' and Crop == True:
			ManualASTRONmk1Mesh(ax)
		#endif

		#Crop image dimensions, use provided dimensions or default if not provided		- REMOVE THIS FROM HERE...
		if isinstance(Crop, (list, np.ndarray) ) == True:
			CropImage(ax,Extent=Crop,Rotate=Rotate)
		elif any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
			if Crop == True:
				CropImage(ax,Rotate=Rotate)
			#endif
		#endif

		#Arrange figure such that labels, legends and titles fit within frame.
		fig.tight_layout()

		return()
	#enddef


	#=========================#
	#=========================#


	def Colourbar(ax=plt.gca(),Label='',cbarbins=5,Lim=[]):
	#Creates and plots a colourbar with given label and binsize.
	#Takes image axis, label string, number of ticks and limits
	#Allows pre-defined colourbar limits in form [min,max].
	#Returns cbar axis if further changes are required.
	#cbar = Colourbar(ax[0],'Label',5,Lim=[0,1])

		#Set default font and spacing options and modify if required
		Rotation,Labelpad = 270,30
		LabelFontSize,TickFontsize = 24,18
		if '\n' in Label: Labelpad += 25		#Pad label for multi-line names

		#Create and define colourbar axis
		try:
			divider = make_axes_locatable(ax)
			cax = divider.append_axes("right", size="2%", pad=0.1)
			cbar = plt.colorbar(im, cax=cax)
		except:
			print('Warning: Colourbar failed to plot, check Colourbar() function')
			return()
		#endtry

		#Set number of ticks, label location and define scientific notation.
		if image_cbarticks == True:
			cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad,fontsize=LabelFontSize)
			cbar.locator = ticker.MaxNLocator(nbins=cbarbins)
			cbar.ax.tick_params(labelsize=TickFontsize)
			cbar.formatter.set_powerlimits((-2,3))
			cbar.update_ticks()

		elif image_cbarticks == False:
			cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad-30,fontsize=LabelFontSize)
			cmin, cmax = im.get_clim()
			cbar.set_ticks([cmin, cmax])
			cbar.set_ticklabels(['min', 'max'])
		#endif

		#Apply colourbar limits if specified.  (lim=[min,max])
		if len(Lim) == 2: im.set_clim(vmin=Lim[0], vmax=Lim[1])

		return(cbar)
	#enddef


	#=========================#
	#=========================#


	def InvisibleColourbar(ax='NaN'):
	#Creates an invisible colourbar to align subplots without colourbars.
	#Takes image axis, returns colourbar axis if further edits are required
	#cax = InvisibleColourbar(ax[0])

		#If no axis is supplied, use current open axis
		if ax == 'NaN': ax = plt.gca()

		#Create colourbar axis, ideally should 'find' values of existing cbar!
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="2%", pad=0.1)

		#Set new cax to zero size and remove ticks.
		try: cax.set_facecolor('none')				#matplotlib v2.x.x method
		except: cax.set_axis_bgcolor('none')		#matplotlib v1.x.x method
		for axis in ['top','bottom','left','right']:
			cax.spines[axis].set_linewidth(0)
		cax.set_xticks([])
		cax.set_yticks([])

		return(cax)
	#enddef



	#=========================#
	#=========================#


	def GenerateAxis(Orientation,Isym=ISYMlist[l],PhaseFrames=range(0,IMOVIE_FRAMES[l])):
	#Generates a 1D SI [cm] axis for plotting, includes radial symmetry.
	#Takes orientation, symmetry and phasecycle options.
	#Returns 1D array in units of [cm] or [omega*t/2pi].
	#Raxis=GenerateAxis('Radial',Isym=ISYMlist[l])

		#Create axis list and extract the number of phaseframes
		PHASEEResolution = len(PhaseFrames)
		axis = list()
		if Orientation == 'Radial':
			if int(Isym) == 1:
				for i in range(-int(R_mesh[l]),int(R_mesh[l])):
					axis.append(i*dr[l])
			#endfor
			elif int(Isym) != 1:
				for i in range(0,int(R_mesh[l])):
					axis.append(i*dr[l])
				#endfor
			#endif
		elif Orientation == 'Axial':
			for i in range(0,int(Z_mesh[l])):
				axis.append(i*dz[l])
			#endfor
		elif Orientation == 'Phase':
			for i in range(0,int(phasecycles*PHASEEResolution)):
				axis.append(  (np.pi*(i*2)/PHASEEResolution)/(2*np.pi)  )
			#endfor
		#endif
		return(axis)
	#enddef


	#=========================#
	#=========================#


	def Normalise(profile,NormFactor=0):
	#Takes 1D or 2D array and returns array Normalised to maximum value.
	#If NormFactor is defined, array will be Normalised to this instead.
	#Returns Normalised image/profile and the max/min normalization factors.
	#NormProfile,Min,Max = Normalise(profile,NormFactor=0)

		NormalisedImage = list()

		#determine dimensionality of profile and select normaliztion method.
		if isinstance(profile[0], (list, np.ndarray) ) == True:

			#Obtain max and min normalization factors for 2D array.
			FlatImage = [item for sublist in profile for item in sublist]
			MaxNormFactor,MinNormFactor = max(FlatImage),min(FlatImage)

			#Fix for division by zero and other infinity related things...
			if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
			if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0
			#endif

			#Normalise 2D array to local maximum.
			if NormFactor == 0: NormFactor = MaxNormFactor
			for i in range(0,len(profile)):
				NormalisedImage.append( [x/NormFactor for x in profile[i]] )
			#endfor
			profile = NormalisedImage
			return(profile,MaxNormFactor,MinNormFactor)

		#Lowest dimention is still list.
		elif isinstance(profile, (list, np.ndarray) ) == True:

			#Obtain max and min normalization factors for 1D profile.
			MaxNormFactor,MinNormFactor = max(profile),min(profile)

			#Fix for division by zero and other infinity related things...
			if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
			if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0

			#Normalise 1D array to local maximum.
			if NormFactor == 0: NormFactor = MaxNormFactor
			for i in range(0,len(profile)):
				profile[i] = profile[i]/NormFactor
			#endfor
		#endif

		return(profile,MinNormFactor,MaxNormFactor)
	#enddef


	#=========================#
	#=========================#


	def DataExtent(folder,aspectratio=image_aspectratio):
	#Takes current image datails and returns extent and rotated aspectratio
	#If mesh uses symmetry, will double radius extent centered on zero.
	#extent,aspectratio = DataExtent(l)

		#Obtain global variables for current folder.
		Isym = float(ISYMlist[folder])
		radius,height = Radius[folder],Height[folder]
		#Rotated Image: [X,Y] = [Height,Radius]
		if image_rotate == True:
			aspectratio = aspectratio[::-1]
			if Isym == 1:
				extent= [0,height, -radius, radius]
			elif Isym == 0:
				extent= [0,height,  0,      radius]
			#endif

		#Default mesh orientation: [X,Y] = [Radius,Height]
		elif image_rotate == False:
			if Isym == 1: extent = [-radius,radius,  height,0]
			elif Isym == 0: extent=[0      ,radius,  height,0]
			#endif
		#endif
		return(extent,aspectratio)
	#enddef


	#=========================#
	#=========================#


	def ImagePlotter1D(axis,profile,aspectratio,fig=111,ax=111):
	#Create figure and plot a 1D graph with associated image plotting requirements.
	#Returns plotted axes and figure if new ones were created.
	#Else plots to existing figure and returns the image object.
	#ImagePlotter1D(Zaxis,Zprofile,image_aspectratio,fig,ax[0]):

		#Generate new figure if required. {kinda hacky...}
		if fig == 111 and ax == 111:
			fig,ax = figure(aspectratio)
		elif fig == 111:
			fig = figure(aspectratio)
		#endif

		#Apply any required numerical changes to the profile.
		if image_logplot == True:
			profile = np.log10(profile)
		if image_normalise == True:
			profile = Normalise(profile)[0]
		#endif

		#Plot profile and return.
		im = ax.plot(axis,profile, lw=2)

		try: return(fig,ax,im)
		except: return()
	#enddef


	#=========================#
	#=========================#


	def ImagePlotter2D(Image,extent,aspectratio=image_aspectratio,variable='N/A',fig=111,ax=111):
	#Create figure and plot a 2D image with associated image plotting requirements.
	#Returns plotted image, axes and figure after applying basic data restructuring.
	#fig,ax,im,Image = ImagePlotter2D(Image,extent,image_aspectratio,VariableStrings[l],fig,ax[0])

		#Generate new figure if required. {kinda hacky...}
		if fig == 111 and ax == 111:
			fig,ax = figure(aspectratio)
		elif fig == 111:
			fig = figure(aspectratio)
		#endif

		#Apply image axis-symmetry, with negative values, if required.
		RadialBool = is_radial_variable(variable)
		Image = SymmetryConverter(Image,RadialBool)

		#Rotate image if required
		if image_rotate == True:
			Image = np.asarray(Image)
			Image = Image.transpose().tolist()
		#endif

		#Apply any required numerical changes to the image.
		if image_logplot == True:
	#		Image = np.log10(Image)								# ln(0) NaN issues in materials
			Image = where(Image != 0, np.log10(Image), 0)		# Avoids ln(0) in images
		#endif

		if image_normalise == True:
			Image = Normalise(Image)[0]
		#endif

		#Plot image with colour fill and contour lines
		if image_plotcolourfill == True and image_plotcontours == True:
			im = ax.imshow(Image,extent=extent,origin='upper',aspect='equal')
	#		im = ax.contourf(Image,extent=extent,origin="upper",levels=999)	# Stretches Pixels Relative to Aspect_Ratio
			im2 = ax.contour(Image,extent=extent,origin="upper",levels=image_contourlvls,cmap='Greys',alpha=0.25)

		#Plot image with only colour fill
		elif image_plotcolourfill == True:
			im = ax.imshow(Image,extent=extent,origin='upper',aspect='equal')
	#		im = ax.contourf(Image,extent=extent,origin="upper",levels=999)	# Stretches Pixels Relative to Aspect_Ratio

		#Plot image with only contour lines
		elif image_plotcontours == True:
			im = ax.contour(Image,extent=extent,origin="upper",levels=image_contourlvls)
		#endif

		return(fig,ax,im,Image)
	#enddef


	#=========================#
	#=========================#


	def TrendPlotter(ax=plt.gca(),TrendArray=[],Xaxis=[],Marker='o-',NormFactor=0):
	#Creates a 1D image from an array of supplied points.
	#Image plotted onto existing axes, figure() should be used.
	#NormFactor = 0 will Normalise to maximum of given profile.
	#TrendPlotter(ax[0],TrendProfiles,Xaxis,'o-',0)

		#Normalise data to provided normalization factor if required.
		if image_normalise == True:
			TrendArray = Normalise(TrendArray,NormFactor)[0]
		#endif

		#Plot results against strings pulled from folder names for batch studies.
		Plot = ax.plot(range(0,numfolders),TrendArray, Marker, lw=2)
		if len(xaxisoverride) > 0:
			ax.set_xticks(np.arange(0,numfolders))
			ax.set_xticklabels(xaxisoverride)
		else:
			ax.set_xticks(np.arange(0,numfolders))
			ax.set_xticklabels(Xaxis)
		#endif

		return(Plot)
	#enddef


	#=========================#
	#=========================#


	def ExtractRadialProfile(Data,process,variable,Profile,Rmesh='NaN',Isym='NaN'):
	#Obtains a radial 1D profile at a requested axial location.
	#Returns a 1D array for plotting and performs unit conversion.

		#If no mesh sizes supplied, collect sizes for current global folder.
		if Rmesh == 'NaN' or Isym == 'NaN':
			Rmesh,Isym = R_mesh[l],ISYMlist[l]
		#endif

		#Obtain start location for requested data and perform SI conversion.
		ZStart = int(Rmesh*Profile)
		ZEnd = int(Rmesh*(Profile+1))

		#Plot lines for each variable at each requested slice, ignoring ones that fail.
		#If mesh is symmetric, copy the data over and make full plot.
		if int(Isym) == 1:
			Zend = len(Data[process])-ZStart
			Zstart = len(Data[process])-ZEnd
			RProfile = Data[process][Zstart:Zend][::-1]
			#If variable is radially symmetric, add negative to symmetry
			if is_radial_variable(variable) == True:
				for m in range(0,len(RProfile)):
					RProfile.append(-Data[process][Zstart:Zend][m])
				#endfor
			#If variable is axially symmetric, add positive.
			elif is_radial_variable(variable) == False:
				for m in range(0,len(RProfile)):
					RProfile.append(Data[process][Zstart:Zend][m])
				#endfor
			#endif
			RProfile = RProfile[::-1]	#Reverse index, negative first then positive values.

		#If the data isn't symmetric, just plot as is.
		elif int(Isym) == 0:
			Zend = len(Data[process])-ZStart
			Zstart = len(Data[process])-ZEnd
			RProfile = Data[process][Zstart:Zend]
		#endif

		#Convert units if required.								!!! RM SJD, CONVERSION PERFORMED AT READ-IN
	#	RProfile = VariableUnitConversion(RProfile,variable, Units, AtomicSpecies)

		return(RProfile)
	#enddef


	#=========================#
	#=========================#


	def ExtractAxialProfile(Data,process,variable,lineout,Rmesh=0,Zmesh=0,Isym=0):
	#Obtains an axial 1D profile at a requested radial location.
	#Returns a 1D array for plotting and performs unit conversion.

		#If no mesh sizes supplied, collect sizes for current global folder.
		if Rmesh == 0 or Zmesh == 0 or Isym == 0:
			Rmesh,Zmesh,ISym = R_mesh[l],Z_mesh[l],ISYMlist[l]
		#endif

		#Pull out Z-data point from each radial line of data and list them.
		Zlineout = list()
		for i in range(0,int(Zmesh)):
			datapoint = Rmesh*i + lineout
			try:
				Zlineout.append(Data[int(process)][int(datapoint)])
			except:
				break
			#endtry
		#endfor

		#Convert units if required.								!!! RM SJD, CONVERSION PERFORMED AT READ-IN
	#	Zlineout = VariableUnitConversion(Zlineout,variable, Units, AtomicSpecies)

		return(Zlineout)
	#enddef


	#=========================#
	#=========================#


	def WaveformLoc(location,origintype):
	#Convert cell location for use with WaveformExtractor function.
	#Returns mesh location based on input string or [R,Z] list.
	#Rcell,Zcell = Waveformloc(electrodeloc,'Phase')

		#if data is from movie1, Zaxis is back to front.
		if len(location) == 2 and origintype == 'Phase':
			RlineoutLoc = (Z_mesh[l]-location[1])-1
			ZlineoutLoc = location[0]
		#if data is from TECPLOT2D, Zaxis is correct orientation.
		elif len(location) == 2 and origintype == '2D':
			RlineoutLoc = location[1]
			ZlineoutLoc = location[0]
		else:
			RlineoutLoc = Z_mesh[l]/2
			ZlineoutLoc = 0
		#endif

		return(RlineoutLoc,ZlineoutLoc)
	#enddef


	#=========================#
	#=========================#


	def WaveformExtractor(PhaseData,PPOT,waveformlocation=electrodeloc):
	#Takes phasedata for current folder and PPOT process number.
	#Returns two arrays: VoltageWaveform at electrode location and average waveform.
	#VoltageWaveform,WaveformBias = WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[0])

		#Create required lists and extract electrode location.
		VoltageWaveform,WaveformBias = list(),list()
		RLoc = WaveformLoc(waveformlocation,'Phase')[0]
		ZLoc = WaveformLoc(waveformlocation,'Phase')[1]

		#Create voltage waveform for requested integer number of phasecycles
		for i in range(0,int(phasecycles*len(PhaseData))):
			Index = i % len(PhaseData)		#Use modulo index for additional phase cycles
			VoltageWaveform.append(ExtractAxialProfile(PhaseData[Index],PPOT,'PPOT',ZLoc)[RLoc])
		#endfor

		#Calculate time averaged waveform bias, i.e. waveform symmetry.
		for m in range(0,len(VoltageWaveform)):
			WaveformBias.append(sum(VoltageWaveform)/len(VoltageWaveform))
		#endfor

		#Calculate maximum positive and negative waveform amplitudes and compute average Vpp
		PositiveAmp,NegativeAmp = max(VoltageWaveform),min(VoltageWaveform)
		PeakToPeakVoltage = abs(PositiveAmp)+abs(NegativeAmp)

		return(VoltageWaveform,WaveformBias,[PositiveAmp,NegativeAmp,PeakToPeakVoltage])
	#enddef


	#=========================#
	#=========================#


	def TrendAtGivenLocation(probeloc,process,variable):
	#Trend analysis for a given point on a 2D Image.
	#Takes global 'probeloc' for safety, process and variable.
	#Returns two arrays: One is the X-axis to plot against
	#Second is the value of variable at location for all simulations.

		#Refresh lists that change per image.
		R,Z = probeloc[0],probeloc[1]
		Trend = list()
		Xaxis = list()

		#For all simulation folders.
		for l in range(0,numfolders):

			#Extract image with given process and variable name.
			Image = ImageExtractor2D(Data[l][process],variable,R_mesh[l],Z_mesh[l])

			#Update X-axis with folder information.
			Xaxis.append(folder_name_trimmer(Dirlist[l]))

			#TREND ANALYSIS - Value at Given Location Comparison.
			try: Trend.append( Image[Z][R] )
			except: Trend.append(float('NaN'))

			#Display Min/Max value trends to terminal if requested.
			if print_generaltrends == True:
				Location = '('+str(round(R*dr[l],1))+'cm,'+str(round(Z*dz[l],1))+'cm)'
				print(folder_name_trimmer(Dirlist[l]))
				print( str(variable)+' @ '+Location+':', round(Trend[-1], 5))
			if print_generaltrends == True and l == numfolders-1:
				print( '')
			#endif
		#endfor

		#Normalise to maximum value in each profile if required.
		if image_normalise == True:
			Trend,Min,Max = Normalise(Image)
		#endif

		return(Xaxis,Trend)
	#enddef


	#=========================#
	#=========================#


	def MinMaxTrends(lineout,Orientation,process):
	#General trend plotting function for use with multiple folders.
	#Takes a lineout location and orientation string as input.
	#And returns the maximum and minimum values For all folders to be analysed.
	#With an associated X-axis composed of the trimmed folder names.

		#Refresh lists that change per profile.
		MaxValueTrend, MinValueTrend = list(), list()
		Xaxis = list()
		p = process

		#For each folder in the directory.
		for l in range(0,numfolders):

			#Create and correct VariableIndices for each folder as required.
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])
			VariableIndices,VariableStrings = variable_interpolator(VariableIndices, VariableStrings, MinSharedVariables, Globalnumvars)

			#Update X-axis with folder information.
			Xaxis.append(folder_name_trimmer(Dirlist[l]))

			#Obtain radial and axial profiles for further processing.
			if Orientation == 'Radial':
				try: Profile = ExtractRadialProfile(Data[l],VariableIndices[p],VariableStrings[p],lineout,R_mesh[l],ISYMlist[l])
				except: Profile = float('NaN')
				#endtry
			elif Orientation == 'Axial':
				try: Profile = ExtractAxialProfile(Data[l],VariableIndices[p],VariableStrings[p],lineout,R_mesh[l],Z_mesh[l],ISYMlist[l])
				except: Profile = float('NaN')
				#endtry
			#endif

			#TREND ANALYSIS - Maximum/Minimum Value Comparison.
			try: MaxValueTrend.append(max(Profile))
			except: MaxValueTrend.append(float('NaN'))
			try: MinValueTrend.append(min(Profile))
			except: MinValueTrend.append(float('NaN'))
			#endtry

			#Display Min/Max value trends to terminal if requested.
			if print_generaltrends == True:
				VariableName = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[p]
				print(folder_name_trimmer(Dirlist[l]))
				print( VariableName+' '+Orientation+'Maximum: ', round(max(MaxValueTrend), 5))
				print( VariableName+' '+Orientation+'Minimum: ', round(min(MinValueTrend), 5))
			if print_generaltrends == True and l == numfolders-1:
				print( '')
			#endif
		#endfor

		#Normalise to maximum value in each profile if required.
		if image_normalise == True:
			MaxValueTrend,MaxMin,MaxMax = Normalise(MaxValueTrend)
			MinValueTrend,MinMin,MinMax = Normalise(MinValueTrend)
		#endif

		return(Xaxis,MaxValueTrend,MinValueTrend)
	#enddef


	#=========================#
	#=========================#

	#TREND ANALYSIS - Speed of Sound
	def CalcSoundSpeed(NeutralDensity,Pressure,Dimension='2D'):
	#Calculates local sound speed via Newton-Laplace equation
	#Takes appropriate neutral density [m-3] and pressure [Torr] (0D,1D or 2D)
	#Returns same dimensionality array of sound speeds in m/s
	#SoundSpeed = CalcSoundSpeed(ArgonDensity,Pressure,Dimension='2D')

		#Initiate required lists and set atomic values
		SoundSpeedArray = list()
		AdiabaticIndex = 5.0/3.0		#Assumes diatomic species (Hydrogen,Helium,Argon, etc...)
		AtomicMass = 39.948*1.66E-27	#[Kg]		#NB HARDCODED FOR ARGON

		#For 0D values:
		if Dimension == '0D':
			ElasticityModulus = AdiabaticIndex*Pressure*133.33					#[Pa] = [kg m-1 s-2]
			MassDensity = NeutralDensity*AtomicMass								#[kg m-3]

			#Calculate local sound speed via Newton-Laplace equation
			try: SoundSpeedArray = np.sqrt( ElasticityModulus/MassDensity )		#[m/s]
			except: SoundSpeedArray = np.nan
		#endif

		#For 1D arrays:
		if Dimension == '1D':
			for i in range(0,len(NeutralDensity)):
				ElasticityModulus = AdiabaticIndex*Pressure[i]*133.33			#[Pa] = [kg m-1 s-2]
				MassDensity = NeutralDensity[i]*AtomicMass						#[kg m-3]

				#Calculate local sound speed via Newton-Laplace equation
				try: SoundSpeed = np.sqrt( ElasticityModulus/MassDensity )		#[m/s]
				except: SoundSpeed = np.nan
				SoundSpeedArray.append( SoundSpeed )
			#endfor
		#endif

		#For 2D arrays:
		if Dimension == '2D':
			for i in range(0,len(NeutralDensity)):
				SoundSpeedArray.append(list())
				for j in range(0,len(NeutralDensity[i])):
					ElasticityModulus = AdiabaticIndex*Pressure[i][j]*133.33	#[Pa]
					MassDensity = NeutralDensity[i][j]*AtomicMass				#[kg m-3]

					#Calculate local sound speed via Newton-Laplace equation
					try: SoundSpeed = np.sqrt( ElasticityModulus/MassDensity )	#[m/s]
					except: SoundSpeed = np.nan
					SoundSpeedArray[i].append( SoundSpeed )
				#endfor
			#endfor
		#endif

		return(SoundSpeedArray)
	#enddef


	#=========================#
	#=========================#

	#TREND ANALYSIS - DCbias
	def DCbiasMagnitude(PPOTlineout):
	#Takes a PPOT profile and calcuates DCbias via difference in voltage drop.
	#Can identify DC-bias for parallel plate discharges and dielectric discharges.

		#Identify if radial or axial.
		if len(PPOTlineout) == Z_mesh[l]:
			electrodelocation = WaveformLoc(electrodeloc,'2D')[1]
		elif len(PPOTlineout) in [R_mesh[l],R_mesh[l]*2]:
			electrodelocation = WaveformLoc(electrodeloc,'2D')[0]
		#endif

		#Identify Min/Max Potential magnitudes and location of max potential.
		MinPPOT,MaxPPOT = min(PPOTlineout),max(PPOTlineout)
		MaxIndex = np.argmax(PPOTlineout)

		#Split PPOT profile into each sheath, pre and post max potential
		PreIndex = PPOTlineout[:MaxIndex]
		PostIndex = PPOTlineout[MaxIndex:]


		##=========================================##

		#Metals have flat PPOT profiles, dielectric/plasma have gradients.
		MetalIndices = list([0])
		DielectricIndices = list()
		for i in range(0,len(PPOTlineout)-1):
			if PPOTlineout[i] == PPOTlineout[i+1]:
				MetalIndices.append(i)
			elif PPOTlineout[i] != PPOTlineout[i+1]:
				DielectricIndices.append(i)
			#endif
		#endfor
		MetalIndices.append(len(PPOTlineout)-1)

		#Grounded metal will have a PPOT of zero -- ##INCORRECT IF DC-BIAS == int(0.0)##
		GMetalIndices = list()
		for i in range(0,len(MetalIndices)):
			if PPOTlineout[MetalIndices[i]] == 0:
				GMetalIndices.append(MetalIndices[i])
			#endif
		#endfor

		#Any metal that is not grounded will be powered -- ##BAD ASSUMPTION FOR ALL MESHES##
		PMetalIndices = list()
		for i in range(0,len(MetalIndices)):
			if MetalIndices[i] not in GMetalIndices:
				PMetalIndices.append(MetalIndices[i])
			#endif
		#endfor

		##=========================================##


		#Identify voltage drop through each sheath from max potential.
		try: PreIndexVoltageDrop = MaxPPOT - min(PreIndex)
		except: PreIndexVoltageDrop = MaxPPOT
		#endtry
		try: PostIndexVoltageDrop = MaxPPOT - min(PostIndex)
		except: PostIndexVoltageDrop = MaxPPOT
		#endtry

		#Minimum voltage is not one of the electrodes - "Dielectric Discharge"
		if min(PPOTlineout) not in [ PPOTlineout[0],PPOTlineout[-1] ]:
			try: DCbias = MinPPOT
			except: DCbias = MaxPPOT
			#endtry

		#Minimum voltage is one of the electrodes - "Parallel Plate Discharge"
		else:
			try: DCbias = PPOTlineout[PMetalIndices[0]]
			except: DCbias = PreIndexVoltageDrop - PostIndexVoltageDrop
		#endif

		if IDEBUG == True:
			X1 = range(0,len(PreIndex))
			X2 = range(len(PreIndex),len(PPOTlineout))

			plt.plot(X1,PreIndex, lw=2)
			plt.plot(X2,PostIndex, lw=2)
			plt.plot(np.argmax(PPOTlineout),max(PPOTlineout), 'go',  ms=12)
			for i in range(0,len(GMetalIndices)):
				plt.plot(GMetalIndices[i],PPOTlineout[GMetalIndices[i]], 'ko',  ms=12)
			#endfor
			for i in range(0,len(PMetalIndices)):
				plt.plot(PMetalIndices[i],PPOTlineout[PMetalIndices[i]], 'ro',  ms=12)
			#endfor

			plt.xlabel('Cell Number')
			plt.ylabel('Voltage [V]')
			plt.legend(['PreBulk','PostBulk','Bulk'])
			plt.title('DCBIAS_DEBUG'+str(l+electrodelocation))
			plt.savefig(DirTrends+'DCBIAS_DEBUG'+str(l+electrodelocation)+'.png')
			clearfigures(fig)
		#endif

		return DCbias
	#enddef


	#=========================#
	#=========================#


	#BRINKMANN SHEATH WIDTH CALCULATOR
	def CalcSheathExtent(folderidx,Orientation='Radial',Phase='NaN',Ne=list(),Ni=list()):
	#Calculates Brinkmann sheath width assuming Child-Langmuir conditions.
	#Calculation Methods: 'AbsDensity', 'IntDensity'
	#Takes current folderidx, current axis, movie1 Phase and sheath calc method.
	#Returns array of sheath distances from symmetry boundary (or from origin)
	#Sx,SxAxis = CalcSheathExtent(folderidx=l,Phase=moviephaselist[k])

		#Initiate required data storage lists
		Eproc,IONproc = list(),list()
		NPos,NNeg = list(),list()
		SxAxis,Sx = list(),list()

		#Import global sheath calculation method, ratio threshold, and charged species names
		Ratio_Threshold = Sheath_IonRatio_Threshold
		SheathMethod=GlobSheathMethod
		if len(SheathIonSpecies) == 0:
			global PosSpecies
			global NegSpecies
		#Force a single sheath species - Legacy Code or for testing purposes
		elif len(SheathIonSpecies) > 0:
			PosSpecies = SheathIonSpecies
			NegSpecies = ['E']
		#endif

		MinimumSpeciesExist = False
		#Check relevant datafile to confirm that positive and negative species exist
		AllSpecies = PosSpecies + NegSpecies
		if all(Species in Header_TEC2D[folderidx] for Species in AllSpecies):
			MinimumSpeciesExist = True
		#endif

		#Return NaN sheath array if no appropriate species are detected
		if MinimumSpeciesExist == False:
			SxAxis = GenerateAxis(Orientation,Isym=ISYMlist[l])
			Sx = np.empty(len(SxAxis))
			[np.nan for x in Sx]

			return(Sx,SxAxis)
		#endif

		#Return NaN sheath array if diagnostic is not requested
		if image_plotsheath not in ['Radial','Axial']:
			SxAxis = GenerateAxis(Orientation,Isym=ISYMlist[l])
			Sx = np.empty(len(SxAxis))
			[np.nan for x in Sx]

			return(Sx,SxAxis)
		#endif

		#Identify charged species and alter names to suit TECPLOT2D nomenclature
		for i in range(0,len(PosSpecies)): PosSpecies[i] = PosSpecies[i] = PosSpecies[i].replace('^','+')
		for i in range(0,len(NegSpecies)): NegSpecies[i] = NegSpecies[i] = NegSpecies[i].replace('^','-')
		if 'E' in NegSpecies: NegSpecies.remove('E')			#Might Cause An Issue With Global!!!

		#=====#=====#

		#=======#	#=======#	#=======#
		#=======#	#=======#	#=======#

		#ISSUE WITH THIS REGARDING THE DATA READ-IN
		#PREVIOUS VERSION SENDS Ne and Ni EXPLICITLY INTO FUNCTION,
		#IDEALLY WOULD HAVE ALL OF THIS RUN ONCE, EXTRACTING THE DATA FROM RAW_PHASEDATA
		#PASSING THE REQUIRED Ne and Neff INTO THE 2ND PART OF THE FUNCTION.
		#IF Ne, Neff don't exist: Run phase-extraction-function
		#Else: run old code below, using global values

		#Obtain current folder ion and electron densities if not already supplied.
		#Default to 2D data format.
	#	if Phase == 'NaN' and len(Ne) == 0:
			#Obtain electron density and extract 2D image for further processing.
	#		Eproc = EnumerateVariables(['E'],Header_TEC2D[folder])[0][0]
	#		Ne = ImageExtractor2D( Data[folder][Eproc] )

			#Obtain all positive and negative ion densities and extract 2D images for further processing
	#		PosSpeciesproc = EnumerateVariables(PosSpecies,Header_TEC2D[folder])[0]
	#		for i in range(0,len(PosSpeciesproc)):
	#			NPos.append( ImageExtractor2D(Data[folder][PosSpeciesproc[i]]) )
			#endfor
	#		NegSpeciesproc = EnumerateVariables(NegSpecies,Header_TEC2D[folder])[0]
	#		for i in range(0,len(NegSpeciesproc)):
	#			NNeg.append( ImageExtractor2D(Data[folder][NegSpeciesproc[i]]) )
			#endfor

		#If phase is supplied, use phase data format.  (Proc=Proc-2 to skip R,Z data in phase data)
	#	elif Phase != 'NaN' and len(Ne) == 0:
	#		Eproc = EnumerateVariables(['E'],Header_movie1[folder])[0][0]
	#		Ne = ImageExtractor2D( PhaseMovieData[folder][Phase][Eproc-2] )

			#Obtain all positive and negative ion densities and extract 2D images for further processing
	#		PosSpeciesproc = EnumerateVariables(PosSpecies,Header_movie1[folder])[0]
	#		for i in range(0,len(PosSpeciesproc)):
	#			NPos.append( ImageExtractor2D(PhaseMovieData[folder][Phase][PosSpeciesproc[i]-2]) )
			#endfor
	#		NegSpeciesproc = EnumerateVariables(NegSpecies,Header_movie1[folder])[0]
	#		for i in range(0,len(NegSpeciesproc)):
	#			NNeg.append( ImageExtractor2D(PhaseMovieData[folder][Phase][NegSpeciesproc[i]-2]) )
			#endfor

		#If specific electron and ion species densities are supplied, use those
	#	elif len(Ne) > 0 or len(Ni) > 0:
	#		Ne = ImageExtractor2D( Ne ) 		#Ne[i][j]
	#		NPos = [ ImageExtractor2D( Ni ) ]	#Put in array []    (NPos[k][i][j])
	#		PosSpeciesproc = ['Ion+']			#Set length to 1
	#		NegSpeciesproc = []					#Set length to 0
		#endif

		#Combine 2D images of all positive ion species densities and all negative ion species densitiies
	#	NPos = [[sum(x) for x in zip(NPos[0][i],NPos[1][i])] for i in range(len(NPos[0]))]
	#	NNeg = [[sum(x) for x in zip(NNeg[0][i],NNeg[1][i])] for i in range(len(NNeg[0]))]
	#							HOW TO ZIP ARBITARY NUMBER OF ARRAYS?
	#	TotNPos = np.zeros( (len(Ne),len(Ne[0])) ).tolist()
	#	for i in range(0,len(TotNPos)):
	#		for j in range(0,len(TotNPos[0])):
	#			for k in range(0,len(PosSpeciesproc)): TotNPos[i][j] += NPos[k][i][j]
				#endfor
			#endfor
		#endfor
	#	TotNNeg = np.zeros( (len(Ne),len(Ne[0])) ).tolist()
	#	for i in range(0,len(TotNNeg)):
	#		for j in range(0,len(TotNNeg[0])):
	#			for k in range(0,len(NegSpeciesproc)): TotNNeg[i][j] += NNeg[k][i][j]
				#endfor
			#endfor
		#endfor

		#Determine effective positive ion density as: Neff = sum(Total NPos)-sum(Total NNeg)
	#	Neff = np.zeros( (len(Ne),len(Ne[0])) ).tolist()
	#	for i in range(0,len(Neff)):
	#		for j in range(0,len(Neff[0])):
	#			Neff[i][j] = TotNPos[i][j] - TotNNeg[i][j]
			#endfor
		#endfor

		#=======#	#=======#	#=======#
		#=======#	#=======#	#=======#

		#!!! OLD METHOD !!!
		#Obtain current folder ion and electron densities if not already supplied.
		#Default to 2D data format.
		if Phase == 'NaN' and len(Ne) == 0:
			IONproc = enumerate_variables(PosSpecies, Header_TEC2D[folderidx])[0][0]
			Eproc = enumerate_variables(['E'], Header_TEC2D[folderidx])[0][0]
			Ne = Data[folderidx][Eproc]
			Ni = Data[folderidx][IONproc]
		#If phase is supplied, use phase data format.				<<< TO BE REMOVED
	#	elif Phase != 'NaN' and len(Ne) == 0:
	#		IONproc = EnumerateVariables(PosSpecies,Header_movie1[folderidx])[0][0]
	#		Eproc = EnumerateVariables(['E'],Header_movie1[folderidx])[0][0]
	#		Ne = PhaseMovieData[folderidx][Phase][Eproc]
	#		Ni = PhaseMovieData[folderidx][Phase][IONproc]
		#endif
		#Extract 2D image for further processing.
		Ne,Neff = ImageExtractor2D(Ne),ImageExtractor2D(Ni)
		#!!! OLD METHOD !!!

		#=======#	#=======#	#=======#
		#=======#	#=======#	#=======#

		#Radial sheath array (Sx) is calculated employing axial integrations for all radial locations
		if Orientation == 'Radial':

			#Determine sheath edge through integration of charge density:
			#Sheath extension: integral_(R0->Rwall) ne dR == integral_(Rwall->R0) ni dR
			if SheathMethod == 'IntDensity':

				#Find radial index where plasma zone ends
				for i in range(0,len(Neff)):
					#Define wall radius to integrate ions into bulk from.
					for j in range(0,len(Neff[i])):

						#if ion density drops to zero, we've hit a material surface.
						if Neff[i][j] == 0.0 and j == 0:
							RadialPlasmaExtent = 0
							break
						elif Neff[i][j] == 0.0 and j > 0:
							RadialPlasmaExtent = j-1
							break
						#endif
					#endfor
		#			RadialPlasmaExtent = len(Neff[i])	#DEBUG OPTION: Sets RadialPlasmaExtent to max for all Z

					#No plasma, all radii are solids, append 'nan' to avoid plotting.
					if RadialPlasmaExtent == 0:
						Sx.append(np.nan)								#[cm]
					#If non-zero plasma extent, determine radial cell satisfying Brinkmann Criterion
					elif RadialPlasmaExtent > 0:
						#Refresh sums after every radial profile.
						Neff_sum,Ne_sum = 0.0,0.0
						for j in range(0,RadialPlasmaExtent):
							#Sum density radially for ions and electrons.
							reversed_j = RadialPlasmaExtent-j-1
							Neff_sum += Neff[i][reversed_j]		#Sum from R=wall to R=0	[reversed_j]
							Ne_sum += Ne[i][j]					#Sum from R=0 to R=wall [j]

							#If ion sum is greater than electron, sheath has begun.
							if Neff_sum/Ne_sum >= Ratio_Threshold:		#~1.00 - 1.03 is good
								Sx.append(j*dr[l])						#[cm]
								break
							#If no sheath found within plasma region, append wall location (i.e. Sx=Rwall)
							if j == (RadialPlasmaExtent-1):
								Sx.append((RadialPlasmaExtent+1)*dr[l])	#[cm]
	#							Sx.append(np.nan)						#[cm] Nice Plots/Breaks Statistics
							#endif
						#endfor
					#endif
				#endfor

			#==========#

			#Determine sheath edge by 'instantaneous' charge density:
			elif SheathMethod == 'AbsDensity':
				#Sheath extension: ni @R >= ne @R, simplified model.
				for Z in range(0,len(Neff)):
					for R in range(0,len(Neff[Z])):

						#Sheath starts when ion density exceeds electron density.
						if Neff[Z][R]/Ne[Z][R] >= Ratio_Threshold:		#~1.00 - 1.03 is good
							Sx.append(R*dr[l])
							break
						#If no sheath found, append 'NaN' to avoid plotting.
						if R == (len(Neff[Z])-1):
	#						Sx.append(0.0)
							Sx.append(np.nan)
						#endif
					#endfor
				#endfor
			#endif

			#Extract axis to plot sheath against
			#RADIAL sheath extent is plotted against AXIAL axis
			SxAxis = GenerateAxis('Axial',Isym=ISYMlist[l])
		#endif

		#=======#	#=======#	#=======#
		#=======#	#=======#	#=======#

		#Axial sheath array (Sx) is calculated exmploying radial integrations for all axial locations
		if Orientation == 'Axial':

			#Determine sheath edge through integration of charge density:
			if SheathMethod == 'IntDensity':

				#Find axial index where plasma zone ends
				for i in range(0,len(Neff)):
					#Define wall height to integrate ions into bulk from.
					for j in range(0,len(Neff[i])):

						#if ion density drops to zero, we've hit a material surface.
						if Neff[i][j] == 0.0 and j == 0:						#### SJD TO BE COMPLETED ####
							AxialPlasmaExtent = 0
							break
						elif Neff[i][j] == 0.0 and j > 0:						#### SJD TO BE COMPLETED ####
							AxialPlasmaExtent = j-1
							break
						#endif
					#endfor
		#			AxialPlasmaExtent = len(Neff)	#DEBUG OPTION: Sets AxialPlasmaExtent to max for all R
				#endfor

				#Sheath extension: integral_(R0->Rwall) ne dR == integral_(Rwall->R0) ni dR
				for i in range(0,len(Neff[0])):									#### SJD TO BE COMPLETED ####
					Sx.append(np.nan)
				#endfor

			#Determine sheath edge by 'instantaneous' charge density:
			elif SheathMethod == 'AbsDensity':
				#Sheath extension: ni @Z >= ne @Z, simplified model.
				for R in range(0,len(Neff[0])):
					for Z in range(0,len(Neff)):

						#Sheath starts at axial index where ion density exceeds electron density.
						# NOTE: Inequality is reversed as axial method calculates from the sheath edge
						#		to the bulk, and therefore Neff/Ne REDUCES from a maximum value.
						if Neff[Z][R]/Ne[Z][R] <= Ratio_Threshold:		# RM NEED TO CHECK IF 1.03 IS GOOD...
							Sx.append(Z*dz[l])							# ... OR IF IT NEEDS TO BE BELOW 1
							break
						#If no sheath found over all axial profiles, append 'NaN' to avoid plotting.
						elif Z == (len(Neff)-1):
	#						Sx.append(0.0)
							Sx.append(np.nan)
						#endif
					#endfor
				#endfor
			#endif

			#Extract axis to plot sheath against
			#AXIAL sheath extent is plotted against RADIAL axis
			SxAxis = GenerateAxis('Radial',0)							#### Isym=ISYMlist[l])
		#endif

		#=======#	#=======#	#=======#
		#=======#	#=======#	#=======#

		#Print sheath characteristics if requested.
		if print_sheath == True:
			#Determine location of interest (loc) for diagnostic
			if Orientation == 'Axial': loc = electrodeloc[0]		#Radial point of interest
			if Orientation == 'Radial': loc = electrodeloc[1]		#Axial point of interest

			#Obtain SheathWidth at electrodeloc
			try: SheathWidth = round(Sx[loc],3)
			except: SheathWidth = 0.0
			print( 'Simulation:', Dirlist[folder])
			print( 'Sheath Location:',SheathWidth*10, 'mm')
			print( 'Sheath Extent:',((sourcewidth[0]*dr[l])-SheathWidth)*10, 'mm')
		#endif

		#Delete any intermediate arrays
		del(Eproc,IONproc,Ne,Ni,Neff)

		#Return sheath axis and sheath extent
		return(Sx,SxAxis)
	#enddef


	#=========================#
	#=========================#


	def PlotSheathExtent(SxAxis,Sx,ax=plt.gca(),ISymmetry=0,Orientation='Radial'):
	#PlotSheathExtent(SxAxis,Sx,ax[0],ISYMlist[l],Orientation=image_plotsheath)

		#Create any require symmetry arrays
		SymSx = list()
		SymSxAxis = list()

		#=====#=====#	#=====#=====#

		if Orientation == 'Axial' and image_rotate == True:

			#Determine mesh symmetry and create symmetric mesh if required
			if int(ISymmetry) == 1:
				#Symmetric sheath creation
				for i in range(0,len(Sx)):
					try: SymSxAxis.append(-SxAxis[i])
					except: SymSxAxis.append(np.nan)
				#endfor
			#endif

			#Plot AXIAL SHEATH EXTENT against the X-AXIS
			X1,Y1 = Sx,SxAxis
			X2,Y2 = Sx,SymSxAxis

		#=====#=====#

		elif Orientation == 'Axial' and image_rotate == False:

			#Determine mesh symmetry and create symmetric mesh if required
			if int(ISymmetry) == 1:
				#Symmetric sheath creation
				for i in range(0,len(Sx)):
					try: SymSxAxis.append(-SxAxis[i])
					except: SymSxAxis.append(np.nan)
				#endfor
			#endif

			#Plot AXIAL SHEATH EXTENT against Y-AXIS
			X1,Y1 = SxAxis,Sx
			X2,Y2 = SymSxAxis,Sx

		#=====#=====#

		elif Orientation == 'Radial' and image_rotate == True:

			#Determine mesh symmetry and create symmetric mesh if required
			if int(ISymmetry) == 1:
				#Symmetric sheath creation
				for i in range(0,len(Sx)):
					try: SymSx.append(-Sx[i])
					except: SymSx.append(np.nan)
				#endfor
			#endif

			#plot RADIAL SHEATH EXTENT against X-AXIS
			X1,Y1 = SxAxis,Sx
			X2,Y2 = SxAxis,SymSx

		#=====#=====#

		elif Orientation == 'Radial' and image_rotate == False:

			#Determine mesh symmetry and create symmetric mesh if required
			if int(ISymmetry) == 1:
				#Symmetric sheath creation
				for i in range(0,len(Sx)):
					try: SymSx.append(-Sx[i])
					except: SymSx.append(np.nan)
				#endfor
			#endif

			#Plot RADIAL SHEATH EXTENT against Y-AXIS
			X1,Y1 = Sx,SxAxis
			X2,Y2 = SymSx,SxAxis
		#endif

		#=====#=====# 	#=====#=====#

		#Plot sheath extent from origin with or without mesh symmetry
		if image_plotsheath in ['Radial','Axial'] and int(ISymmetry) == 0:
			ax.plot(X1,Y1, 'w--', lw=2)
		elif image_plotsheath in ['Radial','Axial'] and int(ISymmetry) == 1:
			ax.plot(X1,Y1, 'w--', lw=2)
			ax.plot(X2,Y2, 'w--', lw=2)
		#endif

		return()
	#enddef


	#=========================#
	#=========================#

	#====================================================================#
	#====================================================================#





































	#====================================================================#
					  #IMAGE PLOTTING AND DATA ANALYSIS#
	#====================================================================#

	#====================================================================#
					   #STEADY STATE TECPLOT2D FIGURES#
	#====================================================================#

	image_rotate = False

	# Generate and save image of required variable for given mesh size.
	if savefig_tecplot2D == True:

		for l in range(0,numfolders):
			# Create new folder to keep output plots.
			Dir2Dplots = create_new_folder(Dirlist[l], 'TECPlot2D')

			# Create VariableIndices for each folder as required.
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])

			# Define image extent for directory 'l'
			extent, aspectratio = DataExtent(l)

			# Reshape specific part of 1D Data array into 2D image for plotting.
			for i in tqdm(range(0,len(VariableIndices))):

				# Extract 2D image from TECPLOT2D "Data" for folder "l"
				Image = ImageExtractor2D(Data[l][VariableIndices[i]],VariableStrings[i])

				# Create figure and plot TECPLOT2D image
				fig,ax = figure(aspectratio)
				fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,VariableStrings[i],fig=fig,ax=ax)

				#=====##=====# IMAGE OVERLAYS #=====##=====#

				# Overlay vector streamplot if exists
				Radial,Axial = enumerate_vectors(VariableStrings[i], Header_TEC2D[l])

				# Confirm Variable[k] has both vector counterparts (e.g. FR-AR3S, FZ-AR3S)
				VectorVariablesExist = True
				if None in Radial or None in Axial: VectorVariablesExist = False

				# Overlay variable[k]'s vector components onto 2D flood plot
				if image_plotvector == True and VectorVariablesExist == True:

					# Create meshgrid
					R_Space = np.linspace(0,Radius[l],R_mesh[l])			# Radial mesh locations [mm]
					Z_Space = np.linspace(0,Height[l],Z_mesh[l])			# Axial mesh locations [mm]
					if image_rotate == True:
						R_Space = np.linspace(0,Height[l],Z_mesh[l])		# Radial location is now axial location
						Z_Space = np.linspace(0,Radius[l],R_mesh[l])		# Axial location is now radial location
					#endif
					R, Z = np.meshgrid(R_Space, Z_Space)
						
					# See definition of "Radial" and "Axial" above
					UR = ImageExtractor2D(Data[l][Radial[1]],Radial[0])  		# Radial vector magnitude
					UZ = ImageExtractor2D(Data[l][Axial[1]] ,Axial[0] )  		# Axial vector magnitude
					if image_rotate == True:
						UR = ImageExtractor2D(Data[l][Axial[1]] ,Axial[0])  	# Radial magnitude DATA is now Axial
						UZ = ImageExtractor2D(Data[l][Radial[1]],Radial[0] )  	# Axial magnitude DATA is now Radial
						if image_plotsymmetry == False:
							UR = UR.transpose()[::-1]							# Without symmetry, Axis Zero is on LHS
							UZ = UZ.transpose()[::-1]
						elif image_plotsymmetry == True:
							UR = UR.transpose()									# With symmetry, Axis Zero is on RHS
							UZ = UZ.transpose()
						#endif
					#endif
					VectorLength = np.sqrt(UR**2 + UZ**2)

					# Streamplot fails if provided "zeros", inform user and skip to next variable
					try:
						Density = image_vectordensity
						Linewidth = image_vectorlw
						ax.streamplot(R, Z, UR, UZ, color=VectorLength, cmap='viridis', density=Density, linewidth=Linewidth)
					except:
						print('Warning: Invalid streamplot for variable: '+VariableStrings[i])
					#endtry
				#endif

				#=====#

				# Compute sheath extent and overlay onto figure if requested
				Sx,SxAxis = CalcSheathExtent(folderidx=l,Orientation=image_plotsheath)
				if image_plotsheath in ['Radial','Axial']:
					PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l],Orientation=image_plotsheath)
				#endif

				#=====#

				# Overlay location of 1D profiles if requested, adjusting for image rotation.
				if image_plotoverlay == True:
					for j in range(0,len(radialprofiles)):
						X1,X2 = extent[0],extent[1]
						Y1,Y2 = radialprofiles[j]*dz[l],radialprofiles[j]*dz[l]
						if image_rotate == True: X1,X2,Y1,Y2 = Y1,Y2,X1,X2
						ax.plot((X1,X2),(Y1,Y2),'k--',lw=2)
					#endfor
					for j in range(0,len(axialprofiles)):
						X1,X2 = axialprofiles[j]*dr[l],axialprofiles[j]*dr[l]
						Y1,Y2 = extent[2],extent[3]
						if image_rotate == True: X1,X2,Y1,Y2 = Y1,Y2,X1,X2
						ax.plot((X1,X2),(Y1,Y2),'k--',lw=2)
					#endfor
				#endif

				#=====##=====# Image Beautification #=====##=====#

				# Define title, labels, etc...
				Title = '2D Steady State Plot of '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
				if image_rotate == True:	Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
				elif image_rotate == False:	Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				#endif

				# Add Colourbar (must happen before image cropping!)
				cbarlabel = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
				cax = Colourbar(ax,cbarlabel[i],image_cbarbins,Lim=CbarMinMax(ax,Image))

				# Crop image dimensions to [image_radialcrop,image_axialcrop]
				# Also resets cbar min/max to cropped region min/max
				if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
					CropImage(ax,Rotate=image_rotate)
				#endif

				# Apply mesh geometry to image
				ImageGeometry(fig,ax,image_plotsymmetry)

				# Apply image options, and enforce overides if requested
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)


				#=====##=====# Image I/O #=====##=====#

				# Save Figure
				plt.savefig(Dir2Dplots+'2DPlot_'+VariableStrings[i]+ext)
				clearfigures(fig)

				#=====#

				# Write data underpinning current figure in .csv format
				if Write_CSV == True:
					CSVDir = create_new_folder(Dir2Dplots, '2DPlots_Data')
					CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
					CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
					CSVFilename = VariableStrings[i]+'.csv'
					CSVTitle = str(Dirlist[l])
					CSVLabel = str(cbarlabel[i])
					CSVISYM = 'ISYM='+str(ISYMlist[l])
					CSVRotate = 'Rotate='+str(image_rotate)
					CSVHeader = [CSVTitle,CSVLabel,CSVISYM,CSVRotate,CSVRMesh,CSVZMesh]

					# Write to .csv
					write_to_csv(Image, CSVDir, CSVFilename, CSVHeader)

					# Write Sheath Data separately as it is not included in the VariableList format
					if image_plotsheath in ['Radial','Axial'] and i == len(VariableIndices)-1:
						np.savetxt(CSVFilename, [p for p in zip(SxAxis,Sx)], delimiter=',', fmt='%s')
						write_data_to_file(Sx, CSVDir + 'Sx.csv',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif
				#endif

				#=====#

				# Write data to ASCII files [OUTDATED FORMAT, TO BE REMOVED]
				if write_ASCII == True:
					DirWrite = create_new_folder(Dir2Dplots, '2DPlots_Data')
					write_data_to_file(Image, DirWrite + VariableStrings[i],
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					if image_plotsheath in ['Radial','Axial'] and i == len(VariableIndices)-1:
						write_data_to_file(Sx, DirWrite + 'Sx-EXT',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif
				#endif

				#=====#

			#endfor		- Variable Loop
		#endfor			- Simulation Directory Loop

		print('-------------------------------------')
		print('# 2D Steady-State Processing Complete')
		print('-------------------------------------')
	#endif

	#====================================================================#
				  #TEMPORALLY RESOLVED MOVIE_ICP.PDT FIGURES#
	#====================================================================#

	#Plot temporally resolved 2D images from movie_icp.pdt
	if savefig_movieicp2D == True:

		#for all folders being processed.
		for l in range(0,numfolders):

			#Create new folder to keep convergence variable folders in.
			DirMovieicp = create_new_folder(Dirlist[l], 'Movieicp/')

			#Enumerate variable strings in the order they appear in movie_icp.pdt
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_movieicp[l])

			# Define image extent for directory 'l'
			extent, aspectratio = DataExtent(l)

			#Extract saved iteration strings and create list for mean convergence trends
			IterArray,IterAxis = list(),list()
			for i in range(0,len(MovieIterlist[l])):
				Iter = list(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))		#List of digits within Iter
				Iter =  int(''.join(Iter))											#Join digits into single integer
				IterArray.append(Iter)												#Append to list of integers

				#Note: len(IterAxis) == len(MovieIterlist)/iterstep
				if i % iterstep == 0: IterAxis.append(Iter)							#X-axis for trend plotting
			#endfor

			#for all variables requested by the user.
			for i in tqdm(range(0,len(VariableIndices))):

				#Create new folder to keep output plots.
				DirMoviePlots = create_new_folder(DirMovieicp, '2DIterPlots_' + VariableStrings[i] + '/')

				#Create empty image array based on mesh size and symmetry options.
				try:
					numrows = len(IterMovieData[l][0][0])/R_mesh[l]
					Image = np.zeros([int(numrows),int(R_mesh[l])])
				except:
					print( 'No Iteration Data Found For '+Dirlist[l])
					break
				#endtry

				#Reshape specific part of 1D Data array into 2D image for plotting.
				for k in range(0,len(MovieIterlist[l]),iterstep):

					# Extract full 2D image for further processing.
					Image = ImageExtractor2D(IterMovieData[l][k][VariableIndices[i]],VariableStrings[i])
					Sx,SxAxis = CalcSheathExtent(folderidx=l)

					# Create figure and plot movie.icp image at iteration 'MovieIterlist[l][k]'
					fig,ax = figure(aspectratio)
					fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,VariableStrings[i],fig=fig,ax=ax)

					#=====##=====# IMAGE OVERLAYS #=====##=====#

					# Compute sheath extent and overlay onto figure if requested
					# # RM SJD: CURRENTLY USES TECPLOT2D DATA # # # # # # # # # # # # # # # # # # # # RM SJD
					Sx,SxAxis = CalcSheathExtent(folderidx=l,Orientation=image_plotsheath)
					if image_plotsheath in ['Radial','Axial']:
						PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l],Orientation=image_plotsheath)
					#endif

					#=====##=====# Image Beautification #=====##=====#

					# Define title, labels, etc...
					Title = str(MovieIterlist[l][k])
					if image_rotate == True:	Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
					elif image_rotate == False:	Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
					#endif

					# Add Colourbar (must happen before image cropping!)
					cbarlabel = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
					cax = Colourbar(ax,cbarlabel[i],image_cbarbins,Lim=CbarMinMax(ax,Image))

					# Crop image dimensions to [image_radialcrop,image_axialcrop]
					# Also resets cbar min/max to cropped region min/max
					if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
						CropImage(ax,Rotate=image_rotate)
					#endif

					# Apply mesh geometry to image
					ImageGeometry(fig,ax,image_plotsymmetry)

					# Apply image options, and enforce overides if requested
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)


					#=====##=====# Image I/O #=====##=====#

					#Save to seperate folders inside simulation folder.
					#N.B. zfill(4) Asumes Iter never exceeds 999 (i.e. max(IterArray) < 1e5)
					IterString = str(IterArray[k]).zfill(4)
					savefig(DirMoviePlots+VariableStrings[i]+'_'+IterString+ext)
					clearfigures(fig)

					#=====#

					# Write data underpinning current figure in .csv format
					if Write_CSV == True:
						CSVDir = create_new_folder(DirMoviePlots, '2DPlots_Data')
						CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
						CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
						CSVCurITER = str( MovieIterlist[l][k].replace(" ", "") )
						CSVFilename = VariableStrings[i]+'_'+IterString+'.csv'
						CSVTitle = str(Dirlist[l])
						CSVLabel = str(cbarlabel[i])
						CSVISYM = 'ISYM='+str(ISYMlist[l])
						CSVRotate = 'Rotate='+str(image_rotate)

						# Write 2D data to CSV file for current variable [i] at current iteration [k]
						CSVHeader = [CSVTitle,CSVLabel,CSVCurITER,CSVISYM,CSVRotate,CSVRMesh,CSVZMesh]
						write_to_csv(Image, CSVDir, CSVFilename, CSVHeader, mode='w')

						# Write Sheath Data separately as it is not included in the VariableList format
						if image_plotsheath in ['Radial','Axial'] and i == len(VariableIndices)-1:
							np.savetxt(CSVFilename, [p for p in zip(SxAxis,Sx)], delimiter=',', fmt='%s')
							write_data_to_file(Sx, CSVDir + 'Sx.csv',
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
					#endif

					#=====#

					#Write data to ASCII files if requested.
					if write_ASCII == True:
						DirWrite = create_new_folder(DirMovieicp, '2DPlots_Data')
						DirWriteVar = create_new_folder(DirWrite, VariableStrings[i] + '_Data')
						write_data_to_file(Image, DirWriteVar + VariableStrings[i] + '_' + str(IterArray[k]).zfill(4),
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						if image_plotsheath in ['Radial','Axial'] and k == len(VariableIndices)-1:
							write_data_to_file(Sx, DirWriteVar + 'Sx-EXT',
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
					#endif
				#endfor

				#=====#

				#Create .mp4 movie from completed images.
				Prefix = folder_name_trimmer(Dirlist[l])
				if ffmpegMovies:
					make_movie(DirMoviePlots, Prefix + '_' + VariableStrings[i])
			#endfor
		#endif

		print('-------------------------------------')
		print('# 2D Temporal Image Plotting Complete')
		print('-------------------------------------')
	#endif

	#=====================================================================#
	#=====================================================================#

	'''			# !!! WRITE 2D DATA INTO ONE LARGE CONCATENATED DATA FILE !!! #
	
					# Write data underpinning current figure in .csv format
					if Write_CSV == True:
						CSVDir = CreateNewFolder(DirMovieicp, '2DPlots_Data')
						CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
						CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
						CSVFilename = VariableStrings[i]+'.csv'
						CSVTitle = str(Dirlist[l])
						CSVLabel = str(cbarlabel[i])
						CSVISYM = 'ISYM='+str(ISYMlist[l])
						CSVRotate = 'Rotate='+str(image_rotate)
						CSVMaxITER = "IMOVIE_FRAMES="+str( MovieIterlist[l][-1].strip("ITER= ") )
						CSVCurITER = str( MovieIterlist[l][k].replace(" ", "") )
						
						print(CSVMaxCYCL)
						print(CSVCurCYCL)
						
						# Write to .csv file, including full header on first CYCLE
						if k == 0:
							CSVHeader = [CSVTitle,CSVLabel,CSVMaxITER,CSVISYM,CSVRotate,CSVRMesh,CSVZMesh]
							WriteToCSV(Image, CSVDir, CSVFilename, CSVHeader, Mode='w')
	
						# Write to .csv file, including only CYCL number for all remaining CYCLEs
						elif k > 0:
							CSVHeader = [CSVCurITER]
							WriteToCSV(Image, CSVDir, CSVFilename, CSVHeader, Mode='a')
						#endif
						
						# Write Sheath Data separately as it is not included in the VariableList format
						if image_plotsheath in ['Radial','Axial'] and i == len(variable_indices)-1:
							np.savetxt(CSVFilename, [p for p in zip(SxAxis,Sx)], delimiter=',', fmt='%s')
							WriteDataToFile(Sx, CSVDir+'Sx.csv')
						#endif
					#endif
	'''








































	#====================================================================#
					#PROFILE PLOTTING AND DATA ANALYSIS#
	#====================================================================#

	#====================================================================#
				#AXIAL AND RADIAL PROFILES FROM SINGLE FOLDERS#
	#====================================================================#

	#Generate and save lineouts of requested variables for given location.
	if savefig_monoprofiles == True:

		for l in range(0,numfolders):

			#Create VariableIndices for each folder as required.
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])

			#Generate SI scale axes for lineout plots and refresh legend.
			Raxis = GenerateAxis('Radial',ISYMlist[l])
			Zaxis = GenerateAxis('Axial',ISYMlist[l])
			Legendlist = list()

			#Generate the radial (horizontal) lineouts for a specific height.
			if len(radialprofiles) > 0:
				#Create folder to keep output plots.
				DirRlineouts = create_new_folder(Dirlist[l], '1DRadial_Profiles/')

				#Loop over all required variables and requested profile locations.
				for i in tqdm(range(0,len(VariableIndices))):

					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					# Create lists to store output for writing to file
					CSVCells = list()
					CSVProfiles = list()

					for j in range(0,len(radialprofiles)):
						#Update legend with location of each lineout.
						if len(Legendlist) < len(radialprofiles):
							Legendlist.append('Z='+str(round((radialprofiles[j])*dz[l], 2))+' cm')
						#endif

						#Plot all requested radial lines on single image per variable.
						Rlineout=ExtractRadialProfile(Data[l],VariableIndices[i],VariableStrings[i],radialprofiles[j])

						# Record profile and cell location for saving to file
						CSVCells.append(radialprofiles[j])
						CSVProfiles.append(Rlineout)

						#Plot lines for each variable at each requested slice.
						ImagePlotter1D(Raxis,Rlineout,image_aspectratio,fig,ax)

						#Write data to ASCII files if requested.
						if write_ASCII == True:
							SaveString = '_R='+str(round((radialprofiles[j])*dz[l], 2))+'cm'
							DirWrite = create_new_folder(DirRlineouts, 'Radial_Data')
							write_data_to_file([Raxis, Rlineout], DirWrite + VariableStrings[i] + SaveString,
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
					#endfor

					#Apply image options and axis labels.
					Title = 'Radial Profiles for '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel = 'Radial Distance R [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

					#Save profiles in previously created folder.
					plt.savefig(DirRlineouts+'1D_Radial_'+VariableStrings[i]+'_Profiles'+ext)
					plt.close(fig)

					# Write data underpinning current figure in .csv format
					if Write_CSV == True:
						CSVDir = create_new_folder(DirRlineouts, '1Dplots_Data')
						CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
						CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
						CSVFilename = VariableStrings[i]+'.csv'
						CSVTitle = str(Dirlist[l])
						CSVLabel = str(Ylabel)
						CSVISYM = 'ISYM='+str(ISYMlist[l])
						CSVRotate = 'Rotate='+str(image_rotate)
						CSVCells = 'Rcells='+str(CSVCells)
						CSVHeader = [CSVTitle,CSVLabel,CSVISYM,CSVRotate,CSVCells,CSVRMesh,CSVZMesh]

						# Write to .csv file
						write_to_csv(CSVProfiles, CSVDir, CSVFilename, CSVHeader)
					#endif
				#endfor
				clearfigures(fig)
			#endif

	#===================##===================#

			#Generate the vertical (height) lineouts for a given radius.
			if len(axialprofiles) > 0:
				#Create folder to keep output plots.
				DirZlineouts = create_new_folder(Dirlist[l], '1DAxial_Profiles/')
				Legendlist = list()

				#Collect and plot required data.
				for i in tqdm(range(0,len(VariableIndices))):

					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					# Create lists to store output for writing to file
					CSVCells = list()
					CSVProfiles = list()

					for j in range(0,len(axialprofiles)):
						#Perform SI conversion and save to legend.
						if len(Legendlist) < len(axialprofiles):
							Legendlist.append('R='+str(round(axialprofiles[j]*dr[l], 2))+' cm')
						#endif

						#Plot all requested radial lines on single image per variable.
						Zlineout=ExtractAxialProfile(Data[l],VariableIndices[i],VariableStrings[i],axialprofiles[j])

						# Record profile and cell location for saving to file
						CSVCells.append(axialprofiles[j])
						CSVProfiles.append(Zlineout)

						#Plot lines for each variable at each requested slice.
						ImagePlotter1D(Zaxis,Zlineout[::-1],image_aspectratio,fig,ax)

						#Write data to ASCII files if requested.
						if write_ASCII == True:
							SaveString = '_Z='+str(round((axialprofiles[j])*dr[l], 2))+'cm'
							DirWrite = create_new_folder(DirZlineouts, 'Axial_Data')
							write_data_to_file([Zaxis, Zlineout[::-1]], DirWrite + VariableStrings[i] + SaveString,
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
					#endfor

					#Apply image options and axis labels.
					Title = 'Height Profiles for '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel = 'Axial Distance Z [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

					#Save profiles in previously created folder.
					plt.savefig(DirZlineouts+'1D_Axial_'+VariableStrings[i]+'_Profiles'+ext)
					plt.close(fig)

					# Write data underpinning current figure in .csv format
					if Write_CSV == True:
						CSVDir = create_new_folder(DirZlineouts, '1Dplots_Data')
						CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
						CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
						CSVFilename = VariableStrings[i]+'.csv'
						CSVTitle = str(Dirlist[l])
						CSVLabel = str(Ylabel)
						CSVISYM = 'ISYM='+str(ISYMlist[l])
						CSVRotate = 'Rotate='+str(image_rotate)
						CSVCells = 'Zcells='+str(CSVCells)
						CSVHeader = [CSVTitle,CSVLabel,CSVISYM,CSVRotate,CSVCells,CSVRMesh,CSVZMesh]

						# Write to .csv file
						write_to_csv(CSVProfiles, CSVDir, CSVFilename, CSVHeader)
					#endif
				#endfor
				clearfigures(fig)
			#endif
		#endfor

		print('--------------------------')
		print('# Single Profiles Complete')
		print('--------------------------')
	#endif



	##====================================================================#
	#			#COMPARITIVE PROFILES FROM MULTI-FOLDERS#
	##====================================================================#

	#Plot comparitive profiles for each variable between folders.
	if savefig_compareprofiles == True:

		#Create folder to keep output plots.
		DirComparisons = create_new_folder(os.getcwd(), '/1D Comparisons')

		#Generate SI scale axes for lineout plots.
		Raxis = GenerateAxis('Radial',ISYMlist[l])
		Zaxis = GenerateAxis('Axial',ISYMlist[l])

		#Perform radial (horizontal) profile comparisons
		for j in range(0,len(radialprofiles)):

			#Create new folder for each axial or radial slice.
			ProfileFolder = 'Z='+str(round((radialprofiles[j])*dz[l], 1))+'cm'
			DirProfile = create_new_folder(DirComparisons, ProfileFolder)

			#For each requested comparison variable.
			for k in tqdm(range(0,len(Variables))):

				#Loop escape if variables that do not exist have been requested.
				if k >= 1 and k > len(VariableStrings)-1:
					break
				#endif

				#Create fig of desired size and refresh legendlist.
				fig,ax = figure(image_aspectratio,1)
				Legendlist = list()

				#For each folder in the directory.
				for l in range(0,numfolders):

					#Create VariableIndices for each folder as required.
					VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])

					#Correct VariableIndices for folders containing different icp.dat.
					VariableIndices,VariableStrings = variable_interpolator(VariableIndices, VariableStrings, MinSharedVariables, Globalnumvars)

					#Update legend with folder information.
					Legendlist.append(folder_name_trimmer(Dirlist[l]))
					Ylabels = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)

					#Plot all radial profiles for all variables in one folder.
					RProfile = ExtractRadialProfile(Data[l],VariableIndices[k],VariableStrings[k], radialprofiles[j],R_mesh[l],ISYMlist[l])

					#Plot radial profile and allow for log y-axis if requested.
					ImagePlotter1D(Raxis,RProfile,image_aspectratio,fig,ax)


					#Write data to ASCII files if requested.
					if write_ASCII == True:
						if l == 0:
							WriteFolder = 'Z='+str(round((radialprofiles[j])*dz[l], 1))+'cm_Data'
							DirWrite = create_new_folder(DirComparisons, WriteFolder)
							write_data_to_file(Raxis + ['\n'], DirWrite + VariableStrings[k], 'w',
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
						write_data_to_file(RProfile + ['\n'], DirWrite + VariableStrings[k], 'a',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif

					#Apply image options and axis labels.
					Title = 'Comparison of '+VariableStrings[k]+' Profiles at Z='+str(round((radialprofiles[j])*dz[l], 2))+'cm for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel,Legend = 'Radial Distance R [cm]',Ylabels[k],Legendlist
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
				#endfor

				#Save one image per variable with data from all simulations.
				plt.savefig(DirProfile+VariableStrings[k]+'@ Z='+str(round((radialprofiles[j])*dz[l], 1))+'cm profiles'+ext)
				clearfigures(fig)
			#endfor
		#endfor

	##===================##===================#

		#Perform vertical (axial) profile comparisons
		for j in range(0,len(axialprofiles)):

			#Create new folder for each axial or radial slice.
			ProfileFolder = 'R='+str(round((axialprofiles[j])*dr[l], 1))+'cm'
			DirProfile = create_new_folder(DirComparisons, ProfileFolder)

			#For each requested comparison variable.
			for k in tqdm(range(0,len(Variables))):

				#Loop escape if variables that do not exist have been requested.
				if k >= 1 and k > len(VariableStrings)-1:
					break
				#endif

				#Create fig of desired size and refresh legendlist.
				fig,ax = figure(image_aspectratio,1)
				Legendlist = list()

				#For each folder in the directory.
				for l in range(0,numfolders):

					#Create VariableIndices for each folder as required.
					VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])

					#Correct VariableIndices for folders containing different icp.dat.
					VariableIndices,VariableStrings = variable_interpolator(VariableIndices, VariableStrings, MinSharedVariables, Globalnumvars)

					#Update legend with folder information.
					Legendlist.append(folder_name_trimmer(Dirlist[l]))
					Ylabels = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)

					#Obtain axial profile for each folder of the current variable.
					ZProfile = ExtractAxialProfile(Data[l],VariableIndices[k],VariableStrings[k], axialprofiles[j],R_mesh[l],Z_mesh[l],ISYMlist[l])

					#Plot axial profile and allow for log y-axis if requested.
					ImagePlotter1D(Zaxis,ZProfile[::-1],image_aspectratio,fig,ax)


					#Write data to ASCII files if requested.
					if write_ASCII == True:
						if l == 0:
							WriteFolder = 'R='+str(round((axialprofiles[j])*dr[l], 1))+'cm_Data'
							DirWrite = create_new_folder(DirComparisons, WriteFolder)
							write_data_to_file(Zaxis + ['\n'], DirWrite + VariableStrings[k], 'w',
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
						write_data_to_file(ZProfile[::-1] + ['\n'], DirWrite + VariableStrings[k], 'a',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif

					#Apply image options and axis labels.
					Title = 'Comparison of '+VariableStrings[k]+' Profiles at R='+str(round((axialprofiles[j])*dr[l], 1))+'cm for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel,Legend = 'Axial Distance Z [cm]',Ylabels[k],Legendlist
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
				#endfor

				#Save one image per variable with data from all simulations.
				plt.savefig(DirProfile+VariableStrings[k]+'@ R='+str(round((axialprofiles[j])*dr[l], 2))+'cm profiles'+ext)
				clearfigures(fig)
			#endfor
		#endfor

		print('-------------------------------')
		print('# Comparitive Profiles Complete')
		print('-------------------------------')
	#endif



	##====================================================================#
	#				  #MULTI-PROFILES FROM SAME FOLDER#
	##====================================================================#

	if savefig_multiprofiles == True:

		#For each folder in turn
		for l in range(0,numfolders):
			#Create global multivar folder.
			Dirlineouts = create_new_folder(Dirlist[l], 'multivar_Profiles/')

			#Create VariableIndices for each folder as required.
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])
			multiVariableIndices,multiVariableStrings = enumerate_variables(multivar, Header_TEC2D[l])

			#Create variable labels with SI unit conversions if required.
			Ylabel = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
			multiYlabel = variable_label_maker(multiVariableStrings, Units, image_logplot, AtomicSpecies)

			#Generate the vertical (height) lineouts for a given radius.
			if len(axialprofiles) > 0:

				#Generate SI scale axes for lineout plots.
				Zaxis = GenerateAxis('Axial',ISYMlist[l])

				#Perform the plotting for all requested variables.
				for i in tqdm(range(0,len(VariableIndices))):

					#Extract the lineout data from the main data array.
					for j in range(0,len(axialprofiles)):
						#Create fig of desired size.
						fig,ax = figure(image_aspectratio,1)

						#Create folder to keep output plots.
						Slice = str(round((axialprofiles[j])*dr[l], 1))
						DirZlineouts = create_new_folder(Dirlineouts, 'R=' + Slice + 'cm/')

						#Create legendlist
						Legendlist = list()
						Legendlist.append(variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i])

						#Plot the initial variable in VariableIndices first.
						ZProfile = ExtractAxialProfile(Data[l],VariableIndices[i],VariableStrings[i], axialprofiles[j],R_mesh[l],Z_mesh[l],ISYMlist[l])
						ImagePlotter1D(Zaxis,ZProfile[::-1],image_aspectratio,fig,ax)

						#Plot all of the requested comparison variables for this plot.
						for m in range(0,len(multiVariableIndices)):
							#Plot profile for multiplot variables in compareVariableIndices.
							ZProfile = ExtractAxialProfile(Data[l],multiVariableIndices[m],multiVariableStrings[m], axialprofiles[j],R_mesh[l],Z_mesh[l],ISYMlist[l])
							ImagePlotter1D(Zaxis,ZProfile[::-1],image_aspectratio,fig,ax)

							#Update legendlist with each variable compared.
							Legendlist.append(variable_label_maker(multiVariableStrings, Units, image_logplot, AtomicSpecies)[m])
						#endfor

						#Apply image options and axis labels.
						Title = str(round((axialprofiles[j])*dr[l], 2))+'cm Height profiles for '+VariableStrings[i]+','' for \n'+Dirlist[l][2:-1]
						Xlabel,Ylabel = 'Axial Distance Z [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
						ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

						#Save figures in original folder.
						R = 'R='+str(round((axialprofiles[j])*dr[l], 2))+'_'
						plt.savefig(DirZlineouts+R+VariableStrings[i]+'_MultiProfiles'+ext)
						clearfigures(fig)
					#endfor
				#endfor
			#endif

	##===================##===================#

			#Generate the horizontal (Radial) lineouts for a given radius.
			if len(radialprofiles) > 0:
				#Create global multivar folder.
				Dirlineouts = create_new_folder(Dirlist[l], 'multivar_Profiles/')

				#Generate SI scale axes for lineout plots.
				Raxis = GenerateAxis('Radial',ISYMlist[l])

				#Perform the plotting for all requested variables.
				for i in tqdm(range(0,len(VariableIndices))):

					#Perform the plotting for all requested variables.
					for j in range(0,len(radialprofiles)):
						#Create fig of desired size.
						fig,ax = figure(image_aspectratio,1)

						#Create folder to keep output plots.
						Slice = str(round((radialprofiles[j])*dz[l], 1))
						DirRlineouts = create_new_folder(Dirlineouts, 'Z=' + Slice + 'cm/')

						#Create legendlist
						Legendlist = list()
						Legendlist.append(variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i])

						#Plot profile for initial variable in VariableIndices.
						RProfile = ExtractRadialProfile(Data[l],VariableIndices[i],VariableStrings[i], radialprofiles[j],R_mesh[l],ISYMlist[l])
						ImagePlotter1D(Raxis,RProfile,image_aspectratio,fig,ax)

						#Plot all of the requested comparison variables for this plot.
						for m in range(0,len(multiVariableIndices)):
							#Plot profile for multiplot variables in compareVariableIndices.
							RProfile = ExtractRadialProfile(Data[l],multiVariableIndices[m], multiVariableStrings[m],radialprofiles[j],R_mesh[l],ISYMlist[l])
							ImagePlotter1D(Raxis,RProfile,image_aspectratio,fig,ax)

							#Update legendlist with each variable compared.
							Legendlist.append(variable_label_maker(multiVariableStrings, Units, image_logplot, AtomicSpecies)[m])
						#endfor


						#Apply image options and axis labels.
						Title = str(round((radialprofiles[j])*dz[l], 2))+'cm Radial Profiles for '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
						Xlabel,Ylabel = 'Radial Distance R [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
						ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

						#Save lines in previously created folder.
						Z = 'Z='+str(round((radialprofiles[j])*dz[l], 2))+'_'
						plt.savefig(DirRlineouts+Z+VariableStrings[i]+'_MultiProfiles'+ext)
						clearfigures(fig)
					#endfor
				#endfor
			#endif
		#endfor

		print('-----------------------------')
		print('# Multiplot Profiles Complete')
		print('-----------------------------')
	#endif



	##====================================================================#
	#			  #ITERMOVIE 1D PROFILES - TEMPORAL ANALYSIS#			  #
	##====================================================================#

	#Plot 1D profile of requested Variables at desired locations
	if savefig_movieicp1D == True:

		for l in range(0,numfolders):

			#Create new diagnostic output folder and initiate required lists.
			DirTemporal = create_new_folder(Dirlist[l], 'Movieicp/')

			#Enumerate variable strings in the order they appear in movie_icp.pdt
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_movieicp[l])

			IterArray = list()
			#Generate a numerical iteration array for save strings
			for i in range(0,len(MovieIterlist[l])):
				Iter = list(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))		#List of digits within Iter
				Iter =  int(''.join(Iter))											#Join digits into single integer
				IterArray.append(Iter)												#Append to list of integers
			#endfor

			#Generate SI scale axes for lineout plots
			Raxis = GenerateAxis('Radial',ISYMlist[l])
			Zaxis = GenerateAxis('Axial',ISYMlist[l])

			#=====#=====#

			#Generate the radial (horizontal) lineouts for a specific height.
			if len(radialprofiles) > 0:
				#Create folder to keep output plots.
				DirRlineouts = create_new_folder(DirTemporal, '1DRadial_Profiles/')

				#Loop over all required variables and requested profile locations.
				for i in tqdm(range(0,len(VariableIndices))):
					for j in range(0,len(radialprofiles)):

						#Create new folder per radial profile
						Slice = str(round((radialprofiles[j])*dz[l], 1))
						DirProfile = create_new_folder(DirRlineouts, VariableStrings[i] + '_Z=' + Slice + 'cm')

						#Refresh lists
						Rlineoutlist = list()

						#For current profile location 'j', loop over all saved iterations
						for k in range(0,len(MovieIterlist[l]),iterstep):

							#Create fig of desired size and refresh legend
							fig,ax = figure(image_aspectratio,1)

							#Plot all requested radial lines on single image per variable.
							Rlineout=ExtractRadialProfile(IterMovieData[l][k],VariableIndices[i],VariableStrings[i],radialprofiles[j])
							#Plot lines for each variable at each requested slice.
							ImagePlotter1D(Raxis,Rlineout,image_aspectratio,fig,ax)

							#Append lineout to lineout list for overall plot
							Rlineoutlist.append(Rlineout)

							#Apply image options and axis labels.
							Title = 'Radial Profiles of '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
							Xlabel,Ylabel = 'Radial Distance R [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
							ImageOptions(fig,ax,Xlabel,Ylabel,Title,[str(MovieIterlist[l][k])],Crop=False)

							#Save profiles in variable and profile location named folder
							plt.savefig(DirProfile+'1D_Z='+Slice+'_'+VariableStrings[i]+'_'+str(IterArray[k]).zfill(4)+ext)
							plt.close(fig)

							#Write data to ASCII files if requested.
							if write_ASCII == True:
								DirWrite = create_new_folder(DirProfile, 'Radial_Data')
								SaveString1 = '_R='+str(round((radialprofiles[j])*dz[l], 2))+'cm'
								SaveString2 = '_'+str(IterArray[k]).zfill(4)
								SaveString = SaveString1+SaveString2
								write_data_to_file([Raxis, Rlineout], DirWrite + VariableStrings[i] + SaveString, 'w',
												   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
							#endif
						#endfor

						#=====#=====#

						#Create fig of desired size and refresh legend
						fig,ax = figure(image_aspectratio,1)

						#Plot all temporal profiles for VariableStrings[i] in one figure
						for k in range(0,len(Rlineoutlist)):
							ImagePlotter1D(Raxis,Rlineoutlist[k],image_aspectratio,fig,ax)
						#endfor

						#Apply image options and axis labels.
						Title = 'Radial Profiles of '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
						Xlabel,Ylabel = 'Radial Distance R [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
						ImageOptions(fig,ax,Xlabel,Ylabel,Title,[],Crop=False)

						#Save profiles Radial_Profiles folder.
						plt.savefig(DirRlineouts+'1D_'+Slice+'_'+VariableStrings[i]+'_Profiles'+ext)
						plt.close(fig)
						#endfor
					#endfor
				clearfigures(fig)
			#endif

			#===================##===================#

			#Generate the vertical (height) lineouts for a given radius.
			if len(axialprofiles) > 0:
				#Create folder to keep output plots.
				DirZlineouts = create_new_folder(DirTemporal, '1DAxial_Profiles/')

				#Collect and plot required data.
				for i in tqdm(range(0,len(VariableIndices))):
					for j in range(0,len(axialprofiles)):

						#Create new folder per axial profile
						Slice = str(round((axialprofiles[j])*dr[l], 1))
						DirProfile = create_new_folder(DirZlineouts, VariableStrings[i] + '_R=' + Slice + 'cm')

						#Refresh lists
						Zlineoutlist = list()

						#For current profile location 'j', loop over all saved iterations
						for k in range(0,len(MovieIterlist[l]),iterstep):

							#Create fig of desired size and refresh legend
							fig,ax = figure(image_aspectratio,1)

							#Plot all requested radial lines on single image per variable.
							Zlineout=ExtractAxialProfile(IterMovieData[l][k],VariableIndices[i],VariableStrings[i],axialprofiles[j])
							#Plot lines for each variable at each requested slice.
							ImagePlotter1D(Zaxis,Zlineout[::-1],image_aspectratio,fig,ax)

							#Append lineout to lineout list for overall plot
							Zlineoutlist.append(Zlineout[::-1])

							#Apply image options and axis labels.
							Title = 'Height Profiles of '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
							Xlabel,Ylabel = 'Axial Distance Z [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
							ImageOptions(fig,ax,Xlabel,Ylabel,Title,[str(MovieIterlist[l][k])],Crop=False)

							#Save profiles in variable and profile location named folder
							plt.savefig(DirProfile+'1D_R='+Slice+'_'+VariableStrings[i]+'_'+str(IterArray[k]).zfill(4)+ext)
							plt.close(fig)

							#Write data to ASCII files if requested.
							if write_ASCII == True:
								DirWrite = create_new_folder(DirProfile, 'Axial_Data')
								SaveString1 = '_R='+str(round((axialprofiles[j])*dr[l], 2))+'cm'
								SaveString2 = '_'+str(IterArray[k]).zfill(4)
								SaveString = SaveString1+SaveString2
								write_data_to_file([Zaxis, Zlineout], DirWrite + VariableStrings[i] + SaveString, 'w',
												   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
							#endif
						#endfor

						#=====#=====#

						#Create fig of desired size and refresh legend
						fig,ax = figure(image_aspectratio,1)

						#Plot all temporal profiles for VariableStrings[i] in one figure
						for k in range(0,len(Zlineoutlist)):
							ImagePlotter1D(Zaxis,Zlineoutlist[k],image_aspectratio,fig,ax)
						#endfor

						#Apply image options and axis labels.
						Title = 'Height Profiles for '+VariableStrings[i]+' for \n'+Dirlist[l][2:-1]
						Xlabel,Ylabel = 'Axial Distance Z [cm]',variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
						ImageOptions(fig,ax,Xlabel,Ylabel,Title,[],Crop=False)

						#Save profiles Axial_Profiles folder.
						plt.savefig(DirZlineouts+'1D_'+Slice+'_'+VariableStrings[i]+'_Profiles'+ext)
						plt.close(fig)
					#endfor
				#endfor
				clearfigures(fig)
			#endif
		#endfor
		print('----------------------------')
		print('# Temporal Profiles Complete')
		print('----------------------------')
	#endif



	##====================================================================#
	#			  #ITERMOVIE TRENDS - TEMPORAL ANALYSIS#
	##====================================================================#

	# Plot 1D profile of requested Variables at desired locations
	if savefig_timeaxis1D == True:

		# for all folders being processed.
		for l in range(0,numfolders):

			# Create new diagnostic output folder and initiate required lists.
			DirTemporal = create_new_folder(Dirlist[l], 'Movieicp/')

			# Create new folder and initiate required lists.
			TemporalTrends,Xaxis = list(),list()
			DirMeshAve = create_new_folder(DirTemporal, '1DTimeAxis_Profiles/')

			# Enumerate variable strings in the order they appear in movie_icp.pdt
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_movieicp[l])

			# Create list and x-axis for convergence trend plotting.
			for i in range(0,len(MovieIterlist[l])):
				IterDigits = list(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))
				IterDigits = float(''.join(IterDigits))
				IterTime = (IterDigits*DTPOS)*1000					# [ms]
				Xaxis.append(IterTime)
			#endfor

			# for all variables requested by the user.
			for i in tqdm(range(0,len(VariableIndices))):

				# Extract 2D image and take mesh averaged value for iteration trend.
				TemporalProfile = list()
				for k in range(0,len(MovieIterlist[l])):
					Image = ImageExtractor2D(IterMovieData[l][k][VariableIndices[i]],VariableStrings[i])

					# RM: DEFAULT TEMPORAL PROFILE REPRESENTS THE MESH AVERAGED VALUE
					TemporalProfile.append( sum(Image.flatten())/len(Image.flatten()) )
				#endfor

				# TemporalTrends contains all variable trends in "VariableIndices" order
				TemporalTrends.append(TemporalProfile)

				# Plot each variable against simulation real-time.
				fig, ax = plt.subplots(1, figsize=(10,10))
				ax.plot(Xaxis,TemporalProfile, lw=2)

				# Image plotting details.
				Title = 'Mesh-Averaged Temporal Profile of '+str(VariableStrings[i])+' for \n'+Dirlist[l][2:-1]
				Xlabel = 'Simulation time [ms]'
				Ylabel =variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)[i]
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

				#=====##=====# Image I/O #=====##=====#

				# Save figure.
				savefig(DirMeshAve+'Temporal_'+VariableStrings[i]+ext)
				clearfigures(fig)

				#=====#

				# Write data underpinning current figure in .csv format
				if Write_CSV == True:
					CSVDir = create_new_folder(DirMeshAve, '1DPlots_Data')
					CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
					CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
					CSVTMax = 'Sim_Time '+str( np.round(Xaxis[-1],9) )+' [ms]'
					CSVdT = 'dT '+str( np.round(Xaxis[-1]-Xaxis[-2],9) )+' [ms]'
					CSVFilename = VariableStrings[i]+'.csv'
					CSVTitle = str(Dirlist[l])
					CSVLabel = str(Ylabel)
					CSVISYM = 'ISYM='+str(ISYMlist[l])
					CSVRotate = 'Rotate='+str(image_rotate)

					# Define Header Contents
					CSVHeader = [CSVTitle,CSVLabel,CSVTMax,CSVdT]

					# Write to .csv file
					write_to_csv(TemporalProfile, CSVDir, CSVFilename, CSVHeader)
				#endif

				#=====#

				# Write data to ASCII files if requested.
				if write_ASCII == True:
					DirWrite = create_new_folder(DirMeshAve, '1DTimeAxis_Data')
					write_data_to_file(Xaxis, DirWrite + VariableStrings[i], 'w',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					write_data_to_file(['\n'] + TemporalProfile, DirWrite + VariableStrings[i], 'a',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif
			#endfor


			#=====##=====# Multi-Variable Figure #=====##=====#

			# Plot mesh averaged value over 'real-time' in simulation.
			Legend = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
			fig, ax = plt.subplots(1, figsize=(14,10))

			# Plot each variable in ConvergenceTrends to single figure.
			for i in range(0,len(TemporalTrends)):
				TemporalTrends[i] = Normalise(TemporalTrends[i])[0]
				ax.plot(Xaxis,TemporalTrends[i], lw=2)
			#endfor

			# Image plotting details.
			Title = 'Mesh-Averaged Temporal Profiles of '+str(VariableStrings)+' for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Simulation time [ms]','Normalised Mesh-Average Value'
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
			ax.set_ylim(0,1)

			# Save figure.
			savefig(DirMeshAve + folder_name_trimmer(Dirlist[l]) + '_Normalised' + ext)
			clearfigures(fig)
		#endfor

		print('--------------------------')
		print('# Temporal Trends Complete')
		print('--------------------------')
	#endif



	##====================================================================#
	#			  #ITERMOVIE TRENDS - CONVERGENCE ANALYSIS#
	##====================================================================#

	#Perform convergence checking from movie_icp.pdt.
	if savefig_convergence == True:

		#for all folders being processed.
		for l in range(0,numfolders):

			#Create new folder to keep convergence variable folders in.
			DirConvergence = create_new_folder(Dirlist[l], 'Movieicp/')

			#Enumerate variable strings in the order they appear in movie_icp.pdt
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_movieicp[l])

			#Extract saved iteration strings and create list for mean convergence trends
			ConvergenceTrends,IterArray,IterAxis = list(),list(),list()
			for i in range(0,len(MovieIterlist[l])):
				Iter = list(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))		#List of digits within Iter
				Iter =  int(''.join(Iter))											#Join digits into single integer
				IterArray.append(Iter)												#Append to list of integers

				#Note: len(IterAxis) == len(MovieIterlist)/iterstep
				if i % iterstep == 0: IterAxis.append(Iter)							#X-axis for trend plotting
			#endfor

			#for all variables requested by the user.
			for i in tqdm(range(0,len(VariableIndices))):

				#Append new list to convergenceTrends for each variable.
				ConvergenceTrends.append(list())

				#for each iteration
				for k in range(0,len(MovieIterlist[l]),iterstep):

					#Extract full 2D image for further processing.
					Image = ImageExtractor2D(IterMovieData[l][k][VariableIndices[i]],VariableStrings[i])

					#If probeloc is provided, take convergence trend from that location
					#Else use the MESH AVERAGED VALUE for general convergence trend
					if len(probeloc) == 2:
						ConvergenceVal = Image[probeloc[0],probeloc[1]]
						labelstring = 'Normalised Values @ R='+str(round(probeloc[0]*dr[l],2))+', Z='+str(round(probeloc[1]*dz[l],2))+' cm'
					else:
						ConvergenceVal = sum(Image.flatten())/len(Image.flatten())
						labelstring = 'Normalised Mesh-Averaged Values'
					#endif

					ConvergenceTrends[i] = np.append(ConvergenceTrends[i], ConvergenceVal)
				#endfor
			#endfor

			# Exit if less than two trend variables have been extracted
			Continue = True
			if len(ConvergenceTrends[0]) < 2:
				print(' ')
				print("#====================================================================#")
				print("WARNING: Convergence Trends Contains Too Few Iterations For Comparison")
				print("                          Aborting Diagnostic                         ")
				print("#====================================================================#")
				Continue = False
			#endif

			#=================#

			if Continue == True:

				#Plot a convergence check for all variables in each folder.
				Legend = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
				fig, ax = plt.subplots(1, figsize=(14,10))

				#Normalise and plot each variable in ConvergenceTrends to single figure.
				for i in range(0,len(ConvergenceTrends)):
					NormFactor=max(max(ConvergenceTrends[i]),abs(min(ConvergenceTrends[i])))
					ConvergenceTrends[i] = Normalise(ConvergenceTrends[i],NormFactor)[0]
					ax.plot(IterAxis,ConvergenceTrends[i], lw=2)
				#endfor

				#Various debug outputs
				if IDEBUG == True:
					print(' ')
					print('Convergence Limits:',Limits)
					print(' ')
				#endif

				#Image plotting details.
				Title = 'Convergence of '+str(VariableStrings)+' for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Simulation Iteration [-]',labelstring
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
				plt.tight_layout()


				# Compute and display convergence trends to terminal
				print('')
				print('Percentage Variation At Final Iteration:')
				print('Var \t Relative % \t Abs Delta \t Abs Sum')
				for i in range(0,len(ConvergenceTrends)):

					# The mesh averaged value 	# NOTE, NORMALISED TO MAX VALUE FOR VARIABLE [i] #
					MeshAveragedVal = abs(ConvergenceTrends[i][-1])

					# Compute absolute delta and absolute sum
					AbsoluteDelta = abs(ConvergenceTrends[i][-1] - ConvergenceTrends[i][-2])
					AbsoluteSum = abs(ConvergenceTrends[i][-1]) + abs(ConvergenceTrends[i][-2])

					# Compute symmetric relative change (2*deltaValue / (Value[n] + Value[n-1]))
					RelativeChange = (2*AbsoluteDelta / AbsoluteSum )
					RelativePercent = round(RelativeChange*100.0, 4)

					# Print convergence data to terminal
					print(
					VariableStrings[i], '\t',
					RelativePercent, '%', '\t',
					'Δ=', round(AbsoluteDelta, 5), '\t',
					'Σ=', round(AbsoluteSum, 5)
					)
				#endfor

				#Save figure.
				savefig(DirConvergence + folder_name_trimmer(Dirlist[l]) + '_Convergence' + ext)
				clearfigures(fig)
			#endfor
		#endif

		print('-------------------------------')
		print('# Convergence Checking Complete')
		print('-------------------------------')
	#endif

	#=====================================================================#
	#=====================================================================#












































	#====================================================================#
					 #GENERAL ENERGY DISTRIBUTION ANALYSIS#
	#====================================================================#

	#====================================================================#
					#ION-NEUTRAL ANGULAR ENERGY DISTRIBUTIONS#
	#====================================================================#

	if savefig_IEDFangular == True:

		#For all simulation folders.
		for l in range(0,numfolders):

			#Create new folder for keeping EEDF/IEDF if required.
			DirIEDF = create_new_folder(Dirlist[l], 'EDFplots')

			#Create VariableIndices for requested EDF species and extract images.
			VariableIndices,VariableStrings = enumerate_variables(IEDFVariables, Header_IEDF[l])

			#For all requested variables.
			for i in tqdm(range(0,len(VariableIndices))):

				#Create any required arrays
				EDFprofile = list()

				#Extract image from required variable and create required profile lists.
				#Flatten angular distribution across all angles to produce energy distribution.
				EDFImage = ImageExtractor2D(DataIEDF[l][VariableIndices[i]],Rmesh=int(EDFangle),Zmesh=int(EDFbins))
				for j in range(0,len(EDFImage)): EDFprofile.append(sum(EDFImage[j]))

				#Smooth kinetic data prior to analysis if requested (Savitzk-Golay filter)
				if PlotKineticFiltering == True:
					WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
					EDFprofile = (savgol_filter(EDFprofile, WindowSize, PolyOrder)).tolist()
					#endfor
				#endif

				#Obtain conversion from energy-bin axis to eV axis and construct energy axis
				deV = (EMAXIPCMC/IEBINSPCMC)
				eVaxis = list()
				for j in range (0,int(IEBINSPCMC)):
					eVaxis.append(j*deV)
				#endfor

				#Transpose Image for plotting and reverse both lists to align with other data.
				EDFImage, EDFprofile = EDFImage[::-1].transpose(), EDFprofile[::-1]

				#Determine region of IEDF to plot based on threshold value from array maximum.
				Threshold = EDF_Threshold*max(EDFprofile)
				index_max = np.argmax(EDFprofile)
				for j in range(int(index_max),len(EDFprofile)):
					if EDFprofile[j] < Threshold and j != 0:
						eVlimit = j*deV
				#break
					elif j == len(EDFprofile)-1:
						eVlimit = EMAXIPCMC
					#endif
				#endfor


				#Plot the angular distribution and EDF of the required species.
				fig,ax = figure([11,9], 2, shareX=True)

				Title = Dirlist[l][2::]+'\n'+VariableStrings[i]+' Angular Energy Distribution Function'
				Extent=[0,int(EMAXIPCMC), -len(EDFImage)/2,len(EDFImage)/2]

				#Angularly resolved IEDF Figure
				im = ax[0].imshow(EDFImage, extent=Extent, aspect='auto')
				cax = Colourbar(ax[0],VariableStrings[i]+' EDF($\\theta$)',5)
				Xlabel,Ylabel = '','Angular Dispersion [$\\theta^{\circ}$]'
				ImageCrop = [[0,int(eVlimit)],[-45,45]]					#[[X1,X2],[Y1,Y2]]
				ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,Crop=ImageCrop,Rotate=False)

				#Angle Integrated IEDF figure
				ax[1].plot(eVaxis,EDFprofile, lw=2)
				cax = InvisibleColourbar(ax[1])
				Xlabel,Ylabel = 'Energy [eV]',VariableStrings[i]+' EDF \n [$\\theta$ Integrated]'
				ImageCrop = [[0,int(eVlimit)],[0,max(EDFprofile)*1.05]]	#[[X1,X2],[Y1,Y2]]
				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=ImageCrop,Rotate=False)

				#=====##=====# Image I/O #=====##=====#

				# Save Figure
				plt.savefig(DirIEDF+VariableStrings[i].replace(' ','')+'_EDF'+ext)
				clearfigures(fig)

				#=====#

				# Write data underpinning current figure in .csv format
				if Write_CSV == True:
					CSVDir = create_new_folder(DirIEDF, 'EDF_Data')
					CSVEnergyAxis = 'Energy_Range_[eV] '+str([eVaxis[0],eVaxis[-1]])+'  :: Bin_Size_[eV] '+str(deV)
					CSVAngleAxis =  'Angle_Range_[Deg] '+str([-45,45])+'  :: Bin_Size_[Deg] '+str(45./len(EDFImage))
					CSVFilename = VariableStrings[i].replace(' ','')+'.csv'
					CSVTitle = str(DirIEDF)
					CSVLabel = str(VariableStrings[i])
					CSVHeader = [CSVTitle,CSVLabel,CSVEnergyAxis,CSVAngleAxis]

					# Write to .csv
					write_to_csv(EDFImage, CSVDir, CSVFilename, CSVHeader)
				#endif

				#=====#

				# Write data underpinning current figure in ASCII format		- OUTDATED, TO REMOVE
				if write_ASCII == True:
					if i == 0:
						DirASCII = create_new_folder(DirIEDF, 'EDF_Data')
						write_data_to_file(eVaxis, DirASCII + VariableStrings[i], 'w',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif
					write_data_to_file(EDFImage, DirASCII + VariableStrings[i] + '_IEDFAngular', 'w',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					write_data_to_file(EDFprofile, DirASCII + VariableStrings[i] + '_IEDFProfile', 'w',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif
			#endfor
		#endfor
	#endif



	#====================================================================#
					#ION-NEUTRAL ANGULAR ENERGY ANALYSIS#
	#====================================================================#

	if savefig_IEDFtrends == True:

		#For all requested IEDF variables
		for i in tqdm(range(0,len(IEDFVariables))):
			#Initiate figure for current variable and any required lists.
			Legendlist,EDFprofiles,EDFEnergyProfile = list(),list(),list()
			Mode_eV,Mean_eV,Median_eV = list(),list(),list()
			Total_eV,Range_eV,FWHM_eV = list(),list(),list()
			GlobRange_eV = [0,0]
			fig,ax = figure()

			#Create new global trend folder if it doesn't exist already.
			TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
			TrendVariable = ''.join(TrendVariable)												#Single string of chars
			DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')
			#Create folder for IEDF trends if it doesn't exist already.
			DirIEDFTrends = create_new_folder(DirTrends, 'IEDF Trends')

			#For all simulation folders.
			for l in range(0,numfolders):
				EDFprofile = list()

				#Create VariableIndices for requested EDF species and extract images.
				VariableIndices,VariableStrings = enumerate_variables(IEDFVariables, Header_IEDF[l])
				Legendlist.append(folder_name_trimmer(Dirlist[l]))

				#Extract image from required variable and flatten angular distribution profile.
				Image = ImageExtractor2D(DataIEDF[l][VariableIndices[i]],Rmesh=EDFangle,Zmesh=EDFbins)
				for j in range(0,len(Image)): EDFprofile.append(sum(Image[j]))
				#Transpose Image for plotting and reverse both lists due to reading error.
				Image, EDFprofile = Image[::-1].transpose(), EDFprofile[::-1]

				#Obtain conversion from energy-bin axis to eV axis and construct energy axis
				deV, eVaxis = (EMAXIPCMC/IEBINSPCMC), list()
				for j in range (0,int(IEBINSPCMC)): eVaxis.append(j*deV)

				#Plot 1D EDF variable profile on open figure for each simulation folder.
				ax.plot(eVaxis,EDFprofile, lw=2)

				#==========#
				#==========#

				#Perform energy trend diagnostics on current folder (variable 'i') IEDF
				#Smooth kinetic data prior to analysis if requested (Savitzk-Golay filter)
				if KineticFiltering == True:
					WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
					EDFprofile = (savgol_filter(EDFprofile, WindowSize, PolyOrder)).tolist()
					#endfor
				#endif

				#Energy extrema analysis: Returns maximum energy where the fraction is above threshold.
				Threshold = EDF_Threshold*max(EDFprofile)
				ThresholdArray = list()
				for j in range(0,len(EDFprofile)):
					if EDFprofile[j] >= Threshold: ThresholdArray.append(j)
				#endfor
				IndexRange = [min(ThresholdArray),max(ThresholdArray)]
				Range_eV = [IndexRange[0]*deV,IndexRange[1]*deV]
				#Save global extrema between folders for plotting range.
				if Range_eV[0] < GlobRange_eV[0]: GlobRange_eV[0] = Range_eV[0]
				if Range_eV[1] > GlobRange_eV[1]: GlobRange_eV[1] = Range_eV[1]

				#Total energy analysis: Returns total energy contained within EDF profile
				EDFEnergyProfile = list()
				for j in range(0,len(EDFprofile)):
					EDFEnergyProfile.append( EDFprofile[j]*(j*deV) )
				#endfor
				Total_eV.append(sum(EDFEnergyProfile))

				#Average energy analysis: Returns mean/mode/median energies from IEDF.
				#Modal energy taken as most common energy fraction (disregarding EDF_threshold)
				Mode_eV.append( EDFprofile.index(max(EDFprofile))*deV )

				#Mean energy obtained as integrated energy fraction most closely matching total energy,
				#which is averaged over effective energy range determined from EDF_threshold.
				if GlobMeanCalculation == 'MeanEnergy':
					BinAveragedEnergy = sum(EDFEnergyProfile)/len(EDFEnergyProfile)		#(Total Energy/Num Bins)
					ResidualArray = list()
					#Calculate Residuals using EDF_threshold as upper energy percentile
					for j in range(IndexRange[0],IndexRange[1]):
						ResidualArray.append(abs(EDFEnergyProfile[j]-BinAveragedEnergy))
					#endfor
					#Capture mean/residual intersections, sort high energy intersection indices first
					NumIntersections = int(len(EDFprofile)/10)
					Intersections = np.argsort(ResidualArray)[:NumIntersections]
					Intersections = sorted(Intersections,reverse=False)
	#				#Single Intersection case - special case, maintained mostly for debugging purposes
	#				Intersection = (np.abs(EDFEnergyProfile-BinAveragedEnergy)).argmin()

					#Index k defines the IEDF array lower energy edge corresponding to the threshold
					IntersectionRegions,k = list(),0
					for j in range(0,len(Intersections)-1):
						RegionThreshold = 0.05*len(EDFEnergyProfile)
						#If threshold jump observed, save current intersect region index (k)
						if Intersections[j]+RegionThreshold < Intersections[j+1]:
							IntersectionRegions.append(Intersections[k:j])
							k = j+1
						#Save final intersect region index
						elif j == len(Intersections)-2:
							IntersectionRegions.append(Intersections[k:j])
						#endif
					#endfor
					Intersections = list()
					#If odd number of intersections, likely that low energy one was missed
					if len(IntersectionRegions) % 2 == 1: Intersections.append(1.0)
					#For all intersection regions identify the central index
					for j in range(0,len(IntersectionRegions)):
						try:
							RegionCentralIndex = int(len(IntersectionRegions[j])/2.0)
							Intersections.append( IntersectionRegions[j][RegionCentralIndex] )
						except:
							#If no intersections, assume a single zero energy intersection
							if j == 0: Intersections.append(0.0)
						#endtry
					#endfor
					#Extrema represent FWHM of EDF, mean energy lies between extrema
					FWHM_eV.append([min(Intersections)*deV,max(Intersections)*deV])
					MeanEnergyIndex = (max(Intersections)+min(Intersections))/2
				#endif

				#Mean energy obtained as ion energy with population fraction most closely matching IEDF average
				#OLD MEAN ENERGY DEFINITION - MEAN DEFINED BY FRACTION NOT BY ENERGY
				if GlobMeanCalculation == 'MeanFraction':
					BinAveragedFraction = sum(EDFprofile)/len(EDFprofile)
					ResidualArray = list()
					#Calculate Residuals using EDF_threshold as upper energy percentile
					for j in range(IndexRange[0],IndexRange[1]):
						ResidualArray.append(abs(EDFprofile[j]-BinAveragedFraction))
					#endfor
					#Capture mean/residual intersections, sort high energy intersection indices first
					NumIntersections = int(len(EDFprofile)/8)
					Intersections = np.argsort(ResidualArray)[:NumIntersections]
					Intersections = sorted(Intersections,reverse=False)

					#Index k defines the IEDF array lower energy edge corresponding to the threshold
					IntersectionRegions,k = list(),0
					for j in range(0,len(Intersections)-1):
						RegionThreshold = 0.05*len(EDFprofile)
						#If threshold jump observed, save current intersect region index (k)
						if Intersections[j]+RegionThreshold < Intersections[j+1]:
							IntersectionRegions.append(Intersections[k:j])
							k = j+1
						#Save final intersect region index
						elif j == len(Intersections)-2:
							IntersectionRegions.append(Intersections[k:j])
						#endif
					#endfor

					Intersections = list()
					#If odd number of intersections, likely that low energy one was missed
					if len(IntersectionRegions) % 2 == 1: Intersections.append(1.0)
					#For all intersection regions identify the central index
					for j in range(0,len(IntersectionRegions)):
						try:
							RegionCentralIndex = int(len(IntersectionRegions[j])/2.0)
							Intersections.append( IntersectionRegions[j][RegionCentralIndex] )
						except:
							#If no intersections, assume a single zero energy intersection
							if j == 0: Intersections.append(0.0)
						#endtry
					#endfor
					#Extrema represent FWHM of EDF, mean energy lies between extrema
					FWHM_eV.append([min(Intersections)*deV,max(Intersections)*deV])
					MeanEnergyIndex = (max(Intersections)+min(Intersections))/2
				#endif
				Mean_eV.append( MeanEnergyIndex*deV )

	#			#Median energy calculated as EDF index representing midpoint of equal integrated energies
	#			RisingSum,FallingSum,AbsDiff = 0.0,0.0,list()
	#			for j in range(0,len(EDFprofile)):
	#				Rising_j, Falling_j = j, (len(EDFprofile)-1-2)-j
	#				RisingSum += EDFprofile[Rising_j]*(Rising_j*deV)
	#				FallingSum += EDFprofile[Falling_j]*(Falling_j*deV)
	#				AbsDiff.append(abs(RisingSum-FallingSum))
	#				#endif
	#			#endfor
	#			MedianIndex = AbsDiff.index(min(AbsDiff))
	#			Median_eV.append( MedianIndex*deV ) 				#### NB: MEDIANS' ALL FUCKED UP BRAH! ####

	#			#Particle energy variance analysis: Returns FWHM of energy distribution.
	#			Take mean and draw line at y = mean
	#			Calculate where y = mean intercepts EDFprofile
	#			If only one intercept, first intercept is x = 0
	#			Integrate EDFprofile indices between intercepts
	#			Save in 1D array, can be used to get energy spread percentage.

				#==========#
				#==========#

				if IDEBUG == True:
					fig2,ax2 = figure()
					try: BinAveragedValue = BinAveragedEnergy		#MeanEnergy Value
					except: BinAveragedValue = BinAveragedFraction	#MeanFraction Value
					Ylims = [0,max(EDFEnergyProfile)]

					fig2,ax2 = figure()
					ax2.plot(EDFEnergyProfile, 'k-', lw=2)
					ax2.plot(ResidualArray, 'r-', lw=2)
					ax2.plot((0,len(EDFEnergyProfile)),(BinAveragedValue,BinAveragedValue),'b--',lw=2)
					ax2.plot((max(Intersections),max(Intersections)),(Ylims[0],Ylims[1]),'b--',lw=2)
					ax2.plot((min(Intersections),min(Intersections)),(Ylims[0],Ylims[1]),'b--',lw=2)
					ax2.plot((MeanEnergyIndex,MeanEnergyIndex),(Ylims[0],Ylims[1]),'m--',lw=2)
					ax2.legend(['Integrated Energy','abs(ResidualArray)','BinAveragedEnergy/Fraction'])
	#				plt.savefig(DirIEDFTrends+'_DebugOutput.png')
				#endif
			#endfor

			#Write data to ASCII format datafile if requested.
			if write_ASCII == True:
				if i == 0:
					DirASCII = create_new_folder(DirTrends, 'Trend_Data')
					DirASCIIIEDF = create_new_folder(DirASCII, 'IEDF_Data')
				#endif
				write_data_to_file(Legendlist + ['\n'] + Total_eV, DirASCIIIEDF + VariableStrings[i] + '_Total', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(Legendlist + ['\n'] + Range_eV, DirASCIIIEDF + VariableStrings[i] + '_Range', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(Legendlist + ['\n'] + Mode_eV, DirASCIIIEDF + VariableStrings[i] + '_Mode', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(Legendlist + ['\n'] + Mean_eV, DirASCIIIEDF + VariableStrings[i] + '_Mean', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif

			##IEDF PROFILES##
			#===============#

			#Apply image options to IEDF plot generated in the above loop.
			Title = Dirlist[l][2::]+'\n'+VariableStrings[i]+' Angular Energy Distribution Function Profiles'
			Xlabel,Ylabel = 'Energy [eV]',VariableStrings[i]+' EDF [$\\theta$ Integrated]'
			ImageCrop = [ [0,GlobRange_eV[1]], [] ]		#[[X1,X2],[Y1,Y2]]
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=ImageCrop,Rotate=False)

			plt.savefig(DirIEDFTrends+VariableStrings[i]+'_EDFprofiles'+ext)
			clearfigures(fig)

			##ENERGY ANALYSIS##
			#=================#

			#Plot IEDF average energies with respect to simulation folder names.
			fig,ax = figure()
			TrendPlotter(ax,Mean_eV,Legendlist,NormFactor=0)
			TrendPlotter(ax,Mode_eV,Legendlist,NormFactor=0)
	##		TrendPlotter(ax,Range_eV[0],Legendlist,NormFactor=0)
	##		TrendPlotter(ax,Range_eV[1],Legendlist,NormFactor=0)

			Title = Dirlist[l][2::]+'\n'+'Average '+VariableStrings[i]+' Energies'
			Legend = ['EDF Mean Energy','EDF Mode Energy','EDF Min Energy','EDF Max Energy']
			Xlabel,Ylabel = 'Varied Property',VariableStrings[i]+' Energy [eV]'
			ImageCrop = [[],[0,max(Mean_eV+Mode_eV)*1.15]]		#[[X1,X2],[Y1,Y2]]
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=ImageCrop,Rotate=False)

			plt.savefig(DirIEDFTrends+VariableStrings[i]+'_AverageEnergies'+ext)
			clearfigures(fig)
		#endfor
	#endif


	if any([savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]) == True:
		print('--------------------------------')
		print('# EEDF/IEDF Processing Complete.')
		print('--------------------------------')
	#endif

	#=====================================================================#
	#=====================================================================#










































	#====================================================================#
					 #GENERAL TREND PLOTTING ANALYSIS#
	#====================================================================#

	#====================================================================#
					#COMPARATIVE TRENDS -- MULTI-FOLDER#
	#====================================================================#

	if savefig_trendphaseaveraged == True or print_generaltrends == True:

		#Create trend folder for outputs - Assumes that scanned variable is within trimmed foler name
		TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

		#For each requested comparison variable.
		for k in tqdm(range(0,len(Variables))):

			#Create VariableIndices for largest output, only compare variables shared between all folders.
			VariableIndices,VariableStrings = enumerate_variables(Variables, Header_TEC2D[l])
			VariableIndices,VariableStrings = variable_interpolator(VariableIndices, VariableStrings, MinSharedVariables, Globalnumvars)

			#Create Y-axis legend for each variable to be plotted.
			YaxisLegend = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)

			#Loop escape if variables that do not exist have been requested.
			if k >= 1 and k > len(VariableStrings)-1:
				break
			#endif

			#Create fig of desired size and refresh legendlist.
			fig,ax = figure(image_aspectratio,1)
			Legendlist = list()


			##AXIAL TRENDS##
			#===============#

			#Perform trend analysis on requested axial profiles.
			for j in range(0,len(axialprofiles)):

				#Create folder for axial trends if needed.
				DirAxialTrends = create_new_folder(DirTrends, 'Axial Trends')

				#Take Trend at Given Location or Default to Min/Max Trends.
				if len(probeloc) == 2:
					#Append requested position to the legendlist.
					R,Z = probeloc[0],probeloc[1]
					Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
					Legendlist.append(Location)
					#Take trend at given location if specified.
					Xaxis,Trend = TrendAtGivenLocation([R,Z],VariableIndices[k],VariableStrings[k])

				elif len(probeloc) == 1:
					#Append requested position to the legendlist.
					R,Z = probeloc[0],axialprofiles[j]
					Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
					Legendlist.append(Location)
					#Take trend at given location if specified.
					Xaxis,Trend = TrendAtGivenLocation([R,Z],VariableIndices[k],VariableStrings[k])

				else:
					#Obtain min/max trend values for requested profile over all folders.
					Xaxis,MaxTrend,MinTrend = MinMaxTrends(axialprofiles[j],'Axial',k)
					Trend = MaxTrend
					#Append the radial position to the legendlist.
					Legendlist.append( 'R='+str(round((axialprofiles[j]*dr[l]), 2))+'cm' )
				#endif

				#Plot trends for each variable over all folders, applying image options.
				TrendPlotter(ax,Trend,Xaxis,NormFactor=0)
				Title='Trend in max '+VariableStrings[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Varied Property','Max '+YaxisLegend[k]
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

				#Write data to ASCII format datafile if requested.
				if write_ASCII == True:
					if j == 0:
						DirASCII = create_new_folder(DirTrends, 'Trend_Data')
						DirASCIIAxial = create_new_folder(DirASCII, 'Axial_Data')
						write_data_to_file(Xaxis + ['\n'], DirASCIIAxial + VariableStrings[k] + '_Trends', 'w',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif
					write_data_to_file(Trend + ['\n'], DirASCIIAxial + VariableStrings[k] + '_Trends', 'a',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif

			#Save one image per variable with data from all simulations.
			if len(axialprofiles) > 0:
				plt.savefig(DirAxialTrends+'Axial Trends in '+VariableStrings[k]+ext)
				clearfigures(fig)
			#endif


			##RADIAL TRENDS##
			#===============#

			#Create fig of desired size and refresh legendlist.
			fig,ax = figure(image_aspectratio,1)
			Legendlist = list()

			#Perform trend analysis on requested radial profiles.
			for j in range(0,len(radialprofiles)):

				#Create folder for axial trends if needed.
				DirRadialTrends = create_new_folder(DirTrends, 'Radial Trends')

				#Take Trend at Given Location or Default to Min/Max Trends.
				if len(probeloc) == 2:
					#Append requested position to the legendlist.
					R,Z = probeloc[0],probeloc[1]
					Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
					Legendlist.append(Location)
					#Take trend at given location if specified.
					Xaxis,Trend = TrendAtGivenLocation([R,Z],VariableIndices[k],VariableStrings[k])

				elif len(probeloc) == 1:
					#Append requested position to the legendlist.
					R,Z = radialprofiles[j],probeloc[0],
					Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
					Legendlist.append(Location)
					#Take trend at given location if specified.
					Xaxis,Trend = TrendAtGivenLocation([R,Z],VariableIndices[k],VariableStrings[k])

				else:
					#Obtain min/max trend values for requested profile over all folders.
					Xaxis,MaxTrend,MinTrend = MinMaxTrends(radialprofiles[j],'Radial',k)
					Trend = MaxTrend
					#Append the axial position to the legendlist.
					Legendlist.append( 'Z='+str(round((radialprofiles[j]*dz[l]), 2))+'cm' )
				#endif

				#Plot trends for each variable over all folders, applying image options.
				TrendPlotter(ax,Trend,Xaxis,NormFactor=0)
				Title='Trend in max '+VariableStrings[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Varied Property','Max '+YaxisLegend[k]
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

				#Write data to ASCII format datafile if requested.
				if write_ASCII == True:
					if j == 0:
						DirASCII = create_new_folder(DirTrends, 'Trend_Data')
						DirASCIIRadial = create_new_folder(DirASCII, 'Radial_Data')
						write_data_to_file(Xaxis + ['\n'], DirASCIIRadial + VariableStrings[k] + '_Trends', 'w',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif
					write_data_to_file(Trend + ['\n'], DirASCIIRadial + VariableStrings[k] + '_Trends', 'a',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif

			#Save one image per variable with data from all simulations.
			if len(radialprofiles) > 0:
				plt.savefig(DirRadialTrends+'Radial Trends in '+VariableStrings[k]+ext)
				clearfigures(fig)
			#endif
		#endfor
	#endif




	#=====================================================================#
							# DC-BIAS CALCULATOR #
	#=====================================================================#

	if savefig_trendphaseaveraged == True or print_DCbias == True:

		#Create Trend folder to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

		#Initiate lists required for storing data.
		Xaxis = list()
		DCbias = list()

		#For all folders.
		for l in range(0,numfolders):

			#Create VariableIndices for each folder as required.
			Process,Variable = enumerate_variables(['P-POT'], Header_TEC2D[l])

			#Update X-axis with folder information.
			Xaxis.append(folder_name_trimmer(Dirlist[l]))

			#Locate powered electrode for bias extraction.
			Rlineoutloc = WaveformLoc(electrodeloc,'2D')[0]
			Zlineoutloc = WaveformLoc(electrodeloc,'2D')[1]

			#Obtain radial and axial profiles for further processing.
			try: Rlineout = ExtractRadialProfile(Data[l],Process[0],Variable[0],Rlineoutloc,R_mesh[l],  ISYMlist[l])
			except: Rlineout = float('NaN')
			#endtry
			try: Zlineout = ExtractAxialProfile(Data[l],Process[0],Variable[0],Zlineoutloc,R_mesh[l],Z_mesh[l],ISYMlist[l])
			except: Zlineout = float('NaN')
			#endtry

			#Obtain DCbias on axis and across the centre radius of the mesh.
			AxialDCbias = DCbiasMagnitude(Zlineout[::-1])
			RadialDCbias = DCbiasMagnitude(Rlineout)

			#Choose axial or radial DCbias based on user input, else autoselect most probable.
			if DCbiasaxis == 'Radial':
				DCbias.append(RadialDCbias)
			elif DCbiasaxis == 'Axial':
				DCbias.append(AxialDCbias)
			elif DCbiasaxis == 'Auto':
				#Compare Axial and Radial DCbias, if same pick Axial, if not pick the largest.
				if AxialDCbias != RadialDCbias:
					if abs(AxialDCbias) > abs(RadialDCbias):
						DCbias.append(AxialDCbias)
					else:
						DCbias.append(RadialDCbias)
					#endif
				else:
					DCbias.append(AxialDCbias)
				#endif
			#endif

			#Display DCbias to terminal if requested.
			if print_DCbias == True:
				print(Dirlist[l])
				print('DC Bias:',round(DCbias[l],5),'V')
			#endif
		#endfor

		#Write data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = create_new_folder(DirTrends, 'Trend_Data')
			DCASCII = [Xaxis,DCbias]
			write_data_to_file(DCASCII, DirASCII + 'DCbias_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
		#endif

		#Plot and beautify the DCbias, applying normalization if requested.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(ax,DCbias,Xaxis,NormFactor=0)

		#Apply image options and axis labels.
		Title = 'Trend in DCbias with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','DC bias [V]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirTrends+'Powered Electrode DCbias'+ext)
		clearfigures(fig)
	#endif



	#====================================================================#
						#POWER DEPOSITED DIAGNOSTIC#
	#====================================================================#

	if savefig_trendphaseaveraged == True or print_totalpower == True:

		#Create Trend folder to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

		#Create required lists.
		RequestedPowers,DepositedPowerList = list(),list()
		Xaxis,Powers = list(),list()

		#Update X-axis with folder information.
		for l in range(0,numfolders): Xaxis.append(folder_name_trimmer(Dirlist[l]))

		#Create list of requested power variables and ensure they also exist in all compared folders
		for i in range(0,len(Variables)):
			if 'POW-' in Variables[i] and Variables[i] in MinSharedVariables:
				RequestedPowers.append(Variables[i])
			#endif
		#endfor

		#For each different power deposition mechanism requested.
		for k in range(0,len(RequestedPowers)):
			#For all folders.
			for l in range(0,numfolders):

				#Create extract data for the neutral flux and neutral velocity.
				VariableIndices,VariableStrings = enumerate_variables(RequestedPowers, Header_TEC2D[l])

				#Extract 2D power density data array [W/m3]
				PowerDensity = ImageExtractor2D(Data[l][VariableIndices[k]])

				#=====#=====#

				#Cylindrical integration of power per unit volume ==> total coupled power.
				if IXZlist[l] == 0:

					Power = 0												#[W]
					#For each axial slice
					for j in range(0,Z_mesh[l]):
						#For each radial slice
						for i in range(0,R_mesh[l]-1):
							#Calculate radial plane volume of a ring at radius [i], correcting for central r=0.
							InnerArea = np.pi*( (i*(dr[l]/100))**2 )		#[m^2]
							OuterArea = np.pi*( ((i+1)*(dr[l]/100))**2 )	#[m^2]
							RingVolume = (OuterArea-InnerArea)*(dz[l]/100)	#[m^3]

							#Calculate Power by multiplying power density for ring[i] by volume of ring[i]
							Power += PowerDensity[j][i]*RingVolume 			#[W]
						#endfor
					#endfor
					DepositedPowerList.append(Power)

					#Display power to terminal if requested.
					if print_totalpower == True:
						print(Dirlist[l])
						print(RequestedPowers[k]+' Deposited:',round(Power,4),'W')
					#endif

				#=====#=====#

				#Cartesian integration of power per unit volume ==> total coupled power.
				elif IXZlist[l] == 1:

					Power = 0												#[W]
					#For each axial slice
					for j in range(0,Z_mesh[l]):
						#For each radial slice
						for i in range(0,R_mesh[l]-1):

							#Calculate cell area, all cells have the same area in a cartesian grid
							CellArea = (dr[l]/100.)*(dz[l]/100.)			#[m^2]
							CellVolume = CellArea*(dy[l]/100.)				#[m^3]

							#Calculate Power by multiplying power density for ring[i] by volume of ring[i]
							Power += PowerDensity[j][i]*CellVolume 			#[W]
						#endfor
					#endfor
					DepositedPowerList.append(Power)

					#Display power to terminal if requested.
					if print_totalpower == True:
						print(Dirlist[l])
						print(RequestedPowers[k]+' Deposited:',round(Power,4),'W')
					#endif
				#endif
			#endfor

			#==========#==========#

			#Plot and beautify each requested power deposition seperately.
			fig,ax = figure(image_aspectratio,1)
			Powers.append( DepositedPowerList[k*numfolders:(k+1)*numfolders] )
			TrendPlotter(ax,Powers[k],Xaxis,NormFactor=0)

			#Apply image options and axis labels.
			Title = 'Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Power Deposited [W]'
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

			plt.savefig(DirTrends+RequestedPowers[k]+' Deposition Trends'+ext)
			clearfigures(fig)
		#endfor

		#Write data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII, TotalPowerASCII = create_new_folder(DirTrends, 'Trend_Data'), [Xaxis]
			for k in range(0,len(RequestedPowers)): TotalPowerASCII.append(Powers[k])
			write_data_to_file(TotalPowerASCII, DirASCII + 'RFPower_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
		#endif

		#Plot a comparison of all power depositions requested.
		fig,ax = figure(image_aspectratio,1)
		for k in range(0,len(RequestedPowers)):TrendPlotter(ax,Powers[k],Xaxis,NormFactor=0)

		#Apply image options and axis labels.
		Title = 'Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Power Deposited [W]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

		plt.savefig(DirTrends+'Power Deposition Comparison'+ext)
		clearfigures(fig)
	#endif



	#====================================================================#
						#ION/NEUTRAL THRUST ANALYSIS#
	#====================================================================#

	#ABORT DIAGNOSTIC UNLESS ARGON IS SUPPLIED, WILL FIX LATER!!!
	if 'AR3S' in list(set(FluidSpecies).intersection(Variables)):
		if savefig_trendphaseaveraged == True or print_thrust == True:

			#Create Trend folder to keep output plots.
			TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
			TrendVariable = ''.join(TrendVariable)												#Single string of chars
			DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

			#Initiate lists required for storing data.
			NeutralThrustlist,IonThrustlist,PresThrustlist,Thrustlist = list(),list(),list(),list()
			NeutralIsplist,IonIsplist,ThrustIsplist = list(),list(),list()
			Xaxis = list()

			#Extract Positive, Negative and Neutral species names (Excluding electrons)
			NeutralSpecies = list(set(FluidSpecies).intersection(Variables))
			PositiveIons = PosSpecies
			NegativeIons = NegSpecies[:-1]

			#For all folders.
			for l in range(0,numfolders):

				#Update X-axis with folder information.
				Xaxis.append(folder_name_trimmer(Dirlist[l]))

				#Extract data required for Thrust calculations, discharge plane (Z) = thrustloc.
				VariableIndices,VariableStrings = enumerate_variables(['VZ-NEUTRAL'], Header_TEC2D[l])
				NeutralVelocity = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc)
				VariableIndices,VariableStrings = enumerate_variables(['VZ-ION+'], Header_TEC2D[l])
				IonVelocity = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc)
				VariableIndices,VariableStrings = enumerate_variables(['FZ-AR3S'], Header_TEC2D[l])
				NeutralAxialFlux = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc)
				VariableIndices,VariableStrings = enumerate_variables(['FZ-AR+'], Header_TEC2D[l])
				IonAxialFlux = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc)
				VariableIndices,VariableStrings = enumerate_variables(['TG-AVE'], Header_TEC2D[l])
				NeutGasTemp = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc)
				VariableIndices,VariableStrings = enumerate_variables(['PRESSURE'], Header_TEC2D[l])
				try:
					Pressure = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc)
					PressureDown = ExtractRadialProfile(Data[l],VariableIndices[0],VariableStrings[0],thrustloc-1)
				except:
					Pressure = np.zeros(R_mesh[l]*2)
					PressureDown = np.zeros(R_mesh[l]*2)
				#endtry

				#Following calculations expect Torr, convert back to Torr from SI if required
				if Units == 'SI':
					for i in range(0,len(Pressure)): Pressure[i] = Pressure[i]/133.33				#[Torr]
					for i in range(0,len(PressureDown)): PressureDown[i] = PressureDown[i]/133.33	#[Torr]
				#endif

				#Define which gas is used and calculate neutral mass per atom.
				NeutralIsp,IonIsp = list(),list()
				Argon,Xenon = 39.948,131.29			 #[amu]
				NeutralMass = Argon*1.67E-27		 #[Kg]

				#Choose which method to solve for thrust: 'ThermalVelocity','AxialMomentum'
				if GlobThrustMethod == 'ThermalVelocity':
					#Technique assumes cylindrical geometry, cartesian geometry will be overestimated.
					#Integrates neutral momentum loss rate based on neutral gas temperature.
					#Assumes angularly symmetric temperature and Maxwellian velocity distribution.
					NeutralThrust = 0
					for i in range(0,R_mesh[l]):

						# Determine Aperture Cell X-sectional Area, based on mesh geometry
						if IXZlist[l] == 0:
							# Calculate cylindrical radial plane area of a ring at radius [i],
							# correcting for central r=0.
							Circumference = 2*np.pi*(i*(dr[l]/100))		#m
							CellArea = Circumference*(dr[l]/100)		#m^2
							if CellArea == 0:
								CellArea = np.pi*((dr[l]/100)**2)		#m^2
							#endif

						elif IXZlist[l] == 1:
							# Calculate cartesian rectilinear planar slice of radius R, and depth dy,
							# slices are orientated "into" the screen as viewed in 2D
							CellArea = (dr[l]/100)*(dy[l]/100)			#m^2
						#endif

						#Calculate most probable neutral velocity based on temperature
						MeanVelocity = np.sqrt( (2*1.38E-23*NeutGasTemp[i])/(NeutralMass) )  	#m/s

						#If current cell is gas phase (Pressure > 0.0), calculate thrust
						if Pressure[i] > 0.0:
							#Calculate Neutral mass flow rate and integrate thrust via F = (dm/dt)Ve.
							NeutralMassFlowRate = NeutralAxialFlux[i]*NeutralMass*CellArea	#Kg/s
							NeutralExitVelocity = NeutralVelocity[i]						#m/s
							NeutralThrust += NeutralMassFlowRate * NeutralExitVelocity 		#N
							if NeutralExitVelocity > 0:
								NeutralIsp.append(NeutralExitVelocity)
							#endif
						#endif
					#endfor

					#Add neutral thrust and Isp to arrays (dummy variables not calculated)
					NeutralThrustlist.append( round(NeutralThrust*1000,5) )		#mN
					Thrustlist.append( round( NeutralThrust*1000,5) )			#mN
					NeutralIsp = (sum(NeutralIsp)/len(NeutralIsp))/9.81			#s
					Thrust,ThrustIsp = NeutralThrust,NeutralIsp					#N,s
					IonThrust,IonIsp = 1E-30,1E-30								#'Not Calculated'
					PresThrust = 1E-30											#'Not Calculated'
				#endif

				#====================#

				elif GlobThrustMethod == 'AxialMomentum':
					#Technique assumes cylindrical geometry, cartesian geometry will be overestimated.
					#Integrates ion/neutral momentum loss rate and differental pressure for concentric rings.
					#Assumes pressure differential, ion/neutral flux equal for all angles at given radii.

					#CellArea increases from central R=0.
					#Ensure pressure index aligns with radial index for correct cell area.
					#ONLY WORKS WHEN SYMMETRY OPTION IS ON, NEED A MORE ROBUST METHOD!
					Pressure,PressureDown = Pressure[0:R_mesh[l]][::-1],PressureDown[0:R_mesh[l]][::-1]
					NeutralVelocity = NeutralVelocity[0:R_mesh[l]][::-1]
					NeutralAxialFlux = NeutralAxialFlux[0:R_mesh[l]][::-1]
					IonVelocity = IonVelocity[0:R_mesh[l]][::-1]
					IonAxialFlux = IonAxialFlux[0:R_mesh[l]][::-1]
					#ONLY WORKS WHEN SYMMETRY OPTION IS ON, NEED A MORE ROBUST METHOD!

					PresThrust,NeutralThrust,IonThrust = 0,0,0
					for i in range(0,R_mesh[l]):

						# Determine Aperture Cell X-sectional Area, based on mesh geometry
						if IXZlist[l] == 0:
							# Calculate cylindrical radial plane area of a ring at radius [i],
							# correcting for central r=0.
							Circumference = 2*np.pi*(i*(dr[l]/100))		#m
							CellArea = Circumference*(dr[l]/100)		#m^2
							if CellArea == 0:
								CellArea = np.pi*((dr[l]/100)**2)		#m^2
							#endif

						elif IXZlist[l] == 1:
							# Calculate cartesian rectilinear planar slice of radius R, and depth dy,
							# slices are orientated "into" the screen as viewed in 2D
							CellArea = (dr[l]/100)*(dy[l]/100)			#m^2
						#endif

						#Calculate differential pressure between thrustloc-(thrustloc+1)
						if Pressure[i] > 0.0:
							DiffPressure = (Pressure[i]-PRESOUT[l])*133.33			#N/m^2
	#						DiffPressure = (Pressure[i]-PressureDown[i])*133.33		#N/m^2
							PresThrust += DiffPressure*CellArea						#N
						else:
							PresThrust += 0.0
						#endif

						#Calculate Neutral mass flow rate and integrate thrust via F = (dm/dt)Ve.
						NeutralMassFlowRate = NeutralAxialFlux[i]*NeutralMass*CellArea	#Kg/s
						NeutralExitVelocity = NeutralVelocity[i]						#m/s
						NeutralThrust += NeutralMassFlowRate * NeutralExitVelocity 		#N
						if NeutralExitVelocity > 0:
							NeutralIsp.append(NeutralExitVelocity)
						#endif

						#Calculate Ion mass flow rate and integrate thrust via F = (dm/dt)Ve.
						IonMassFlowRate = IonAxialFlux[i]*NeutralMass*CellArea	#Kg/s
						IonExitVelocity = IonVelocity[i]*1000					#m/s
						IonThrust += IonMassFlowRate * IonExitVelocity 			#N
						if IonExitVelocity > 0:
							IonIsp.append(IonExitVelocity)
						#endif
					#endfor
					if len(IonIsp) == 0: IonIsp.append(np.nan)
					if len(NeutralIsp) == 0: NeutralIsp.append(np.nan)

					#Add total thrust and calculate Isp of each component
					Thrust = PresThrust + NeutralThrust + IonThrust				#N
					NeutralFraction = NeutralThrust/(Thrust-PresThrust)			#Ignore dP/dz
					IonFraction = IonThrust/(Thrust-PresThrust)					#Ignore dP/dz

					IonIsp = (sum(IonIsp)/len(IonIsp))/9.81						#s
					NeutralIsp = (sum(NeutralIsp)/len(NeutralIsp))/9.81			#s
					ThrustIsp = NeutralFraction*NeutralIsp+IonFraction*IonIsp 	#s

					NeutralThrustlist.append( round(NeutralThrust*1000,5) )		#mN
					IonThrustlist.append( round(IonThrust*1000,5) )				#mN
					PresThrustlist.append( round(PresThrust*1000,5) )			#mN
					Thrustlist.append( round(Thrust*1000,5) )					#mN
					NeutralIsplist.append( round(NeutralIsp,5) )				#s
					IonIsplist.append( round(IonIsp,5) )						#s
					ThrustIsplist.append( round(ThrustIsp,5) )					#s
				#endif

				#====================#

				#Display thrust to terminal if requested.
				if print_thrust == True:
					print(Dirlist[l], '@ Z=',round(thrustloc*dz[l],2),'cm')
					print('NeutralThrust', round(NeutralThrust*1000,2), 'mN @ ', round(NeutralIsp,2),'s')
					print('IonThrust:', round(IonThrust*1000,4), 'mN @ ', round(IonIsp,2),'s')
					print('D-Pressure:', round(PresThrust*1000,4), 'mN')
					print('Thrust:',round(Thrust*1000,4),'mN @ ', round(ThrustIsp,2),'s')
					print('')
				#endif
			#endfor

			#Write data to ASCII format datafile if requested.
			if write_ASCII == True:
				DirASCII = create_new_folder(DirTrends, 'Trend_Data')
				write_data_to_file(Xaxis + ['\n'], DirASCII + 'Thrust_Trends', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(Thrustlist + ['\n'], DirASCII + 'Thrust_Trends', 'a',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(ThrustIsplist, DirASCII + 'Thrust_Trends', 'a',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif

			#=====#=====#

			#Plot total thrust and ion/neutral components.
			fig,ax1 = figure(image_aspectratio,1)
			TrendPlotter(ax1,Thrustlist,Xaxis,Marker='ko-',NormFactor=0)
			TrendPlotter(ax1,NeutralThrustlist,Xaxis,Marker='r^-',NormFactor=0)
			TrendPlotter(ax1,PresThrustlist,Xaxis,Marker='bs-',NormFactor=0)
			TrendPlotter(ax1,IonThrustlist,Xaxis,Marker='mv-',NormFactor=0)

			#Apply image options and save figure.
			Title='Thrust at Z='+str(round(thrustloc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Thrust F$_{T}$ [mN]'
			ax1.legend(['Total Thrust','Neutral Thrust','Pressure Thrust','Ion Thrust'], fontsize=18, frameon=False)
			ImageOptions(fig,ax1,Xlabel,Ylabel,Title,Crop=False)

			plt.savefig(DirTrends+'Thrust Trends'+ext)
			clearfigures(fig)

			#=====#=====#

			#Plot Specific Impulse for total thrust and ion/neutral components.
			fig,ax1 = figure(image_aspectratio,1)
			TrendPlotter(ax1,ThrustIsplist,Xaxis,Marker='ko-',NormFactor=0)
			TrendPlotter(ax1,NeutralIsplist,Xaxis,Marker='r^-',NormFactor=0)
	#		TrendPlotter(ax1,IonIsplist,Xaxis,Marker='bs-',NormFactor=0)

			#Apply image options and save figure.
			Title = 'Specific Impulse at Z='+str(round(thrustloc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Specific Impulse I$_{sp}$ [s]'
			ax1.legend(['Total I$_{sp}$','Neutral Component','Ion Component'], fontsize=18, frameon=False)
			ImageOptions(fig,ax1,Xlabel,Ylabel,Title,Crop=False)

			plt.savefig(DirTrends+'Isp Trends'+ext)
			clearfigures(fig)
		#endif
	#endif



	#====================================================================#
						#PHASE-AVERAGED SHEATH TRENDS#
	#====================================================================#

	if savefig_trendphaseaveraged == True or print_sheath == True:
	#NB: 	This diagnostic is very out of date, particularily the sheathROI treatment...
	#		Could be easily worked into the savefig_movieicp1D diagnostic or simply re-written
	#		Might be worth considering an overhaul... but it works for now.

		#Create Trend folder to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

		#Initialize any required lists.
		Xaxis,SxLocExtent,SxMaxExtent = list(),list(),list()

		#Obtain sheathROI and sourcewidth automatically if none are supplied.
		if len(sheathROI) != 2:
			#image_radialcrop Convert to Cells
			#image_axialcrop Convert to Cells
			#Use axialcrop or radialcrop to set automatic ROI!
			Start,End = 34,72				#AUTOMATIC ROUTINE REQUIRED#
			sheathROI = [Start,End]			#AUTOMATIC ROUTINE REQUIRED#
		#endif
		if len(sourcewidth) == 0:
			#Take Variable that is zero in metals (Density?)
			#Take Axial/Radial slice depending on sheath direction.
			#Find Cell distance from zero to 'wall' at electrodeloc.
			#Convert to SI [cm], set to automatic width.
			sourcewidth = [0.21]			#AUTOMATIC ROUTINE REQUIRED#
		#endif

		SxMeanExtent,SxMeanExtentArray = list(),list()
		#For all selected simulations, obtain Xaxis, sheath value and save to array.
		for l in range(0,numfolders):
			Xaxis.append(folder_name_trimmer(Dirlist[l]))

			#Obtain sheath thickness array for current folder
			Sx,SxAxis = CalcSheathExtent(folderidx=l)

			#Calculate mean sheath extent across ROI. On failure provide null point for sheath thickness.
			try:
				SxMeanExtentArray = list()
				for i in range(sheathROI[0],sheathROI[1]):	SxMeanExtentArray.append(Sx[i])
				SxMeanExtent.append(sum(SxMeanExtentArray)/len(SxMeanExtentArray))
			except:
				SxMeanExtent.append( np.nan )
			#endtry

			#Extract maximum sheath thickness from within region of interest
			try: SxMaxExtent.append( ((sourcewidth[0]*dr[l])-max(Sx[sheathROI[0]:sheathROI[1]]))*10 )
			except: SxMaxExtent.append( np.nan )

			#Extract sheath width adjacent to powered electrode
			#loc = electrodeloc[0]		#Radial
			loc = electrodeloc[1] 		#Axial
			try: SxLocExtent.append( ((sourcewidth[0]*dr[l])-Sx[loc])*10 )
			except:	SxLocExtent.append( np.nan )
		#endfor

		#===============================#

		#Write trend data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = create_new_folder(DirTrends, 'Trend_Data')
			write_data_to_file(Xaxis, DirASCII + 'Sx-Avg_Trends', 'w',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			write_data_to_file('\n', DirASCII + 'Sx-Avg_Trends', 'w',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			write_data_to_file(SxMaxExtent, DirASCII + 'Sx-Avg_Trends', 'a',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
		#endif

		#Generate figure and plot trends.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(ax,SxMeanExtent,Xaxis,NormFactor=0)

		#Apply image options and axis labels.
		Title = 'Maximum Sheath Extension With Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Sheath Extension [mm]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

		plt.savefig(DirTrends+'Sheath Extension (Phase-Averaged)'+ext)
		clearfigures(fig)
	#endif


	#====================================================================#
							#KNUDSEN NUMBER ANALYSIS#
	#====================================================================#

	#Only perform on bulk fluid dynamics relevent species.
	if bool(set(FluidSpecies).intersection(Variables)) == True:
		if savefig_trendphaseaveraged == True or print_Knudsennumber == True:

			#Create Trend folder to keep output plots.
			TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
			TrendVariable = ''.join(TrendVariable)												#Single string of chars
			DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

			#Initiate lists required for storing data.
			KnudsenAverage,Xaxis = list(),list()

			#For all folders - calculate Knudsen number for each cell of TECPLOT2D data
			for l in range(0,numfolders):

				#Using effective radius of argon in this calculation.
				Dimentionality = 2*(Radius[l]/100)		#[m]
				CrossSection = np.pi*((7.1E-11)**2)		#[m^2]

				#Extract data for the neutral flux and neutral velocity.
				VariableIndices,VariableStrings = enumerate_variables(FluidSpecies, Header_TEC2D[l])
				Sx,SxAxis = CalcSheathExtent(folderidx=l)

				#Update X-axis with folder information.
				Xaxis.append(folder_name_trimmer(Dirlist[l]))

				#Create empty image array based on mesh size and symmetry options.
				numrows = len(Data[l][0])/R_mesh[l]
				Image = np.zeros([Z_mesh[l],R_mesh[l]])

				#Produce Knudsen number 2D image using density image.
				for j in range(0,Z_mesh[l]):
					for i in range(0,R_mesh[l]):
						Start = R_mesh[l]*j
						Row = Z_mesh[l]-1-j

						#Convert units to SI
						LocalDensity = (Data[l][VariableIndices[0]][Start+i])*1E6					#[m-3]

						try:
							KnudsenNumber = (1/(LocalDensity*CrossSection*Dimentionality))		#[-]
						except:
							KnudsenNumber = 0													#[-]
						#endtry
						Image[Row,i] = KnudsenNumber
					#endfor
				#endfor

				#Display average Knudsen number to terminal if requested.
				KnudsenAverage.append( sum(Image)/(len(Image[0])*len(Image)) )
				if print_Knudsennumber == True:
					print(Dirlist[l])
					print('Average Knudsen Number:', KnudsenAverage[l])
				#endif

				#Create new folder to keep 2D output plots.
				Dir2Dplots = create_new_folder(Dirlist[l], 'TECPlot2D')
				#Write image data to ASCII format datafile if requested.
				if write_ASCII == True:
					DirASCII = create_new_folder(Dir2Dplots, '2DPlots_Data')
					write_data_to_file(Image, DirASCII + 'Kn', 'w',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif

				#Label and save the 2D Plots.
				extent,aspectratio = DataExtent(l)
				fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio)
				if image_plotsheath in ['Radial','Axial']:
					PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l],Orientation=image_plotsheath)
				#endif

				#Image plotting details
				Title = 'Knudsen Number Image for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				cax = Colourbar(ax,'Knudsen Number $K_{n}$',5,Lim=CbarMinMax(ax,Image))
				ImageOptions(fig,ax,Xlabel,Ylabel,Title)

				#Save Figure
				plt.savefig(Dir2Dplots+'2DPlot_Kn'+ext)
				clearfigures(fig)
			#endfor


			#Write trend data to ASCII format datafile if requested.
			if write_ASCII == True:
				DirASCII = create_new_folder(DirTrends, 'Trend_Data')
				write_data_to_file(Xaxis, DirASCII + 'Kn_Trends', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file('\n', DirASCII + 'Kn_Trends', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(KnudsenAverage, DirASCII + 'Kn_Trends', 'a',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif

			#Plot a comparison of all average Knudsen numbers.
			fig,ax = figure(image_aspectratio,1)
			TrendPlotter(ax,KnudsenAverage,Xaxis,NormFactor=0)

			#Image plotting details.
			Title = 'Average Knudsen Number with Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Average Knudsen Number $K_{n}$'
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

			#Save figure.
			plt.savefig(DirTrends+'KnudsenNumber_Comparison'+ext)
			clearfigures(fig)
		#endif
	#endif



	#====================================================================#
					#REYNOLDS NUMBER / SOUND SPEED ANALYSIS#
	#====================================================================#

	#Only perform on bulk fluid dynamics relevent species.
	if bool(set(FluidSpecies).intersection(Variables)) == True:
		if savefig_trendphaseaveraged == True or print_Reynolds == True:

			#Create Trend folder to keep output plots.
			TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
			TrendVariable = ''.join(TrendVariable)												#Single string of chars
			DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')

			#Initiate lists required for storing data.
			AverageSoundSpeed,Xaxis = list(),list()
			NeutralDensities = list()

			#For all folders.
			for l in range(0,numfolders):

				#Extract spatially resolved pressure and neutral densities.
				VariableIndices,VariableStrings = enumerate_variables(['PRESSURE'], Header_TEC2D[l])
				Pressure = ImageExtractor2D(Data[l][VariableIndices[0]],VariableStrings[0])
				Sx,SxAxis = CalcSheathExtent(folderidx=l)

				#If only single neutral species - extract that density
				VariableIndices,VariableStrings = enumerate_variables(FluidSpecies, Header_TEC2D[l])
				if len(VariableIndices) == 1:
					NeutralDensity = ImageExtractor2D(Data[l][VariableIndices[0]],VariableStrings[0])
				#If there are multiple neutral species, combine them to get total neutral density
				elif len(VariableIndices) > 1:
					for i in range(0,len(VariableIndices)):
						NeutralDensities.append( ImageExtractor2D(Data[l][VariableIndices[i]],VariableStrings[i]) )
					#endfor

					#Create empty neutral density array based on mesh size and symmetry options.
					numrows = len(Data[l][0])/R_mesh[l]
					NeutralDensity = np.zeros([Z_mesh[l],R_mesh[l]])

					#Combine all neutral densities to get total neutral density - if required.
					for i in range(0,len(NeutralDensities)):
						for j in range(0,len(NeutralDensities[i])):
							for k in range(0,len(NeutralDensities[i][j])):
								NeutralDensity[j][k] += NeutralDensities[i][j][k]
							#endfor
						#endfor
					#endfor
				#endif

				#Update X-axis with folder information.
				Xaxis.append(folder_name_trimmer(Dirlist[l]))

				#Calculate 2D sound speed image using neutral density and pressure
				Image = CalcSoundSpeed(NeutralDensity,Pressure,Dimension='2D')

				#Display mesh-averaged sound speed to terminal if requested.
				AverageSoundSpeed.append( sum(Image)/(len(Image[0])*len(Image)) )
				if print_Reynolds == True:
					print(Dirlist[l])
					print('Average Sound Speed:', AverageSoundSpeed[l])
				#endif

				#Create new folder to keep 2D output plots.
				Dir2Dplots = create_new_folder(Dirlist[l], 'TECPlot2D')
				#Write image data to ASCII format datafile if requested.
				if write_ASCII == True:
					DirASCII = create_new_folder(Dir2Dplots, '2DPlots_Data')
					write_data_to_file(Image, DirASCII + 'Cs', 'w',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif

				#Label and save the 2D Plots.
				extent,aspectratio = DataExtent(l)
				fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio)
				if image_plotsheath in ['Radial','Axial']:
					PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l],Orientation=image_plotsheath)
				#endif

				#Image plotting details
				#ERROR WITH IMAGE LIMIT - LIKELY DUE TO NANS - #Lim=CbarMinMax(ax,Image)
				Title = 'Sound Speed Image for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				cax = Colourbar(ax,'Sound Speed $C_{s}$ [m/s]',5,Lim=[])
				ImageOptions(fig,ax,Xlabel,Ylabel,Title)

				#Save Figure
				plt.savefig(Dir2Dplots+'2DPlot_Cs'+ext)
				clearfigures(fig)
			#endfor


			#Write trend data to ASCII format datafile if requested.
			if write_ASCII == True:
				DirASCII = create_new_folder(DirTrends, 'Trend_Data')
				write_data_to_file(Xaxis, DirASCII + 'Cs_Trends', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file('\n', DirASCII + 'Cs_Trends', 'w',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				write_data_to_file(AverageSoundSpeed, DirASCII + 'Cs_Trends', 'a',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif

			#Plot a comparison of all average Knudsen numbers.
			fig,ax = figure(image_aspectratio,1)
			TrendPlotter(ax,AverageSoundSpeed,Xaxis,NormFactor=0)

			#Image plotting details.
			Title = 'Average Sound Speed with Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Average Sound Speed'
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

			#Save figure.
			plt.savefig(DirTrends+'SoundSpeed_Comparison'+ext)
			clearfigures(fig)
		#endif
	#endif

	#===============================#
	#===============================#

	if any([savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_Reynolds, print_totalpower, print_DCbias, print_thrust, print_sheath]) == True:
		print('---------------------------')
		print('# Trend Processing Complete')
		print('---------------------------')
	#endif

	#=====================================================================#
	#=====================================================================#

















































	#====================================================================#
					#PHASE RESOLVED DIAGNOSTICS (REQ MOVIE1)#
	#====================================================================#

	#====================================================================#
							#1D PHASE RESOLVED MOVIES#
	#====================================================================#

	#Plot Phase-Resolved profiles with electrode voltage and requested variables.
	if savefig_phaseresolve1D == True:

		#Tnitiate any required lists.
		VoltageWaveforms,WaveformBiases,VariedValuelist = list(),list(),list()

		#for all folders.
		for l in range(0,numfolders):

			#Create global folders to keep output plots and collect graph title.
			Dirphaseresolved = create_new_folder(Dirlist[l], '1DPhase')
			VariedValuelist.append(folder_name_trimmer(Dirlist[l]))

			#Create VariableIndices for each folder as required. (Always get 'E','AR+','PPOT')
			PhaseData,Phaselist,proclist,VariableStrings = read_TEC2D_phase(l, PhaseVariables, Dir,
																			movie1, Header_movie1,
																			MinSharedVariables, Globalnumvars,
																			R_mesh, Z_mesh,
																			Units, AtomicSpecies)
			SxData,SxPhase,Sxproc,Sxvar = read_TEC2D_phase(l, ['E', 'AR+'], Dir,
														   movie1, Header_movie1,
														   MinSharedVariables, Globalnumvars,
														   R_mesh, Z_mesh,
														   Units, AtomicSpecies)
			PPOT = read_TEC2D_phase(l, ['PPOT'], Dir,
									movie1, Header_movie1,
									MinSharedVariables, Globalnumvars,
									R_mesh, Z_mesh,
									Units, AtomicSpecies)[2][0]

			#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
			Phaseaxis = GenerateAxis('Phase',ISYMlist[l],Phaselist)
			Raxis = GenerateAxis('Radial',ISYMlist[l])
			Zaxis = GenerateAxis('Axial',ISYMlist[l])

			#=============#

			#Extract waveforms from desired electrode locations.
			for j in range(0,len(waveformlocs)):
				VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
				WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
			#endfor
			ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(PhaseData,PPOT)

			#Plot the phase-resolved waveform.
			fig,ax = figure(image_aspectratio,1)

			ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
			for j in range(0,len(waveformlocs)):
				ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
				#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
			#endfor
			#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

			Title = 'Phase-Resolved Voltage Waveforms for ' + folder_name_trimmer(Dirlist[l])
			Legend = ['Waveform Vpp: '+str(round(ElectrodeVpp[2],2))+'V']
			Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
			ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

			plt.savefig(Dirphaseresolved+VariedValuelist[l]+' Waveform'+ext)
			clearfigures(fig)

			#Write waveform data in ASCII format if required.
			if write_ASCII == True:
				ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
				for j in range(0,len(waveformlocs)):
					ASCIIWaveforms.append(VoltageWaveforms[j])
				#endfor
				DirASCIIPhase = create_new_folder(Dirphaseresolved, '1DPhase_Data')
				write_data_to_file(ASCIIWaveforms, DirASCIIPhase + 'VoltageWaveforms',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif

			#==============#

			#for all requested variables.
			for i in tqdm(range(0,len(proclist))):

				#Create new folder to keep specific plots.
				DirMovieplots = create_new_folder(Dirphaseresolved, VariableStrings[i] + '_1Dphaseresolved')

				#Refresh lineout lists between variables.
				Lineouts,ProfileOrientation = list(),list()

				#Concatinate all requested lineouts together, keeping seperate orientation.
				for m in range(0,len(radialprofiles)):
					Lineouts.append(radialprofiles[m])
					ProfileOrientation.append('Radial')
				#endfor
				for m in range(0,len(axialprofiles)):
					Lineouts.append(axialprofiles[m])
					ProfileOrientation.append('Axial')
				#endfor

				#For all requested lineouts and orientations.
				for k in range(0,len(Lineouts)):

					#Refresh required lists.
					VariableMax,VariableMin = list(),list()

					#Create folders to keep output plots for each variable.
					if ProfileOrientation[k] == 'Axial':
						NameString = VariableStrings[i]+'_'+str(round(Lineouts[k]*dr[l],2))+'cm[R]'
					if ProfileOrientation[k] == 'Radial':
						NameString = VariableStrings[i]+'_'+str(round(Lineouts[k]*dz[l],2))+'cm[Z]'
					if savefig_phaseresolve1D == True:
						Dir1DProfiles = create_new_folder(DirMovieplots, NameString)
					#endif

					#Collect Normalization data for plotting.
					for j in range(0,len(Phaselist)):
						#Record local maximum and minimum for each phase.
						if ProfileOrientation[k] == 'Axial':
							Profile = ExtractAxialProfile(PhaseData[j],proclist[i],VariableStrings[i],Lineouts[k])
						elif ProfileOrientation[k] == 'Radial':
							Profile = ExtractRadialProfile(PhaseData[j],proclist[i],VariableStrings[i],Lineouts[k])
						#endif
						Profile,Minimum,Maximum = Normalise(Profile)
						VariableMax.append(Maximum)
						VariableMin.append(Minimum)
					#endfor
					#Find global maximum and minimum for all phases.
					VariableMax = max(VariableMax)
					VariableMin = min(VariableMin)

					#for all recorded phases, plot spatially varying variable and waveform.
					for j in range(0,len(Phaselist)):
						Phase = int( round(Phaseaxis[j]*360.0,3) )		#[Deg]

						if ProfileOrientation[k] == 'Axial':
							ZlineoutLoc,axis = Lineouts[k],Zaxis
							# SJD NEED TO REVERSE 1D DATA AT READIN AND REMOVE [::-1] FROM HERE				!!! SJD
							# NOTE, ALSO NEED TO REMOVE FROM DATA SAVE ROUTINES BELOW (csv and ACII) 		!!! SJD
							phaseresolvedProfile = ExtractAxialProfile(PhaseData[j],proclist[i],VariableStrings[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],ISYMlist[l])[::-1]
							ProfileString = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm \n'
							Xlabel = 'Axial Distance Z [cm]'
						elif ProfileOrientation[k] == 'Radial':
							RlineoutLoc,axis = Lineouts[k],Raxis
							phaseresolvedProfile = ExtractRadialProfile(PhaseData[j],proclist[i],VariableStrings[i],RlineoutLoc,R_mesh[l],ISYMlist[l])
							ProfileString = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm \n'
							Xlabel = 'Radial Distance R [cm]'
						#endif

						#Create figures and plot the 1D profiles. (ax[0]=variable, ax[1]=waveform)
						if image_plotphasewaveform == True:
							fig,ax = figure(image_aspectratio,2)
							ax0 = ax[0]							# Image Sub-Fig (top fig)
							ax1 = ax[1]							# Waveform Sub-Fig (bottom fig)
						else:
							fix,ax = figure(image_aspectratio,1)
							ax0 = ax							# Image Sub-Fig (top fig)
						#endif
						Ylabel = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
						fig.suptitle('Phase-Resolved '+VariableStrings[i]+' for '+VariedValuelist[l]+ProfileString+str(Phaselist[j]), y=0.97, fontsize=16)

						#Plot profile and apply image options.
						ax0.plot(axis, phaseresolvedProfile, lw=2)
						ImageOptions(fig,ax0,Xlabel,Ylabel[i],Crop=False)
						ax0.set_ylim(VariableMin,VariableMax*1.02)

						if image_plotphasewaveform == True:
							#Plot waveform and apply image options.
							ax1.plot(Phaseaxis, ElectrodeWaveform, lw=2)
							ax1.axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
							Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
							ImageOptions(fig,ax1,Xlabel,Ylabel,Crop=False)

							#Clean up image subfigures
							fig.tight_layout()
							plt.subplots_adjust(top=0.90)
						#endif

						#NOTE:	zfill assumes phase < 999 degrees	(i.e. < 1e5)
						plt.savefig(Dir1DProfiles+NameString+'_'+str(Phase).zfill(4)+ext)
						clearfigures(fig)



						# Write data underpinning current figure in .csv format
						if Write_CSV == True:
							CSVDir = create_new_folder(Dir1DProfiles, '1DPhase_Data')
							CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
							CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
							CSVFilename = VariableStrings[i]+'.csv'
							CSVTitle = str(Dirlist[l])
							CSVLabel = str(Ylabel)
							CSVISYM = 'ISYM='+str(ISYMlist[l])
							CSVRotate = 'Rotate='+str(image_rotate)
							CSVMaxCYCL = "IMOVIE_FRAMES="+str( Phaselist[-1].strip("CYCL= ") )
							CSVCurCYCL = str( Phaselist[j].replace(" ", "") )

							# Write to .csv file, including full header on first CYCLE
							if j == 0:
								CSVHeader = [CSVTitle,CSVLabel,CSVMaxCYCL,CSVISYM,CSVRotate,CSVRMesh,CSVZMesh]
								write_to_csv(phaseresolvedProfile[::-1], CSVDir, CSVFilename, CSVHeader, mode='w')

							# Write to .csv file, including only CYCL number for all remaining CYCLEs
							elif j > 0:
								CSVHeader = [CSVCurCYCL]
								write_to_csv(phaseresolvedProfile[::-1], CSVDir, CSVFilename, CSVHeader, mode='a')
							#endif
						#endif

						#Write Phase data in ASCII format if required.
						if write_ASCII == True:
							DirASCIIPhase = create_new_folder(Dirphaseresolved, '1DPhase_Data')
							DirASCIIPhaseloc = create_new_folder(DirASCIIPhase, ProfileString[3:-2])
							Cycle = str( Phaselist[j].replace(" ", "") )
							SaveString = DirASCIIPhaseloc+VariableStrings[i]+'_'+Cycle
							write_data_to_file(phaseresolvedProfile[::-1], SaveString,
											   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
						#endif
					#endfor

					#Create .mp4 movie from completed images.
					Prefix = folder_name_trimmer(Dirlist[l]) + '_' + NameString
					if ffmpegMovies:
						make_movie(Dir1DProfiles, Prefix)
				#endfor
			#endfor
		#endfor

		print('---------------------------------------')
		print('# 1D Phase-Resolved Processing Complete')
		print('---------------------------------------')
	#endif




	#====================================================================#
						#2D PHASE RESOLVED MOVIES#
	#====================================================================#

	#Plot 2D images over all saved phase cycles with included wavevform guide.
	if savefig_phaseresolve2D == True:

		#Initialize required lists.
		VoltageWaveforms,WaveformBiases,VariedValuelist = list(),list(),list()
		PPOT = list()

		#for all folders being processed.
		for l in range(0,numfolders):

			#Create global folder to keep output plots and collect graph title.
			Dirphaseresolved = create_new_folder(Dirlist[l], 'Movie1/')
			VariedValuelist.append(folder_name_trimmer(Dirlist[l]))

			#Create VariableIndices for each folder as required
			PhaseData,Phaselist,proclist,VariableStrings = read_TEC2D_phase(l, PhaseVariables, Dir,
																			movie1, Header_movie1,
																			MinSharedVariables, Globalnumvars,
																			R_mesh, Z_mesh,
																			Units, AtomicSpecies)

			#Preload 'E','AR+' for sheath plotting
			if image_plotsheath in ['Radial','Axial']:
				SxData,SxPhase,Sxproc,Sxvar = read_TEC2D_phase(l, ['E', 'AR+'], Dir,
															   movie1, Header_movie1,
															   MinSharedVariables, Globalnumvars,
															   R_mesh, Z_mesh,
															   Units, AtomicSpecies)
			#endif

			#Preload 'PPOT' for voltage waveform plotting
			try:
				PPOT = read_TEC2D_phase(l, ['PPOT'], Dir,
										movie1, Header_movie1,
										MinSharedVariables, Globalnumvars,
										R_mesh, Z_mesh,
										Units, AtomicSpecies)[2][0]
				PPOTexists = True
			except:
				print('Warning: PPOT not in movie1.pdt, skipping waveform plotting')
				PPOTexists = False
			#endtry

			#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
			Phaseaxis = GenerateAxis('Phase',ISYMlist[l],Phaselist)
			Raxis = GenerateAxis('Radial',ISYMlist[l])
			Zaxis = GenerateAxis('Axial',ISYMlist[l])

			# Define image extent for directory 'l'
			extent, aspectratio = DataExtent(l)

			#=====##=====# STANDALONE WAVEFORM IMAGE #=====##=====#

			if PPOTexists == True:
				#Extract waveforms from desired electrode locations.
				for j in range(0,len(waveformlocs)):
					VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
					WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
				#endfor
				ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(PhaseData,PPOT)

				#Plot the phase-resolved waveform.
				fig,ax = figure(image_aspectratio,1)

				ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
				for j in range(0,len(waveformlocs)):
					ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
					#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
				#endfor
				#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

				Title = 'Phase-Resolved Voltage Waveforms for ' + folder_name_trimmer(Dirlist[l])
				Legend = ['Waveform Vpp: '+str(round(ElectrodeVpp[2],2))+'V']
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
				ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

				plt.savefig(Dirphaseresolved+VariedValuelist[l]+' Waveform'+ext)
				clearfigures(fig)

				#Write PROES data in ASCII format if required.
				if write_ASCII == True:
					ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
					for j in range(0,len(waveformlocs)):
						ASCIIWaveforms.append(VoltageWaveforms[j])
					#endfor
					DirASCIIPhase = create_new_folder(Dirphaseresolved, '2DPhase_Data')
					write_data_to_file(ASCIIWaveforms, DirASCIIPhase + 'VoltageWaveforms',
									   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
				#endif
			#endif
			fig,ax = figure(image_aspectratio,1)

			#=====##=====# PHASE-RESOLVED IMAGE1.PDT VARIABLES #=====##=====#

			#for all variables requested by the user.
			for i in tqdm(range(0,len(proclist))):

				#Create new folder to keep specific plots.
				DirMovieplots = create_new_folder(Dirphaseresolved, VariableStrings[i] + '_2Dphaseresolved/')

				#Obtain maximum and minimum values of current variable over all phases.
				MinLim,MaxLim = list(),list()
				for j in range(0,len(Phaselist)):
					Image = ImageExtractor2D(PhaseData[j][proclist[i]],VariableStrings[i])
					MinLim.append( CbarMinMax(ax,Image,Symmetry=False)[0] )
					MaxLim.append( CbarMinMax(ax,Image,Symmetry=False)[1] )
				#endfor
				CbarLimits = [min(MinLim),max(MaxLim)]

				#Reshape specific part of 1D Data array into 2D image for plotting.
				for j in range(0,len(Phaselist)):

					#Convert phase 'CYCL' to degrees
					Phase = int( round(Phaseaxis[j]*360.0,3) )		#[Deg]

					#Extract full 2D image for further processing.
					Image = ImageExtractor2D(PhaseData[j][proclist[i]],VariableStrings[i])
					#Extract Ni and Ne variables for sheath processing.
					if image_plotsheath in ['Radial','Axial']:
						Ne = SxData[j][Sxproc[Sxvar.index('E')]]
						Ni = SxData[j][Sxproc[Sxvar.index('AR+')]]
						Sx,SxAxis = CalcSheathExtent(folderidx=l,Phase=j,Ne=Ne,Ni=Ni)
					#endif

					#Create figure and axes, plot image on top and waveform underneath.
					if image_plotphasewaveform == True and PPOTexists == True:
						fig,ax = figure(aspectratio,2)
						ax0 = ax[0]							# Image Sub-Fig (top fig)
						ax1 = ax[1]							# Waveform Sub-Fig (bottom fig)
					else:
						fix,ax = figure(aspectratio,1)
						ax0 = ax							# Image Sub-Fig (top fig)
					#endif

					#Plot 2D image of variable[i] at phase[j]
					fig,ax0,im,Image = ImagePlotter2D(Image,extent,aspectratio,VariableStrings[i],fig,ax0)

					#=====##=====# IMAGE OVERLAYS #=====##=====#

					# Overlay vector streamplot if exists
					Radial,Axial = enumerate_vectors(VariableStrings[i], Header_movie1[l])

					# Confirm Variable[k] has both vector counterparts (e.g. FR-AR3S, FZ-AR3S)
					VectorVariablesExist = True
					if None in Radial or None in Axial: VectorVariablesExist = False

					# Overlay variable[k]'s vector components onto 2D flood plot
					if image_plotvector == True and VectorVariablesExist == True:

						# Create meshgrid
						R_Space = np.linspace(0,Radius[l],R_mesh[l])			# Radial mesh locations [mm]
						Z_Space = np.linspace(0,Height[l],Z_mesh[l])			# Axial mesh locations [mm]
						if image_rotate == True:
							R_Space = np.linspace(0,Height[l],Z_mesh[l])		# Radial location is now axial location
							Z_Space = np.linspace(0,Radius[l],R_mesh[l])		# Axial location is now radial location
						#endif
						R, Z = np.meshgrid(R_Space, Z_Space)
							
						# See definition of "Radial" and "Axial" above
						UR = ImageExtractor2D(PhaseData[j][Radial[1]],Radial[0])  		# Radial vector magnitude
						UZ = ImageExtractor2D(PhaseData[j][Axial[1]] ,Axial[0] )  		# Axial vector magnitude
						if image_rotate == True:
							UR = ImageExtractor2D(PhaseData[j][Axial[1]] ,Axial[0])  	# Radial magnitude DATA is now Axial
							UZ = ImageExtractor2D(PhaseData[j][Radial[1]],Radial[0] )  	# Axial magnitude DATA is now Radial
							if image_plotsymmetry == False:
								UR = UR.transpose()[::-1]								# Without symmetry, Axis Zero is on LHS
								UZ = UZ.transpose()[::-1]
							elif image_plotsymmetry == True:
								UR = UR.transpose()										# With symmetry, Axis Zero is on RHS
								UZ = UZ.transpose()
							#endif
						#endif
						VectorLength = np.sqrt(UR**2 + UZ**2)

						# Streamplot fails if provided "zeros", inform user and skip to next variable
						try:
							Density = image_vectordensity
							Linewidth = image_vectorlw
							ax.streamplot(R, Z, UR, UZ, color=VectorLength, cmap='viridis', density=Density, linewidth=Linewidth)
						except:
							print('Warning: Invalid streamplot for variable: '+VariableStrings[i])
						#endtry
					#endif

					#=====#

					#Plot sheath onto ax1 if requested
					if image_plotsheath in ['Radial','Axial']:
						PlotSheathExtent(SxAxis,Sx,ax0,ISYMlist[l],Orientation=image_plotsheath)
					#endif

					#=====#

					# Plot waveform onto ax1 if requested
					if image_plotphasewaveform == True and PPOTexists == True:

						#Plot phase-resolved 'PPOT' waveform at 'electrodeloc'
						ax1.plot(Phaseaxis, ElectrodeWaveform, lw=2)
						ax1.axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
						ax1Xlabel,ax1Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'

						#Apply image options
						ImageOptions(fig,ax1,ax1Xlabel,ax1Ylabel,Crop=False)
						ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

						#Add Invisible Colourbar to sync X-axis
						InvisibleColourbar(ax0)

						#Clean up image subfigures
						fig.tight_layout()
						plt.subplots_adjust(top=0.90)
					#endif

					#=====##=====# Image Beautification #=====##=====#

					# Define title, labels, etc...
					Title = 'Phase-Resolved '+VariableStrings[i]+'\n'+str(Phaselist[j])
					if image_rotate == True:	Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
					elif image_rotate == False:	Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
					#endif

					# Add Colourbar (must happen before image cropping!)
					Cbarlabel = variable_label_maker(VariableStrings, Units, image_logplot, AtomicSpecies)
					cax = Colourbar(ax0,Cbarlabel[i],image_cbarbins,Lim=CbarLimits)

					# Crop image dimensions to [image_radialcrop,image_axialcrop]
					# Also resets cbar min/max to cropped region min/max
					if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
						CropImage(ax0,Rotate=image_rotate)
					#endif

					# Apply mesh geometry to image
					ImageGeometry(fig,ax0,image_plotsymmetry)

					# Apply image options, and enforce overides if requested
					ImageOptions(fig,ax0,Xlabel,Ylabel,Title,Crop=False)


					#=====##=====# Image I/O #=====##=====#

					#NOTE:	zfill assumes phase < 9999 degrees	(i.e. < 1e5)
					savefig(DirMovieplots+VariableStrings[i]+'_'+str(Phase).zfill(4)+ext)
					clearfigures(fig)

					#=====#

					# Write data underpinning current figure in .csv format
					if Write_CSV == True:
						CSVDir = create_new_folder(DirMovieplots, '2DPhase_Data')
						CSVRMesh = 'R_Mesh [Cells] '+str(R_mesh[l])+'  :: dR [cm/cell] '+str(dr[l])
						CSVZMesh = 'Z_Mesh [Cells] '+str(Z_mesh[l])+'  :: dZ [cm/cell] '+str(dz[l])
						CSVFilename = VariableStrings[i]+'.csv'
						CSVTitle = str(Dirlist[l])
						CSVLabel = str(Cbarlabel)
						CSVISYM = 'ISYM='+str(ISYMlist[l])
						CSVRotate = 'Rotate='+str(image_rotate)
						CSVMaxCYCL = "IMOVIE_FRAMES="+str( Phaselist[-1].strip("CYCL= ") )
						CSVCurCYCL = str( Phaselist[j].replace(" ", "") )

						# Write to .csv file, including full header on first CYCLE
						if j == 0:
							CSVHeader = [CSVTitle,CSVLabel,CSVMaxCYCL,CSVISYM,CSVRotate,CSVRMesh,CSVZMesh]
							write_to_csv(Image, CSVDir, CSVFilename, CSVHeader, mode='w')

						# Write to .csv file, including only CYCL number for all remaining CYCLEs
						elif j > 0:
							CSVHeader = [CSVCurCYCL]
							write_to_csv(Image, CSVDir, CSVFilename, CSVHeader, mode='a')
						#endif
					#endif

					#=====#

					#Write Phase data in ASCII format
					if write_ASCII == True:
						DirASCIIPhase = create_new_folder(DirMovieplots, '2DPhase_Data')
						DirASCIIPhaseVar = create_new_folder(DirASCIIPhase, VariableStrings[i])
						Cycle = str( Phaselist[j].replace(" ", "") )
						SaveString = DirASCIIPhaseVar+VariableStrings[i]+'_'+Cycle
						write_data_to_file(Image, SaveString,
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif
				#endfor

				#=====#

				#Create .mp4 movie from completed images.
				Prefix = folder_name_trimmer(Dirlist[l])
				if ffmpegMovies:
					make_movie(DirMovieplots, Prefix + '_' + VariableStrings[i])
			#endfor
		#endfor

		print('---------------------------------------')
		print('# 2D Phase-Resolved Processing Complete')
		print('---------------------------------------')
	#endif




	#====================================================================#
					#PHASE TRENDS & SHEATH DYNAMICS#
	#====================================================================#

	#Process phase resolved data from multiple folders to extract trends.
	if savefig_sheathdynamics == True:

		#Initiate arrays between folders
		VoltageWaveforms,WaveformBiases = list(),list()
		SxMaxExtTrend,SxMeanExtTrend= list(),list()
		SxMaxVelTrend,SxMeanVelTrend = list(),list()
		SxDynRangeTrend = list()
		VariedValuelist = list()

		#Read data from each simulation folder individually, saves on RAM.
		for l in tqdm(range(0,numfolders)):

			#Create global folders to keep output plots.
			TrendVariable = list(filter(lambda x: x.isalpha(), folder_name_trimmer(Dirlist[0])))	#List of discrete chars
			TrendVariable = ''.join(TrendVariable)												#Single string of chars
			Dirphaseresolved = create_new_folder(Dirlist[l], '2DPhase/')
			DirTrends = create_new_folder(os.getcwd() + '/', TrendVariable + ' Trends')
			DirSheath = create_new_folder(DirTrends, 'Sheath Trends')

			#Create VariableIndices for each folder as required. (Always get 'E','AR+','PPOT')
			SxData,SxPhase,Sxproc,Sxvar = read_TEC2D_phase(l, PhaseVariables + ['E', 'AR+', 'PPOT'], Dir,
														   movie1, Header_movie1,
														   MinSharedVariables, Globalnumvars,
														   R_mesh, Z_mesh,
														   Units, AtomicSpecies)
			VariedValuelist.append(folder_name_trimmer(Dirlist[l]))

			#Extract waveforms from desired electrode locations.
			PPOT = Sxproc[Sxvar.index('PPOT')]
			for j in range(0,len(waveformlocs)):
				VoltageWaveforms.append(WaveformExtractor(SxData,PPOT,waveformlocs[j])[0])
				WaveformBiases.append(WaveformExtractor(SxData,PPOT,waveformlocs[j])[1])
			#endfor
			ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(SxData,PPOT)

			### CURRENTLY ONLY AXIAL METHOD IS EMPLOYED ###
			#Axial sheath array (Sx) is calculated using radial integrations.
			#Radial sheath array (Sx) is calculated using axial integrations.
			Orientation = 'Axial'
			if Orientation == 'Axial': loc = electrodeloc[1]		#Axial depth where sheath plotted
			elif Orientation == 'Radial': loc = electrodeloc[0]		#Radial depth where sheath plotted
			Phaseaxis = GenerateAxis('Phase',Isym=ISYMlist[l])

			#=============#

			SxLoc = list()
			#For all phases, calculate sheath width and record sheath width at electrodeloc
			for k in range(0,len(SxPhase)):
				#Extract Ni and Ne variables for sheath processing.
				Ne = SxData[k][Sxproc[Sxvar.index('E')]]
				Ni = SxData[k][Sxproc[Sxvar.index('AR+')]]

				#calculate sheath width employing 'E' and 'AR+'
				Sx,SxAxis = CalcSheathExtent(folderidx=l,Orientation=image_plotsheath,Ne=Ne,Ni=Ni)
				for j in range(0,len(Sx)):
					Sx[j] = Sx[j]*10					#Convert to mm
				#endfor

				SxLoc.append(Sx[loc])
			#endfor

			#Determine phase-averaged sheath proprties, removing 'nans' unless all data is 'nan'
			SxLocNoNaN = [x for x in SxLoc if np.isnan(x) == False]
			if len(SxLocNoNaN) > 0:
				SxDynRangeTrend.append(max(SxLocNoNaN)-min(SxLocNoNaN))		#Dynamic Range
				SxMeanExtTrend.append(sum(SxLocNoNaN)/len(SxLocNoNaN))		#Mean Extension
				SxMaxExtTrend.append(max(SxLocNoNaN))						#Max Extension
			else:
				SxDynRangeTrend.append( np.nan )
				SxMeanExtTrend.append( np.nan )
				SxMaxExtTrend.append( np.nan )
			#endif

			#Calculate phase-averaged (mean) sheath velocity.
			#Assumes one sheath collapse and one sheath expansion per rf-cycle.
			RFPeriod = 1.0/FREQMIN[l]										#[s]
			MeanSheathExtent = (sum(SxLoc)/len(SxLoc))/1E6					#[km]
			MeanSheathVelocity = (2*MeanSheathExtent)/RFPeriod				#[km/s]
			SxMeanVelTrend.append( MeanSheathVelocity )						#[km/s]

			#Calculate maximum instantaneous sheath velocity.
			#Assumes sheath collapse velocity > sheath expansion velocity
			Collapsed,CollapsedPhase = min(SxLoc),SxLoc.index(min(SxLoc))
			Extended,ExtendedPhase = max(SxLoc),SxLoc.index(max(SxLoc))

			SheathExtension = (Extended-Collapsed)/1000.0  					#[m]
			PHASEEResolution = 1.0/(FREQMIN[l]*len(Phaseaxis))				#[s]
			SheathTime = (ExtendedPhase-CollapsedPhase)*PHASEEResolution		#[s]
			try:
				MaxSheathVelocity = SheathExtension/SheathTime				#[m/s]
				SxMaxVelTrend.append(MaxSheathVelocity/1000.0)				#[km/s]
			except:
				SxMaxVelTrend.append(0.0)									#[km/s]
			#endtry

			#=============#

			#Scale sheath extension by required number of phasecycles.
			ScaledSxLoc = list()
			for n in range(0,int(phasecycles*len(SxLoc))):
				Index = n % len(SxLoc)		#Modulo index for multiple phasecycles
				ScaledSxLoc.append(SxLoc[Index])
			#endfor
			SxLoc=ScaledSxLoc

			#Print results to terminal if requested.
			if print_sheath == True:
				print(Dirlist[l][2:-1])
				print('Sheath Extension:',round(SheathExtension*1E3,2),'[mm]')
				print('Sheath Collapse Time:',round(SheathTime*1E9,2),'[ns]')
				print('Mean Sheath Velocity:',round(MeanSheathVelocity/1E3,2),'[km/s]')
				print('Max Sheath Velocity:',round(MaxSheathVelocity/1E3,2),'[km/s]')
				print('')
			#endif

			#=============#

			#Plot phase-resolved sheath extension [ax0] and voltage waveform [ax1] for current folder
			fig,ax = figure(image_aspectratio,2,shareX=True)
			ax[0].plot(Phaseaxis,SxLoc, lw=2)
			Title = 'Phase-Resolved Sheath Extension for '+VariedValuelist[l]
			Ylabel = 'Sheath Extension [mm]'
			ImageOptions(fig,ax[0],Title=Title,Ylabel=Ylabel,Crop=False)

			ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
			ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
			Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
			ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)

			plt.savefig(Dirphaseresolved+VariedValuelist[l]+' SheathDynamics'+ext)
			clearfigures(fig)

			#Write phase-resolved sheath dynamics to ASCII format datafile if requested.
			if write_ASCII == True:
				DirASCII = create_new_folder(Dirphaseresolved, '2DPhase_Data')
				DirASCIISheath = create_new_folder(DirASCII, 'SheathDynamics')
				write_data_to_file(Phaseaxis + ['\n'] + SxLoc, DirASCIISheath + VariedValuelist[l] + 'SheathDynamics',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif
		#endfor

		#Write sheath trends to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = create_new_folder(DirTrends, 'Trend_Data')
			DirASCIISheath = create_new_folder(DirASCII, 'Sheath_Data')
			write_data_to_file(VariedValuelist + ['\n'] + SxMaxExtTrend, DirASCIISheath + 'MaxExtent_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			write_data_to_file(VariedValuelist + ['\n'] + SxMeanExtTrend, DirASCIISheath + 'MeanExtent_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			write_data_to_file(VariedValuelist + ['\n'] + SxDynRangeTrend, DirASCIISheath + 'DynamicRange_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			write_data_to_file(VariedValuelist + ['\n'] + SxMaxVelTrend, DirASCIISheath + 'MaxVelocity_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			write_data_to_file(VariedValuelist + ['\n'] + SxMeanVelTrend, DirASCIISheath + 'MeanVelocity_Trends',
							   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
		#endif


		#Plot maximum sheath extension trend for all folders
		fig,ax = figure(image_aspectratio)
		TrendPlotter(ax,SxMaxExtTrend,VariedValuelist,NormFactor=0)
		Title='Maximum Sheath Extension W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Max Sheath Extension [mm]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirSheath+'Max Sheath Extension Trends'+ext)
		clearfigures(fig)

		#Plot mean sheath extension trend for all folders
		fig,ax = figure(image_aspectratio)
		TrendPlotter(ax,SxMeanExtTrend,VariedValuelist,NormFactor=0)
		Title='Mean Sheath Extension W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Mean Sheath Extension [mm]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirSheath+'Mean Sheath Extension Trends'+ext)
		clearfigures(fig)

		#Plot sheath dynamic range (max-min extension) trend for all folders
		fig,ax = figure(image_aspectratio)
		TrendPlotter(ax,SxDynRangeTrend,VariedValuelist,NormFactor=0)
		Title='Sheath Dynamic Range W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Sheath Dynamic Range [mm]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirSheath+'Sheath Dynamic Range Trends'+ext)
		clearfigures(fig)

		#Plot maximum sheath velocity trend for all folders
		fig,ax = figure(image_aspectratio)
		TrendPlotter(ax,SxMaxVelTrend,VariedValuelist,NormFactor=0)
		Title='Maximum Sheath Velocity W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Sheath Velocity [km/s]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirSheath+'Max Sheath Velocity Trends'+ext)
		clearfigures(fig)

		#Plot mean sheath velocity trend for all folders
		fig,ax = figure(image_aspectratio)
		TrendPlotter(ax,SxMeanVelTrend,VariedValuelist,NormFactor=0)
		Title='Phase-averaged Sheath Velocity W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Sheath Velocity [km/s]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirSheath+'Mean Sheath Velocity Trends'+ext)
		clearfigures(fig)


		print('----------------------------------------')
		print('# Phase-Resolved Trend Analysis Complete')
		print('----------------------------------------')
	#endif



	#====================================================================#
						#SIMULATED PROES DIAGNOSTIC#
	#====================================================================#

	#Plot Phase-Resolved profiles with electrode voltage and requested variables.
	if savefig_PROES == True:
		VariedValuelist = list()

		#for all folders.
		for l in range(0,numfolders):

			#Create global folders to keep output plots and collect graph title.
			VoltageWaveforms,WaveformBiases = list(),list()
			DirPROES = create_new_folder(Dirlist[l], 'PROES')

			#Update X-axis with folder information.
			VariedValuelist.append(folder_name_trimmer(Dirlist[l]))

			#Create VariableIndices for each folder as required. (Always get 'E','AR+','PPOT')
			PhaseData,Phaselist,proclist,varlist = read_TEC2D_phase(l, PhaseVariables, Dir,
																	movie1, Header_movie1,
																	MinSharedVariables, Globalnumvars,
																	R_mesh, Z_mesh,
																	Units, AtomicSpecies)
			try:
				PPOT = read_TEC2D_phase(l, ['PPOT'], Dir,
										movie1, Header_movie1,
										MinSharedVariables, Globalnumvars,
										R_mesh, Z_mesh,
										Units, AtomicSpecies)[2][0]
			except:
				PPOT = list()
				for i in range(0,180): PPOT.append(i)
			#endtry
			if image_plotsheath in ['Radial','Axial']:
				SxData,SxPhase,Sxproc,Sxvar = read_TEC2D_phase(l, ['E', 'AR+'], Dir,
															   movie1, Header_movie1,
															   MinSharedVariables, Globalnumvars,
															   R_mesh, Z_mesh,
															   Units, AtomicSpecies)
			#endif

			#Generate SI scale axes for lineout plots.
			Phaseaxis = GenerateAxis('Phase',ISYMlist[l],Phaselist)		#[omega*t/2pi]
			Raxis = GenerateAxis('Radial',ISYMlist[l])					#[cm]
			Zaxis = GenerateAxis('Axial',ISYMlist[l])					#[cm]

			#=============#

			#Extract waveforms from desired electrode locations.
			for j in range(0,len(waveformlocs)):
				VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
				WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
			#endfor
			try:
				ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(PhaseData,PPOT)
			except:
				ElectrodeWaveform = list()
				for i in range(0,180): ElectrodeWaveform.append(i)
			#endtry


			#Plot the phase-resolved waveform.
			fig,ax = figure(image_aspectratio,1)

			ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
			for j in range(0,len(waveformlocs)):
				ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
				#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
			#endfor
			#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

			Title = 'Phase-Resolved Voltage Waveforms for ' + folder_name_trimmer(Dirlist[l])
			Legend = ['Waveform Vpp: '+str(round(ElectrodeVpp[2],2))+'V']
			Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
			ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

			plt.savefig(DirPROES+VariedValuelist[l]+' Waveform'+ext)
			clearfigures(fig)

			#Write PROES data in ASCII format if required.
			if write_ASCII == True:
				ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
				for j in range(0,len(waveformlocs)):
					ASCIIWaveforms.append(VoltageWaveforms[j])
				#endfor
				DirASCIIPROES = create_new_folder(DirPROES, 'PROES_Data')
				write_data_to_file(ASCIIWaveforms, DirASCIIPROES + 'VoltageWaveforms',
								   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
			#endif

			#==============#


			#for all requested variables.
			for i in tqdm(range(0,len(proclist))):
				#Refresh lineout lists between variables.
				Lineouts,ProfileOrientation = list(),list()

				#Concatinate all requested lineouts together, keeping seperate orientation.
				for m in range(0,len(radialprofiles)):
					Lineouts.append(radialprofiles[m])
					ProfileOrientation.append('Radial')
				#endfor
				for m in range(0,len(axialprofiles)):
					Lineouts.append(axialprofiles[m])
					ProfileOrientation.append('Axial')
				#endfor

				#For all requested profiles and orientations.
				for k in range(0,len(Lineouts)):
					#Refresh required lists.
					VariableMax,VariableMin = list(),list()
					PhaseSx = list()
					PROES = list()

					#for all recorded phases, plot spatially varying variable and waveform.
					for j in range(0,len(Phaselist)):
						#Refresh lists between each phasecycle.
						IntegratedDoFArray,DoFArrays = list(),list()

						#Collect each profile for stitching into a PROES image if required.
						if DoFwidth > 0:
							#Determine range of lineouts within the depth of field.
							DOFRegion = [(Lineouts[k]-DoFwidth),(Lineouts[k]+DoFwidth)]
							#If DOF extends beyond mesh, alert user and abort diagnostic.
							if any(DOFRegion) < 0:
								print('----------------------------------')
								print('Depth-of-Field Extends Beyond Mesh')
								print('----------------------------------')
								break
							#endif

							#Collect profiles from DOF region and transpose to allow easy integration.
							for LineoutLoc in range(DOFRegion[0],DOFRegion[1]):
								if ProfileOrientation[k] == 'Radial': DoFArrays.append(ExtractRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc))
								elif ProfileOrientation[k] == 'Axial': DoFArrays.append(ExtractAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
								#endif
							#endfor
							DoFArrays = np.asarray(DoFArrays).transpose().tolist()

							#Integrate DoF profiles spatially, obtaining PROES profile for phase 'j'
							for n in range(0,len(DoFArrays)):
								IntegratedDoFArray.append( sum(DoFArrays[n])/(DoFwidth*2+1) )
							#endif
							PROES.append(IntegratedDoFArray)

						#If no DoF then simply collect lineout from required location.
						elif DoFwidth == 0:
							LineoutLoc = Lineouts[k]
							if ProfileOrientation[k] == 'Radial':
								PROES.append(ExtractRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc))
							elif ProfileOrientation[k] == 'Axial':
								PROES.append(ExtractAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
							#endif
						#endif

						if image_plotsheath in ['Radial','Axial']:
							#Extract Ni and Ne variables and perform sheath processing.
							Ne = SxData[j][Sxproc[Sxvar.index('E')]]
							Ni = SxData[j][Sxproc[Sxvar.index('AR+')]]
							#Extract the full spatial sheath extent for phase 'j'
							#and save sheath extent at PROES profile location
							Sx,SxAxis = CalcSheathExtent(folderidx=l,Orientation=ProfileOrientation[k],Phase=j,Ne=Ne,Ni=Ni)
							PhaseSx.append(Sx[Lineouts[k]])
						#endif
					#endfor

					#Scale PROES image and Sx arrays by required number of phasecycles.
					ScaledPROES,ScaledPhaseSx,ScaledPhaseSxAxis = list(),list(),list()
					for n in range(0,int(phasecycles*len(PROES))):
						Index = n % len(PROES)		#Modulo index for multiple phasecycles
						ScaledPROES.append(PROES[Index])
						if image_plotsheath in ['Radial','Axial']:
							ScaledPhaseSx.append(PhaseSx[Index])
						#endif
					#endfor
					try: PROES,PhaseSx = ScaledPROES,ScaledPhaseSx
					except: PROES = ScaledPROES

					#Create figure and rotate PROES such that phaseaxis aligns with waveform.
					fig,ax = figure(image_aspectratio,2,shareX=True)
					PROES = ndimage.rotate(PROES, 90)

					#Set phase axis and determine labels depending on image orientation
					x1,x2 = Phaseaxis[0],Phaseaxis[-1]
					if ProfileOrientation[k] == 'Axial':
						ProfileString = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm'
						NameString = varlist[i]+'_'+ProfileString[2::]
						Ylabel = 'Axial Distance Z [cm]'
						Crop = [image_axialcrop[::-1],image_radialcrop] #Reversed accounting for rotation.
						y1,y2 = Zaxis[-1],Zaxis[0]						#Reversed accounting for top origin.
					elif ProfileOrientation[k] == 'Radial' and int(ISYMlist[l]) == 1:
						ProfileString = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
						NameString = varlist[i]+ProfileString[2::]
						Ylabel = 'Radial Distance R [cm]'
						Crop = [image_radialcrop,image_axialcrop]
						y1,y2 = Raxis[-1],-Raxis[-1]
					elif ProfileOrientation[k] == 'Radial' and int(ISYMlist[l]) == 0:
						ProfileString = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
						NameString = varlist[i]+ProfileString[2::]
						Ylabel = 'Radial Distance R [cm]'
						Crop = [image_radialcrop,image_axialcrop]
						y1,y2 = Raxis[-1],0
					#endif

					#Create output folder for current profile
					DirPROESloc = create_new_folder(DirPROES, ProfileString[3::])

					#Create PROES image along line of sight with phase-locked waveform.
					fig.suptitle( 'Simulated '+varlist[i]+' PROES for '+VariedValuelist[l]+ProfileString+'\n DoF = '+str(round(((2*DoFwidth)+1)*dz[l],2))+' cm', y=0.95, fontsize=18)
					im = ax[0].contour(PROES,extent=[x1,x2,y1,y2],origin='lower')
					im = ax[0].imshow(PROES,extent=[x1,x2,y1,y2],origin='lower',aspect='auto')
					if image_plotsheath in ['Radial','Axial']:
						PlotSheathExtent(Phaseaxis,PhaseSx,ax[0],ISYMlist[l],Orientation=image_plotsheath)
					#endif

					#Apply image options and set axes size
					ImageOptions(fig,ax[0],Xlabel='',Ylabel=Ylabel,Crop=Crop)
					ax[0].set_xticks([])
					ax[0].set_xlim(x1,x2)
					ax[0].set_ylim(y1,y2)
					#Add Colourbar
					cbarlabel = variable_label_maker(varlist, Units, image_logplot, AtomicSpecies)
					Colourbar(ax[0],cbarlabel[i],image_cbarbins,Lim=CbarMinMax(ax,Image=PROES,PROES=ProfileOrientation[k]))

					#Plot Voltage Waveform.
					ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
					ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
					Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
					ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
					ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
					ax[1].set_xlim(x1,x2)
					#Add Invisible Colourbar to sync X-axis
					InvisibleColourbar(ax[1])

					#Cleanup layout and save images.
					fig.tight_layout()
					plt.subplots_adjust(top=0.85)
					plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' PROES'+ext)
					clearfigures(fig)

					#Write PROES data in ASCII format if required.
					if write_ASCII == True:
						DirASCIIPROES = create_new_folder(DirPROES, 'PROES_Data')
						DirASCIIPROESloc = create_new_folder(DirASCIIPROES, ProfileString[3::])
						write_data_to_file(PROES, DirASCIIPROESloc + varlist[i] + '_PROES',
										   R_mesh_i=R_mesh[l], dr_i=dr[l], Z_mesh_i=Z_mesh[l], dz_i=dz[l])
					#endif

					#==============##==============#
					#==============##==============#

					#Spatially Collapse 2D PROES image along the line of sight.
					PROES,TemporalPROES = PROES.transpose().tolist(),list()
					for m in range(0,len(PROES) ):
						TemporalPROES.append( (sum(PROES[m][::]))/(len(PROES)/2) )
					#endfor

					#Plot Spatially Collapsed PROES with phase axis.
					fig,ax = figure(image_aspectratio,2,shareX=True)
					ax[0].plot(Phaseaxis,TemporalPROES, lw=2)
					Title = 'Spatially Integrated '+varlist[i]+' for '+VariedValuelist[l]+ProfileString+'\n DoF = '+str(round(((2*DoFwidth)+1)*dz[l],2))+' cm'
					Ylabel = 'Spatially Integrated '+varlist[i]
					ImageOptions(fig,ax[0],Title=Title,Ylabel=Ylabel,Crop=False)
	#				ax[0].set_xlim(x1,x2)

					#Plot Waveform onto Temporally collapsed PROES.
					ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
					ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
					Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
					ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
	#				ax[1].set_xlim(x1,x2)

					plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' TemporalPROES'+ext)
					clearfigures(fig)

					#=====#=====#

					#					!!!NOTE!!!
					#RM: 	SpatialPROES array length mismatch with R,Z Axes
					#		try/except Hack no longer working as both axes are incorrect
					#		Not sure why data axes is different length, something wrong with sum?

					#Temporally collapse 2D PROES image through defined phase fraction.
	#				SpatialPROES = list()
	#				for m in range(0,len(PROES)):
	#					SpatialPROES.append( (sum(PROES[m][::])) ) 		#Needs 'maths'-ing
					#endfor

					#Plot Temporally Collapsed PROES with required axis.
	#				fig,ax = figure(image_aspectratio,2,shareX=True)
	#				try: ax[0].plot(Raxis,SpatialPROES, lw=2)		### HACKY ###
	#				except: ax[0].plot(Zaxis,SpatialPROES, lw=2)	### HACKY ###
	#				Xlabel = 'Phase [$\omega$t/2$\pi$]'
	#				Ylabel = 'Temporally Integrated '+varlist[i]
	#				ImageOptions(fig,ax[0],Xlabel,Ylabel,Crop=False)

					#Plot Waveform onto Spatially collapsed PROES.
	#				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
	#				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
	#				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
	#				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)

	#				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' SpatialPROES'+ext)
	#				clearfigures(fig)

					#==============##==============#

				if l == numfolders and k == len(Lineouts):
					print('-------------------------------')
					print('# PROES Image analysis Complete')
					print('-------------------------------')
				#endif
			#endfor
		#endfor
	#endif

	#===============================#

	if any([savefig_sheathdynamics, savefig_phaseresolve1D, savefig_phaseresolve2D, savefig_PROES]) == True:
		print('----------------------------------')
		print('# Phase Resolved Profiles Complete')
		print('----------------------------------')
	#endif

	#=====================================================================#
	#=====================================================================#



