# Transfers feature-barcode matrices for every library of the subproject from cluster to local

# Load packages
import os
import argparse
import subprocess

# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to transfer feature-barcodes matrices from cluster to lcal")
parser.add_argument("-l", "--libraries",
		    dest = "libraries",
		    action = "store",
		    nargs= "+",
		    default = None,
		    help = "List of libraries to transfer")
parser.add_argument("-p", "--subproject",
		    dest = "subproject",
		    action = "store",
		    default = None,
		    help = "Subproject id (e.j BCLLATLAS_01")
options = parser.parse_args()
libraries = options.libraries
subproject = options.subproject
project = subproject[:-3]

# Transfer files
for lib in libraries:
	lib_path_cluster = "/scratch/devel/rmassoni/{}/{}/jobs/{}/{}/outs".format(project, subproject, lib, lib)
	data_path_local = "data/{}/{}".format(subproject, lib)
	os.mkdir(data_path_local)
	subprocess.run(["scp", "-r", "rmassoni@172.16.10.20:{}/filtered_feature_bc_matrix".format(lib_path_cluster), data_path_local])
	subprocess.run(["scp", "rmassoni@172.16.10.20:{}/web_summary.html".format(lib_path_cluster), data_path_local])