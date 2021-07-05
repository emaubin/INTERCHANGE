#!/usr/bin/python3
# -*-coding:Utf-8 -*

import os
import sys
import subprocess
import time, datetime
import pandas as pd
import shutil
import importlib
from itertools import permutations
import argparse

__version__ = '1.0'

current_directory = os.getcwd()
final_directory = f'{current_directory}/INTERCHANGE_output'

##########################################################################################################################################
####################################################### Define the arguments #############################################################
##########################################################################################################################################

parser=argparse.ArgumentParser(description="Step 2 of INTERCHANGE pipeline : generate Index", add_help=False,usage="python3 index_param.py [options] -p")

required_args = parser.add_argument_group("Positional arguments")
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')

setting_args = parser.add_argument_group('Settings')
setting_args.add_argument('-k',dest='kmersize',type=str,default=30,help='Size of k-mers. Default [30]')

output_args = parser.add_argument_group("Output options")
output_args.add_argument('-o',dest='output_directory',type=str,default=final_directory, help='Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory')

other_args = parser.add_argument_group('Other')
other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help='Show this help message and exit')
other_args.add_argument('-v','--version', action='version',version='INTERCHANGE v' + __version__,help="Show program's version number and exit")


args=parser.parse_args()

output_dir=args.output_directory
paths=args.paths
kmersize=args.kmersize

if not os.path.exists(output_dir):
        os.makedirs(output_dir)

output_sp=f'{output_dir}/species'
if not os.path.exists(output_sp):
        os.makedirs(output_sp)

with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]
gt_path = param[3]

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

dict_index={}
with open(f"{output_sp}/fasta_tab.txt", 'r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split("\t")
        sp=line[0]
        fa=line[1]
        if sp not in dict_index:
            dict_index[sp]=fa

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

for key, value in dict_index.items():
    fa=value
    species=key
    suffixindex= f"{output_sp}/{species}/{species}"
    indexmer=f"{output_sp}/{species}/{species}"
    paramfile=f"{output_sp}/{species}/{species}_step2.param"
    with open(paramfile, 'w') as filout:
        filout.write(f'{fa}\n{suffixindex}\n{indexmer}\n{kmersize}\n{gt_path}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/INTERCHANGE-V.1.0/scripts/2.Index/index.py -p {paramfile} &')
    print(f"Job running : {key}")
