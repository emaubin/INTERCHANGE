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
final_directory = f'{current_directory}/INTERCHANGE_results'

##########################################################################################################################################
####################################################### Define the arguments #############################################################
##########################################################################################################################################

parser=argparse.ArgumentParser(description="Step 4 of INTERCHANGE pipeline : Cleaning and Assembly of reads", add_help=False,usage="python3 assembly_param.py [options] -P")

required_args = parser.add_argument_group("Positional arguments")
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')

setting_args = parser.add_argument_group('Settings')
setting_args.add_argument('-k',dest='kspades',type=str,default="auto",help='Size of k-mers for SPAdes assembly. Default [auto]')
setting_args.add_argument('-t',dest='threads',type=str,default=2,help='Number of CPU for SPAdes Assembly. Default [2]')

output_args = parser.add_argument_group("Output options")
output_args.add_argument('-o',dest='output_directory',type=str,default=final_directory, help='Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory')

other_args = parser.add_argument_group('Other')
other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help='Show this help message and exit')
other_args.add_argument('-v','--version', action='version',version='INTERCHANGE v' + __version__,help="Show program's version number and exit")

args=parser.parse_args()
output_dir=args.output_directory
paths=args.paths
threads=args.threads
kspades=args.kspades

args=parser.parse_args()

#Variable for species directory with contain k-mer index for each species
output_sp=output_dir+'/species'
if not os.path.exists(output_sp):
        os.makedirs(output_sp)

#Directory for Comparisons
output_comp=f'{output_dir}/Comparisons'
if not os.path.exists(output_comp):
        os.makedirs(output_comp)

#Directory for Assemblies
output_assembly=f"{output_dir}/Assembly"

if not os.path.exists(output_assembly):
        os.makedirs(output_assembly)

with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]
prinseq_path=param[5]
spades_path = param[7]


##########################################################################################################################################
################################################## Built Lists and Dictionaries ##########################################################
##########################################################################################################################################

dico_compare={}
with open(f'{output_sp}/all_species_comparisons.tab','r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split('\t')
        spread=line[1]
        spindex=line[0]
        if spindex and spread not in dico_compare:
            dico_compare[spindex,spread]=[f'{output_comp}/{spread}-reads_{spindex}-index/{spread}_reads_{spindex}_index',f'{output_sp}/{spread}/{spread}_reads.fa']

print(dico_compare)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

for key, value in dico_compare.items():
    kmerreads=f'{value[0]}.fa'
    good_reads=f'{value[0]}_cleaned'
    reads=value[1]
    output=f'{output_assembly}/SPAdes_{key[1]}-reads_{key[0]}-index_cleaned'
    fa1=f'{good_reads}_1.fa'
    fa2=f'{good_reads}_2.fa'
    paramfile=f'{output_comp}/{key[1]}-reads_{key[0]}-index_Step4.param'
    with open(paramfile, 'w') as filout:
        filout.write(f'{prog_dir}\n{kmerreads}\n{good_reads}\n{reads}\n{threads}\n{output}\n{fa1}\n{fa2}\n{kspades}\n{prinseq_path}\n{spades_path}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/Pipeline_HT/scripts/4.Assembly/assembly.py -p {paramfile} &')
    print(f"Job running : index = {key[0]} ; reads = {key[1]}")
