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

parser=argparse.ArgumentParser(description="Step 6 of INTERCHANGE pipeline : Search of homologous scaffolds", add_help=False,usage="python3 homologous_scfd_param.py [options] -p")

required_args = parser.add_argument_group("Positional arguments")
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')

setting_args = parser.add_argument_group('Settings')
setting_args.add_argument('-t',dest='threads',type=str,default=2,help='Number of CPU for alignments. Default [2]')

output_args = parser.add_argument_group("Output options")
output_args.add_argument('-o',dest='output_directory',type=str,default=final_directory, help='Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory')

other_args = parser.add_argument_group('Other')
other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help='Show this help message and exit')
other_args.add_argument('-v','--version', action='version',version='INTERCHANGE v' + __version__,help="Show program's version number and exit")

args=parser.parse_args()
output_dir=args.output_directory
paths=args.paths
threads=args.threads


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

# Annotation Directories
output_annot=f"{output_dir}/Annotation"
if not os.path.exists(output_annot):
        os.makedirs(output_annot)

#Homologs Directory
output_blast=f"{output_dir}/Homologous_scaffolds"
if not os.path.exists(output_blast):
        os.makedirs(output_blast)

with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]
blast_path = param[11]

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
            dico_compare[spindex,spread]=[spread,spindex]


##########################################################################################################################################
########################################## Write parameter in new script to run KmerAnalysis class #######################################
##########################################################################################################################################

for key, value in dico_compare.items():
    db=f"{output_assembly}/SPAdes_{key[1]}-reads_{key[0]}-index_cleaned/{key[1]}-reads_{key[0]}-index_scaffolds.fa_sup300.fa"
    query=f"{output_assembly}/SPAdes_{value[1]}-reads_{value[0]}-index_cleaned/{value[1]}-reads_{value[0]}-index_scaffolds.fa_sup300.fa"
    comp_blast=f"{output_blast}/{key[1]}-DB_{value[1]}-query"
    if not os.path.exists(comp_blast):
        os.makedirs(comp_blast)
    output=f"{comp_blast}/{key[1]}-DB_{value[1]}-query.bl6"
    paramfile=f'{comp_blast}/{key[1]}-DB_{value[1]}-query_step6.param'
    with open(paramfile, 'w') as filout:
        filout.write(f'{db}\n{query}\n{output}\n{threads}\n{output_sp}/all_species_comparisons.tab\n{output_blast}\n{blast_path}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/Pipeline_HT/scripts/6.Homologous_scaffolds/homologous_scfd.py -p {paramfile} &')
    print(f"Job running : index = {key[0]} ; reads = {key[1]}")

##########################################################################################################################################
