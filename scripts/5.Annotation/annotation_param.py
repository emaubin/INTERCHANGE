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

parser=argparse.ArgumentParser(description="Step 5 of INTERCHANGE pipeline : Annotation of scaffolds", add_help=False,usage="python3 annotation_param.py [options] -p")

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

output_cdd=f"{output_annot}/cdd_delta"
if not os.path.exists(output_cdd):
        os.makedirs(output_cdd)

output_mcr=f"{output_annot}/mcr"
if not os.path.exists(output_mcr):
        os.makedirs(output_mcr)

output_Repbase=f"{output_annot}/repbase"
if not os.path.exists(output_Repbase):
        os.makedirs(output_Repbase)

output_otherDB=f"{output_annot}/otherDB"
if not os.path.exists(output_otherDB):
        os.makedirs(output_otherDB)


with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]
diamond_path = param[9]
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
        assembly_file=f"{output_assembly}/SPAdes_{spread}-reads_{spindex}-index_cleaned/scaffolds.fasta"
        if spindex and spread not in dico_compare:
            dico_compare[spindex,spread]=[assembly_file]

##########################################################################################################################################
########################################## Write parameter in new script to run KmerAnalysis class #######################################
##########################################################################################################################################

for key, value in dico_compare.items():
    species=key[1]+'-reads_'+key[0]+'-index'
    fasta=value[0]
    assembly=f"{output_assembly}/SPAdes_{key[1]}-reads_{key[0]}-index_cleaned/{species}_scaffolds.fa"
    paramfile=f'{output_assembly}/SPAdes_{key[1]}-reads_{key[0]}-index_cleaned/{key[1]}-reads_{key[0]}-index_step5.param'
    with open(paramfile, 'w') as filout:
        filout.write(f'{species}\n{fasta}\n{assembly}\n{output_Repbase}\n{output_mcr}\n{output_cdd}\n{output_otherDB}\n{threads}\n{prog_dir}\n{diamond_path}\n{blast_path}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/Pipeline_HT/scripts/5.Annotation/annotation.py -p {paramfile} &')
    print(f"Job running : index = {key[0]} ; reads = {key[1]}")

##########################################################################################################################################
