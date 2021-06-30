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


parser=argparse.ArgumentParser(description="Script to prepare and format input data for INTERCHANGE", add_help=False,usage="python3 format_param.py -i -p")

required_args = parser.add_argument_group("Positional arguments")
required_args.add_argument("-i",dest='Table',type=str,required=True, help='Input file containing Table of species.')
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')

setting_args = parser.add_argument_group('Settings')
setting_args.add_argument('-t',dest='thread',type=int,default=2, help='Number of CPU for gzip/gunzip. Default [2]')

output_args = parser.add_argument_group("Output options")
output_args.add_argument('-o',dest='output_directory',type=str,default=final_directory, help='Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory')

other_args = parser.add_argument_group('Other')
other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help='Show this help message and exit.')
other_args.add_argument('-v','--version', action='version',version='INTERCHANGE v' + __version__,help="Show program's version number and exit.")

args=parser.parse_args()

Table=args.Table
paths=args.paths
output_dir=args.output_directory
output_reads=f"{output_dir}/Reads"
output_sp=f'{output_dir}/species'
thread=args.thread

with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]

##########################################################################################################################################
############################## Ask if output directory already exists else create this directory #########################################
##########################################################################################################################################

if os.path.exists(output_dir):
    print('Ouput directory already exists')
else:
    os.makedirs(output_dir)

if not os.path.exists(output_sp):
    os.makedirs(output_sp)

if not os.path.exists(output_reads):
    os.makedirs(output_reads)

##########################################################################################################################################
###################################################### Import CSV file ###################################################################
##########################################################################################################################################

data=pd.read_csv(Table, sep='\t')
specie_name=data['tag']
for specie in specie_name:
    i=0
    for letter in specie:
        if letter.isspace()==True: #Test if there is space in different specie names
            i+=1
    if i!=0:
        specie_name=data['species'].str.replace(" ", "") # If there are space, script replaces spaces by no space

col1=data["fastq1"]
col2=data["fastq2"]


##########################################################################################################################################
################################################ Add items in diferrent lists ############################################################
##########################################################################################################################################

species_list=[]
fastq1_list=[]
fastq2_list=[]

for name in specie_name:
    if name not in species_list:
        species_list.append(name)

for fq1 in col1:
    if fq1 not in fastq1_list:
        fastq1_list.append(fq1)
for fq2 in col2:
    if fq2 not in fastq2_list:
        fastq2_list.append(fq2)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


i=0
fasta_dict={}
while i < len(species_list):
    fq1= f"{fastq1_list[i]}"
    fq2= f"{fastq2_list[i]}"
    outfq1=f"{output_reads}/{species_list[i]}_1.fq"
    outfq2=f"{output_reads}/{species_list[i]}_2.fq"
    specie_dir=f"{output_sp}/{species_list[i]}"
    if not os.path.exists(specie_dir):
        os.makedirs(specie_dir)
    fa=f"{specie_dir}/{species_list[i]}"
    if species_list[i] not in fasta_dict:
        fasta_dict[species_list[i]]=fa
    paramfile=f"{specie_dir}/{species_list[i]}_step1.param"
    with open(paramfile, 'w') as filout:
        filout.write(f'{fq1}\n{fq2}\n{outfq1}\n{outfq2}\n{fa}\n{thread}\n{prog_dir}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/INTERCHANGE-master/scripts/1.Genome_format/format.py -p {paramfile} &')
    print(f"Job running : {species_list[i]}")
    i+=1


fasta_tab=open(f"{output_sp}/fasta_tab.txt", "w")
for k, v in fasta_dict.items():
    fasta_tab.write(f"{k}\t{v}_reads.fa\n")

fasta_tab.close()
