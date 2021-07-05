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

parser=argparse.ArgumentParser(description="Step 3 of INTERCHANGE pipeline : Search identical k-mers", add_help=False,usage="python3 identical_kmers_param.py [options] -p ")

required_args = parser.add_argument_group("INTERCHANGE directory")
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')

setting_args = parser.add_argument_group('Settings')
setting_args.add_argument('-k',dest='kmersize',type=str,default=30,help='Size of k-mers. Default [30]')
setting_args.add_argument('-t',dest='threshold',type=str,default=50,help='Threshold for the k-mer filter. Default [50]')
setting_args.add_argument('-r',dest='readsize',type=str,default=150,help='Size of reads. Default [150]')

output_args = parser.add_argument_group("Output options")
output_args.add_argument('-o',dest='output_directory',type=str,default=final_directory, help='Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory')

other_args = parser.add_argument_group('Other')
other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help='Show this help message and exit')
other_args.add_argument('-v','--version', action='version',version='INTERCHANGE v' + __version__,help="Show program's version number and exit")

args=parser.parse_args()

output_dir=args.output_directory
kmersize=args.kmersize
readsize=args.readsize
threshold=args.threshold
paths=args.paths

#Variable for species directory with contain k-mer index for each species
output_sp=output_dir+'/species'
#Directory of Comparisons
output_comp=f'{output_dir}/Comparisons'
#Directory with files of coordinates of each kmers
coord_dir=f'{output_comp}/kmers_info'

# Creation of new output directory for comparisons
if not os.path.exists(output_comp):
    os.makedirs(output_comp)
if not os.path.exists(coord_dir):
    os.makedirs(coord_dir)

with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]
gt_path = param[3]



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

species_list=[]
comparisons_list=[]
with open(f"{output_sp}/fasta_tab.txt", 'r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split("\t")
        specie_name=line[0]
        reads=line[1]
        if specie_name not in species_list:
            species_list.append(specie_name)
        compare=permutations(species_list,2) # list of all pairwise comparisons
        for pairs in list(compare):
            if pairs not in comparisons_list:
                comparisons_list.append(pairs)


#Creat a file containing the list of all pairwise comparisons
comparisons_file=open(f'{output_sp}/all_species_comparisons.tab', 'w')
for pair in comparisons_list:
    pair1=pair[0]
    pair2=pair[1]
    if pair1!=pair2:
        comparisons_file.write(f'{pair[0]}\t{pair[1]}\n')
comparisons_file.close()

dico_compare={}
with open(f'{output_sp}/all_species_comparisons.tab','r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split('\t')
        spread=line[1]
        spindex=line[0]
        if spindex and spread not in dico_compare:
            dico_compare[spindex,spread]=[f'{output_sp}/{spindex}/{spindex}_{kmersize}mers_index',f'{output_sp}/{spread}/{spread}_reads.fa']


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################


for key, value in dico_compare.items():
    indexmer=value[0]
    query=value[1]
    output=f'{coord_dir}/{key[0]}-index_{key[1]}-reads_{kmersize}mers'
    kmercoord=f'{coord_dir}/{key[0]}-index_{key[1]}-reads.{kmersize}mers.info'
    outdir_KR=f'{output_comp}/{key[1]}-reads_{key[0]}-index'
    if not os.path.exists(outdir_KR):
        os.makedirs(outdir_KR)
    output_filter=f'{outdir_KR}/{key[1]}_reads_{key[0]}_index'
    kmerreads=f"{output_filter}.fa"
    paramfile=f'{outdir_KR}/{key[1]}-reads_{key[0]}-index.param'
    with open(paramfile, 'w') as filout:
        filout.write(f'{indexmer}\n{query}\n{output}\n{kmersize}\n{kmercoord}\n{threshold}\n{readsize}\n{output_filter}\n{kmerreads}\n{prog_dir}\n{gt_path}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/INTERCHANGE-V.1.0/scripts/3.Search_identical_kmers/identical_kmers.py -p {paramfile} &')
    print(f"Job running")
