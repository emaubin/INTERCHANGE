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


parser=argparse.ArgumentParser(description="Step 8 of INTERCHANGE pipeline : Identification of Busco genes", add_help=False,usage="busco_identification_param.py -i -p")

required_args = parser.add_argument_group("Positional arguments")
required_args.add_argument("-i",dest='Table',type=str,required=True, help='Input file containing Table of species.')
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')

setting_args = parser.add_argument_group('Settings')
setting_args.add_argument('-m',dest='tmini',type=int,default=2, help='Number of CPU for Minimap2 mapping. Default [2]')
setting_args.add_argument('-s',dest='tspades',type=int,default=2, help='Number of CPU for SPAdes assembly. Default [2]')
setting_args.add_argument('-b',dest='tblast',type=int,default=2, help='Number of CPU for BLAST alignments. Default [2]')
setting_args.add_argument('-k',dest='kspades',type=str,default="auto",help='Size of k-mers for SPAdes assembly. Default [auto]')

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
busco_dir=f'{output_dir}/Busco_genes'
tmini=args.tmini
tspades=args.tspades
k=args.kspades
tblast=args.tblast

with open(paths, "r") as param:
    param= param.read().split("\n")

prog_dir = param[1]
spades_path = param[7]
minimap2_path=param[13]
samtools_path=param[15]
busco_path=param[17]
blast_path = param[11]


target=f"{busco_path}/BUSCO_viridiplantae_odb10_2021_homologs_400_plant_genomes.nt.fasta"

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

if not os.path.exists(busco_dir):
    os.makedirs(busco_dir)

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

print(species_list)
print(fastq1_list)
print(fastq2_list)

###########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


i=0
fasta_dict={}
while i < len(species_list):
    fq1= f"{fastq1_list[i]}"
    fq2= f"{fastq2_list[i]}"
    busco_sp_dir=f"{busco_dir}/{species_list[i]}"
    if not os.path.exists(busco_sp_dir):
        os.makedirs(busco_sp_dir)
    out_mapped=f"{busco_sp_dir}/{species_list[i]}_busco_mapped.bam"
    sp_reads_dir=f"{busco_sp_dir}/Reads"
    if not os.path.exists(sp_reads_dir):
        os.makedirs(sp_reads_dir)
    mapped_id=f"{sp_reads_dir}/{species_list[i]}_busco_mapped.id"
    reads_busco_id1=f"{sp_reads_dir}/{species_list[i]}_busco_reads_1.id"
    reads_busco_id2=f"{sp_reads_dir}/{species_list[i]}_busco_reads_2.id"
    reads_busco1=f"{sp_reads_dir}/{species_list[i]}_busco_reads_1.fa"
    reads_busco2=f"{sp_reads_dir}/{species_list[i]}_busco_reads_2.fa"
    reads_sp=f"{output_sp}/{species_list[i]}/{species_list[i]}_reads.fa"
    output_Assembly=f"{sp_reads_dir}/SPADes_{species_list[i]}_busco_genes"
    sp_scaffolds_busco=f"{output_Assembly}/{species_list[i]}_scaffolds_busco_genes.fasta"
    species=species_list[i]
    output_blast=f"{output_Assembly}/{species_list[i]}_scaffolds_busco_genes"
    paramfile=f"{sp_reads_dir}/{species_list[i]}_step8.param"
    with open(paramfile, 'w') as filout:
        filout.write(f'{minimap2_path}\n{tmini}\n{target}\n{fq1}\n{fq2}\n{out_mapped}\n{samtools_path}\n{mapped_id}\n{reads_busco_id1}\n{reads_busco_id2}\n{reads_sp}\n{reads_busco1}\n{reads_busco2}\n{spades_path}\n{tspades}\n{output_Assembly}\n{k}\n{sp_scaffolds_busco}\n{species}\n{prog_dir}\n{output_blast}\n{blast_path}\n{tblast}\n')
        filout.close()
    os.system(f'python3 {prog_dir}/INTERCHANGE-V.1.0/scripts/8.Busco_genes/busco_identification.py -p {paramfile} &')
    print(f"Job running : {species_list[i]}")
    i+=1


#########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
