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

parser=argparse.ArgumentParser(description="Step 9 of INTERCHANGE pipeline : High similarity", add_help=False,usage="python3 high_similarity_param [options] -p")

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


output_sp=f'{output_dir}/species'
busco_dir=f'{output_dir}/Busco_genes'
output_assembly=f'{output_dir}/Assembly'
output_annot=f'{output_dir}/Annotation'
output_HS_candidates=f'{output_dir}/HS_candidates'


if not os.path.exists(output_sp):
        os.makedirs(output_sp)

if not os.path.exists(busco_dir):
        os.makedirs(busco_dir)

if not os.path.exists(output_assembly):
        os.makedirs(output_assembly)

if not os.path.exists(output_HS_candidates):
        os.makedirs(output_HS_candidates)

with open(paths, "r") as param:
    param= param.read().split("\n")
prog_dir = param[1]
blast_path = param[11]

##########################################################################################################################################
################################################## Lists and Dictionaries ################################################################
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

#print(dico_compare)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

for key, value in dico_compare.items():
    busco_sp_dir_q=f"{busco_dir}/{key[0]}"
    sp_reads_dir_q=f"{busco_sp_dir_q}/Reads"
    output_Assembly_q=f"{sp_reads_dir_q}/SPADes_{key[0]}_busco_genes"
    busco_sp_dir_db=f"{busco_dir}/{value[0]}"
    sp_reads_dir_db=f"{busco_sp_dir_db}/Reads"
    output_Assembly_db=f"{sp_reads_dir_db}/SPADes_{value[0]}_busco_genes"
    db=f"{output_Assembly_db}/{value[0]}_scaffolds_busco_genes_vf.fa"
    query=f"{output_Assembly_q}/{key[0]}_scaffolds_busco_genes_vf.fa"
    output_blast=f"{busco_dir}/Busco_comparisons"
    if not os.path.exists(output_blast):
        os.makedirs(output_blast)
    comp_blast=f"{output_blast}/{value[0]}-DB_{key[0]}-query.bln6"
    paramfile=f'{output_blast}/{value[0]}-DB_{key[0]}-query_step9.param'
    species=key[1]+'-reads_'+key[0]+'-index'
    scfd_fasta=f"{output_assembly}/SPAdes_{key[1]}-reads_{key[0]}-index_cleaned/{species}_scaffolds.fa_sup300.fa"
    TE_annot=f"{output_annot}/Annotation_tab.txt_candidates_TE.txt"
    TE_validated=f"{output_HS_candidates}/TE_HSvalidation.txt"
    TE_ID=f"{output_HS_candidates}/TE_HSvalidation.id"
    TE_fasta=f"{output_HS_candidates}/TE_HSvalidation.fa"
    GENE_annot=f"{output_annot}/Annotation_tab.txt_candidates_GENE.txt"
    GENE_validated=f"{output_HS_candidates}/GENE_HSvalidation.txt"
    GENE_ID=f"{output_HS_candidates}/GENE_HSvalidation.id"
    GENE_fasta=f"{output_HS_candidates}/GENE_HSvalidation.fa"
    with open(paramfile, 'w') as filout:
        filout.write(f'{db}\n{query}\n{comp_blast}\n{threads}\n{output_sp}/all_species_comparisons.tab\n{output_blast}\n{blast_path}\n{species}\n{scfd_fasta}\n{TE_annot}\n{TE_validated}\n{TE_ID}\n{TE_fasta}\n{GENE_annot}\n{GENE_validated}\n{GENE_ID}\n{GENE_fasta}')
        filout.close()
    os.system(f'python3 {prog_dir}/INTERCHANGE-V.1.0/scripts/9_High_similarity/high_similarity.py -p {paramfile} &')
    print(f"Job running : index = {key[0]} ; reads = {value[0]}")

##########################################################################################################################################
