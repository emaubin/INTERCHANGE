#!/usr/bin/python3
# -*-coding:Utf-8 -*

import os
import sys
import subprocess
import time, datetime
from datetime import date
import pandas as pd
import shutil
import importlib
from itertools import permutations
import argparse

__version__ = '1.0'

current_directory = os.getcwd()
final_directory = f'{current_directory}/INTERCHANGE_results'

start=time.time()
today = date.today()
day = today.strftime("%b_%d_%Y")

##########################################################################################################################################
####################################################### Define the arguments #############################################################
##########################################################################################################################################

parser=argparse.ArgumentParser(description="Step 7 of INTERCHANGE pipeline : Generate an annotation table", add_help=False,usage="python3 annotation_table.py [options] -i -p")

required_args = parser.add_argument_group("Positional arguments")
required_args.add_argument("-i",dest='Table',type=str,required=True, help='Input file containing Table of species.')
required_args.add_argument('-p',dest='paths',type=str,required=True, help='File of tools paths.')


output_args = parser.add_argument_group("Output options")
output_args.add_argument('-o',dest='output_directory',type=str,default=final_directory, help='Output directory for INTERCHANGE results. Default: /INTERCHANGE_results in current directory')

other_args = parser.add_argument_group('Other')
other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help='Show this help message and exit')
other_args.add_argument('-v','--version', action='version',version='INTERCHANGE v' + __version__,help="Show program's version number and exit")


args=parser.parse_args()
Table=args.Table
paths=args.paths
output_dir=args.output_directory

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
            dico_compare[spindex,spread]=[f"{spread}-reads_{spindex}-index"]


for key, value in dico_compare.items():
    IDsup300=f"{output_assembly}/SPAdes_{key[1]}-reads_{key[0]}-index_cleaned/{value[0]}_scaffolds.fa.IDsup300"
    TABreciprocal=f"{output_blast}/Homologs_table.csv"
    BESTrepbase=f"{output_annot}/repbase/{value[0]}_sup300_Repbase.blx6.best"
    BESTcmr=f"{output_annot}/mcr/{value[0]}_scaffolds_DBMCR.bln6.best"
    BESTcdd=f"{output_annot}/cdd_delta/{value[0]}_scaffolds_cddDelta.blx6.best"
    BESTother_DB=f"{output_annot}/otherDB/{value[0]}_scaffolds_otherDB.bln6.best"
    TABannotation=f"{output_annot}/Annotation_{day}_tab.txt"

    dict_species={}
    dict_repbase={}
    dict_cdddelta={}
    dict_dbmcr={}
    dict_reciprocal1={}
    dict_reciprocal2={}
    dict_otherDB={}

    list_ID=[]
    with open(IDsup300,'r') as filin:
        for line in filin:
            line= line[:-1]
            id =line
            species=line.split("_")
            sp1=species[0]
            sp1=sp1.split("-")
            sp1=sp1[0]
            sp2=species[1]
            sp2=sp2.split("-")
            sp2=sp2[0]
            if id not in list_ID:
                list_ID.append(id)
            if sp1 not in dict_species:
                dict_species[sp1]=sp2
        filin.close()

    dict_sp1={}
    dict_sp2={}
    data=pd.read_csv(Table, sep='\t')
    genus=data['genus']
    species=data['species']
    type=data['type']
    tag=data['tag']
    genus_list=[]
    species_list=[]
    type_list=[]
    tag_list=[]
    for g in genus:
        genus_list.append(g)
    for s in species:
        species_list.append(s)
    for ty in type:
        type_list.append(ty)
    for ta in tag:
        tag_list.append(ta)

    i=0
    while i < len(tag_list):
        for k, v in dict_species.items():
            if tag[i]==k:
                if k not in dict_sp1:
                    dict_sp1[k]=[species_list[i],genus_list[i],type_list[i]]
            if tag[i]==v:
                if v not in dict_sp2:
                    dict_sp2[v]=[species_list[i],genus_list[i],type_list[i]]
        i+=1

    list_reciprocal1=[]
    list_reciprocal2=[]
    with open(TABreciprocal,'r') as filin:
        for line in filin:
            line=line[:-1]
            line=line.split('\t')
            ID1=line[0]
            ID2=line[1]
            pid1=line[2]
            len1=line[3]
            score1=line[4]
            scalen1=line[5]
            pid2=line[7]
            len2=line[8]
            score2=line[9]
            scalen2=line[10]
            if ID1 not in dict_reciprocal1:
                dict_reciprocal1[ID1]=[pid1,len1,score1,scalen1,scalen2,ID2]
            if ID2 not in dict_reciprocal2:
                dict_reciprocal2[ID2]=[pid2,len2,score2,scalen2,scalen1,ID1]
            if ID1 not in list_reciprocal1:
                list_reciprocal1.append(ID1)
            if ID2 not in list_reciprocal2:
                list_reciprocal2.append(ID2)
        filin.close()

    list_repbase=[]
    with open(BESTrepbase,'r') as filin:
        for line in filin:
            line=line[:-1]
            line=line.split("\t")
            ID=line[0]
            pid=float(line[2])
            length=int(line[3])
            score=float(line[11])
            annot=line[1]
            if ID not in dict_repbase:
                dict_repbase[ID]=[annot,pid,length,score]
            if ID not in list_repbase:
                list_repbase.append(ID)
        filin.close()

    list_cmr=[]
    with open(BESTcmr,'r') as filin:
        for line in filin:
            line=line[:-1]
            line=line.split("\t")
            ID=line[0]
            annot=line[1]
            pid=float(line[2])
            length=int(line[3])
            score=float(line[11])
            if ID not in dict_dbmcr:
                dict_dbmcr[ID]=[annot,pid,length, score]
            if ID not in list_cmr:
                list_cmr.append(ID)

    list_cdd=[]
    with open(BESTcdd,'r') as filin:
        for line in filin:
            line=line[:-1]
            line=line.split("\t")
            ID=line[0]
            pid=float(line[2])
            length=int(line[3])
            score=float(line[11])
            annot=line[1]
            if ID not in dict_cdddelta:
                dict_cdddelta[ID]=[annot,pid,length,score]
            if ID not in list_cdd:
                list_cdd.append(ID)
        filin.close()

    list_otherDB=[]
    with open(BESTother_DB,'r') as filin:
        for line in filin:
            line=line[:-1]
            line=line.split("\t")
            ID=line[0]
            pid=float(line[2])
            length=int(line[3])
            score=float(line[11])
            annot=line[1]
            if ID not in dict_otherDB:
                dict_otherDB[ID]=[annot,pid,length,score]
            if ID not in list_otherDB:
                list_otherDB.append(ID)
        filin.close()

    with open(TABannotation,'a') as filout:
        i =0
        for k, v in dict_species.items():
            sp1=k
            sp2=v
        while i < len(list_ID):

            if list_ID[i] in list_repbase:
                score1=dict_repbase.get(list_ID[i])[3]
                alilen1=dict_repbase.get(list_ID[i])[2]
                filout.write(f'{sp1}\t{dict_sp1.get(sp1)[0]}\t{dict_sp1.get(sp1)[1]}\t{dict_sp1.get(sp1)[2]}\t{sp2}\t{dict_sp2.get(sp2)[0]}\t{dict_sp2.get(sp2)[1]}\t{dict_sp2.get(sp2)[2]}\t{list_ID[i]}\t{dict_repbase.get(list_ID[i])[0]}\t{dict_repbase.get(list_ID[i])[1]}\t{dict_repbase.get(list_ID[i])[2]}\t{dict_repbase.get(list_ID[i])[3]}\t')
            else:
                score1=0
                alilen1=0
                filout.write(f'{sp1}\t{dict_sp1.get(sp1)[0]}\t{dict_sp1.get(sp1)[1]}\t{dict_sp1.get(sp1)[2]}\t{sp2}\t{dict_sp2.get(sp2)[0]}\t{dict_sp2.get(sp2)[1]}\t{dict_sp2.get(sp2)[2]}\t{list_ID[i]}\tN/A\tN/A\tN/A\tN/A\t')

            if list_ID[i] in list_cmr:
                score2=dict_dbmcr.get(list_ID[i])[3]
                alilen2=dict_dbmcr.get(list_ID[i])[2]
                filout.write(f'{dict_dbmcr.get(list_ID[i])[0]}\t{dict_dbmcr.get(list_ID[i])[1]}\t{dict_dbmcr.get(list_ID[i])[2]}\t{dict_dbmcr.get(list_ID[i])[3]}\t')
            else:
                score2=0
                alilen2=0
                filout.write(f'N/A\tN/A\tN/A\tN/A\t')

            if list_ID[i] in list_cdd:
                score3=dict_cdddelta.get(list_ID[i])[3]
                alilen3=dict_cdddelta.get(list_ID[i])[2]
                filout.write(f'{dict_cdddelta.get(list_ID[i])[0]}\t{dict_cdddelta.get(list_ID[i])[1]}\t{dict_cdddelta.get(list_ID[i])[2]}\t{dict_cdddelta.get(list_ID[i])[3]}\t')
            else:
                score3=0
                alilen3=0
                filout.write(f'N/A\tN/A\tN/A\tN/A\t')

            if list_ID[i] in list_otherDB:
                score4=dict_otherDB.get(list_ID[i])[3]
                alilen4=dict_otherDB.get(list_ID[i])[2]
                filout.write(f'{dict_otherDB.get(list_ID[i])[0]}\t{dict_otherDB.get(list_ID[i])[1]}\t{dict_otherDB.get(list_ID[i])[2]}\t{dict_otherDB.get(list_ID[i])[3]}\t')
            else:
                score4=0
                alilen4=0
                filout.write(f'N/A\tN/A\tN/A\tN/A\t')

            if list_ID[i] in list_reciprocal1:
                filout.write(f'{dict_reciprocal1.get(list_ID[i])[5]}\t{dict_reciprocal1.get(list_ID[i])[0]}\t{dict_reciprocal1.get(list_ID[i])[1]}\t{dict_reciprocal1.get(list_ID[i])[2]}\t{dict_reciprocal1.get(list_ID[i])[3]}\t{dict_reciprocal1.get(list_ID[i])[4]}\t')
            elif list_ID[i] in list_reciprocal2:
                filout.write(f'{dict_reciprocal2.get(list_ID[i])[5]}\t{dict_reciprocal2.get(list_ID[i])[0]}\t{dict_reciprocal2.get(list_ID[i])[1]}\t{dict_reciprocal2.get(list_ID[i])[2]}\t{dict_reciprocal2.get(list_ID[i])[3]}\t{dict_reciprocal2.get(list_ID[i])[4]}\t')
            else:
                filout.write(f'N/A\tN/A\tN/A\tN/A\tN/A\tN/A\t')

            if score1==0 and score2==0 and score3==0 and score4==0:
                filout.write(f'Unknown\n')
            elif score1 > score2 and score1 > score3 and score1 > score4:
                filout.write(f'TE\n')
            elif score2 > score1 and score2 > score3 and score2 > score4:
                filout.write(f'MCR\n')
            elif score3 > score1 and score3 > score2 and score3 > score4:
                filout.write(f'GENE\n')
            elif score4 > score1 and score4 > score2 and score4 > score3:
                filout.write(f'Other\n')
            elif score1==score2 and score3 < score1 or score2==score3 and score1 < score2 or score1==score3 and score2 < score1 or score1==score2==score3:
                if alilen1 > alilen2 and alilen1 > alilen3:
                    filout.write(f'TE\n')
                if alilen2 > alilen1 and alilen2 > alilen3:
                    filout.write(f'MCR\n')
                if alilen3 > alilen1 and alilen3 > alilen2:
                    filout.write(f'GENE\n')
            i+=1
        filout.close()
"""
candiate_dict={}
with open(TABannotation) as tab:
    for line in tab:
        line=line[:-1]
        line=line.split("\t")
        genus1=line[1]
        sp1=line[2]
        type1=line[3]
        scf1=line[8]
        genus2=line[5]
        sp2=line[6]
        type2=line[7]
        scf2=line[25]
        pid=line[26]
        alilen=line[27]
        tag_annot=line[31]
        if tag_annot == "TE":
            print(tag_annot)
            if pid != "N/A":
                print(pid)
                if pid >= "80":
                    if alilen >= "1000":
                        candiate_dict[scf1]=[genus1,sp1,type1,scf1,genus2,sp2,type2,scf2,pid,alilen,tag_annot,annot]
        elif tag_annot == "GENE":
            annot=line[17]
            if pid != "N/A":
                if pid >= "80":
                    if alilen >= "1000":
                        candiate_dict[scf1]=[genus1,sp1,type1,scf1,genus2,sp2,type2,scf2,pid,alilen,tag_annot,annot]

    tab.close()

candidate_tab=open(f"{output_annot}/Candidates_{day}_tab.txt", "w")
for k, v in candiate_dict.items():
    candidate_tab.write(f"{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\t{v[5]}\t{v[6]}\t{v[7]}\t{v[8]}\t{v[9]}\t{v[10]}\t{v[11]}\n")

candidate_tab.close()
"""
