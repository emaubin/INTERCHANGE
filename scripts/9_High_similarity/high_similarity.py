#!/usr/bin/python3
# -*-coding:Utf-8 -*


import os, sys
import subprocess
import time, datetime
import pandas as pd
from itertools import permutations
import argparse
import numpy as np

start=time.time()

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

class Bestblast():
    def best(self,bl1,bl_best1,bl2,bl_best2):
        self.bl1=bl1
        self.bl_best1=bl_best1
        self.bl2=bl2
        self.bl_best2=bl_best2

        os.system(f"sort -k12,12nr {bl1} | sort -u -k1,1 > {bl_best1}")
        os.system(f"sort -k12,12nr {bl2} | sort -u -k1,1 > {bl_best2}")
        return

class Blast():

    def makeblastdb(self,db,blast_path):
        self.db=db
        self.blast_path=blast_path

        return os.system(f"{blast_path}/makeblastdb -in {db} -dbtype nucl")


    def blastn(self,db,query,output,threads,blast_path):
        self.db=db
        self.query=query
        self.output=output
        self.threads=threads
        self.blast_path=blast_path

        return os.system(f"{blast_path}/blastn -task blastn -query {query} -db {db} -outfmt 6 -evalue 1e-20 -out {output} -num_threads {threads}")

class ComparisonBlast():

    def besthit(self,blast1,blast2,outputhit):
        self.blast1=blast1
        self.blast2=blast2
        self.outputhit=outputhit

        os.system('''awk -F "\t" '{split($1,a,"_");print $0"\t"a[5]"\t"$4/a[5]*100}' '''+blast1+''' > '''+blast1+'''.cov''')
        os.system('''awk -F "\t" '{split($1,a,"_");print $0"\t"a[5]"\t"$4/a[5]*100}' '''+blast2+''' > '''+blast2+'''.cov''')

        COMMAND1='''sort  -k12,12nr '''+blast1+'''.cov | awk '!seen[$1]++' | awk '{print $1"#"$2"\t"$3"\t"$4"\t"$12"\t"$13"\t"$14}' | sort -k1,1 > '''+blast1+'''.best'''
        COMMAND2='''sort  -k12,12nr '''+blast2+'''.cov | awk '!seen[$1]++' | awk '{print $2"#"$1"\t"$3"\t"$4"\t"$12"\t"$13"\t"$14}' | sort -k1,1 > '''+blast2+'''.best'''
        COMMAND3='''join -1 1 -2 1 '''+blast1+'''.best '''+blast2+'''.best | sed 's/ /\t/g' | sed 's/#/\t/g' > '''+outputhit
        subprocess.call(COMMAND1, shell=True)
        subprocess.call(COMMAND2, shell=True)
        subprocess.call(COMMAND3, shell=True)

        return

def getSeq(self,fasta_file,id_list,output_fasta):
    self.fasta_file=fasta_file
    self.id_list=id_list
    self.output_fasta=output_fasta
    return os.system(f'{prog_dir}/INTERCHANGE-V.1.0/Tools/getSeq.py -f {fasta_file} -l {id_list} -o {output_fasta}')

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

parser=argparse.ArgumentParser()
parser.add_argument('-p',dest='parameters',type=str,required=True, help='Define parameters file')
args=parser.parse_args()
paramfile=args.parameters


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

with open(paramfile, 'r') as param:
    param= param.read().split("\n")


db=param[0]
query=param[1]
output=param[2]
threads=param[3]
comparisons=param[4]
output_blast=param[5]
blast_path=param[6]
species=param[7]
scfd_fasta=param[8]
TE_annot=param[9]
TE_validated=param[10]
TE_ID=param[11]
TE_fasta=param[12]
GENE_annot=param[13]
GENE_validated=param[14]
GENE_ID=param[15]
GENE_fasta=param[16]

reciprocal=Blast()
reciprocal.makeblastdb(db,blast_path)
reciprocal.blastn(db,query,output,threads,blast_path)

###########################################################################################################################################
############################################## Calculation of HS ##########################################################################
###########################################################################################################################################


list_specie1=[]
with open(comparisons,'r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split('\t')
        spread=line[1]
        spindex=line[0]
        if spindex not in list_specie1:
            list_specie1.append(spindex)
    filin.close()

list_file=[]
i=0
while i<len(list_specie1):
    d=0+i
    while d<len(list_specie1):
        if list_specie1[i]!=list_specie1[d]:
            blast1=f"{output_blast}/{list_specie1[i]}-DB_{list_specie1[d]}-query.bln6"
            blast2=f"{output_blast}/{list_specie1[d]}-DB_{list_specie1[i]}-query.bln6"
            outputhit=f"{output_blast}/{list_specie1[i]}_{list_specie1[d]}_busco_genes.txt"
            if outputhit not in list_file:
                list_file.append(outputhit)
            my_analysis=ComparisonBlast()
            my_analysis.besthit(blast1,blast2,outputhit)
        d+=1
    i+=1

#print(list_file)


list_pid=[]

hs_file = open(f"{output_blast}/Busco_genes_HS.txt", 'a')
for file in list_file:
    splt=file.split("/")
    species=splt[3]
    splt2=species.split("_")
    sp1=splt2[0]
    sp2=splt2[1]
    with open(file, "r") as filin:
        for line in filin:
            line=line[:-1]
            line=line.split("\t")
            pid=float(line[2])
            if pid not in list_pid:
                list_pid.append(pid)
    filin.close()
    pid_file=open(f"{file}_busco_genes_pid.txt", 'w')
    for val in list_pid:
        pid_file.write(f"{val}\n")
    pid_file.close()
    Q1=np.percentile(list_pid, 25, interpolation = 'midpoint')
    Q3=np.percentile(list_pid, 75, interpolation = 'midpoint')
    IQR = Q3 - Q1
    HS = Q3 + (IQR / 2)
    hs_file.write(f"{sp1}\t{sp2}\t{Q1}\t{Q3}\t{HS}\n")

hs_file.close()

###############################################################################################################
######################################### Building HS dictionnary #############################################
###############################################################################################################

HS_dict={}
with open(f'{output_blast}/Busco_genes_HS.txt', 'r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split('\t')
        sp1=line[0]
        sp2=line[1]
        sp12=f"{sp1}-{sp2}"
        sp21=f"{sp2}-{sp1}"
        HS=float(line[4])
        if sp12 not in HS_dict:
            HS_dict[sp12]=HS
        if sp21 not in HS_dict:
            HS_dict[sp21]=HS
    filin.close()

#print(HS_dict)

###############################################################################################################
######################################### Selection of genes validating HS ####################################
###############################################################################################################

Gene_HS=open(f"{GENE_annot}.temp", 'w')
with open(f"{GENE_annot}", 'r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split("\t")
        sp1=line[0]
        sp2=line[1]
        sp=f"{sp1}-{sp2}"
        gene=line[2]
        PID=float(line[3])
        sp1_id_scfd=line[4]
        sp2_id_scfd=line[5]
        if sp in HS_dict.keys():
            Gene_HS.write(f"{sp1}\t{sp2}\t{gene}\t{PID}\t{sp1_id_scfd}\t{sp2_id_scfd}\t{HS_dict.get(sp)}\n")
    filin.close()
Gene_HS.close()

Gene_HSvalidation=open(GENE_validated,'w')
out_id_gene=open(GENE_ID, 'w')
id_gene_list=[]
with open(f"{GENE_annot}.temp", "r") as filin:
    for line in filin:
        line = line [:-1]
        line = line.split("\t")
        sp1=line[0]
        sp2=line[1]
        gene=line[2]
        PID=float(line[3])
        sp1_id_scfd=line[4]
        sp2_id_scfd=line[5]
        HS=float(line[6])
        if PID >= HS:
            print(f"{gene} validates HS")
            Gene_HSvalidation.write(f"{sp1}\t{sp2}\t{gene}\t{PID}\t{sp1_id_scfd}\t{sp2_id_scfd}\t{HS}\n")
            if sp1_id_scfd not in id_gene_list:
                id_gene_list.append(sp1_id_scfd)
            if sp2_id_scfd not in id_gene_list:
                id_gene_list.append(sp2_id_scfd)
        else:
            print(f"{gene} not validates HS")
    filin.close()

Gene_HSvalidation.close()

for id in id_gene_list:
    out_id_gene.write(id)
    out_id_gene.close()

getSeq(scfd_fasta,GENE_ID,GENE_fasta)

###############################################################################################################
######################################### Selection of TEs validating HS ######################################
###############################################################################################################

TE_HS=open(f"{TE_annot}.temp", 'w')
with open(f"{TE_annot}", 'r') as filin:
    for line in filin:
        line=line[:-1]
        line=line.split("\t")
        sp1=line[0]
        sp2=line[1]
        sp=f"{sp1}-{sp2}"
        te=line[2]
        PID=float(line[3])
        sp1_id_scfd=line[4]
        sp2_id_scfd=line[5]
        if sp in HS_dict.keys():
            TE_HS.write(f"{sp1}\t{sp2}\t{te}\t{PID}\t{sp1_id_scfd}\t{sp2_id_scfd}\t{HS_dict.get(sp)}\n")
    filin.close()
TE_HS.close()

TE_HSvalidation=open(TE_validated,'w')
out_id_te=open(TE_ID, 'w')
id_te_list=[]
with open(f"{TE_annot}.temp", "r") as filin:
    for line in filin:
        line = line [:-1]
        line = line.split("\t")
        sp1=line[0]
        sp2=line[1]
        te=line[2]
        PID=float(line[3])
        sp1_id_scfd=line[4]
        sp2_id_scfd=line[5]
        HS=float(line[6])
        if PID >= HS:
            print(f"{te} validates HS")
            TE_HSvalidation.write(f"{sp1}\t{sp2}\t{te}\t{PID}\t{sp1_id_scfd}\t{sp2_id_scfd}\t{HS}\n")
            if sp1_id_scfd not in id_te_list:
                id_te_list.append(sp1_id_scfd)
            if sp2_id_scfd not in id_te_list:
                id_te_list.append(sp2_id_scfd)
        else:
            print(f"{te} not validates HS")
    filin.close()

TE_HSvalidation.close()

for id in id_te_list:
    out_id_te.write(id)
    out_id_te.close()

getSeq(scfd_fasta,TE_ID,TE_fasta)


###############################################################################################################
######################################### Deleting .temp files ################################################
###############################################################################################################

if os.path.isfile(f"{Gene_annot}.temp"):
  os.remove(f"{Gene_annot}.temp")
  print(f"{Gene_annot}.temp has been deleted")
else:
  print(f"{Gene_annot}.temp does not exist")

if os.path.isfile(f"{TE_annot}.temp"):
  os.remove(f"{TE_annot}.temp")
  print(f"{TE_annot}.temp has been deleted")
else:
  print(f"{TE_annot}.temp does not exist")

if os.path.isfile(f"{Unknown_annot}.temp"):
  os.remove(f"{Unknown_annot}.temp")
  print(f"{Unknown_annot}.temp has been deleted")
else:
  print(f"{Unknown_annot}.temp does not exist")

###############################################################################################################
###############################################################################################################
###############################################################################################################

print(f"Test of High Similarity done")
