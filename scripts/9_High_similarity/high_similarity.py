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


reciprocal=Blast()
reciprocal.makeblastdb(db,blast_path)
reciprocal.blastn(db,query,output,threads,blast_path)

###########################################################################################################################################
###########################################################################################################################################
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
            print(blast1)
            print(blast2)
            print(outputhit)
            print(f"{list_specie1[i]}\t{list_specie1[d]}")
            print(f"{i} and {d}")
            my_analysis=ComparisonBlast()
            my_analysis.besthit(blast1,blast2,outputhit)
        d+=1
    i+=1

print(list_file)


list_pid=[]

for file in list_file:
    splt=file.split("/")
    species=splt[7]
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
    hs_file = open(f"{output_blast}/Busco_genes_HS_{sp1}_{sp2}.txt", 'w')
    hs_file.write(f"{sp1}\t{sp2}\t{Q1}\t{Q3}\t{HS}\n")


os.system(f"cat {output_blast}/Busco_genes_HS_*.txt")

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
