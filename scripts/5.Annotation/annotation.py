#!/usr/bin/python3
# -*-coding:Utf-8 -*


import os, sys
import subprocess
import time, datetime
import pandas as pd
from itertools import permutations
import argparse

start=time.time()

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################


class Annotation():

    def fastaheader(self,species,fasta,assembly):
        self.species=species
        self.fasta=fasta
        self.assembly=assembly
        """
        Function allowing to rename the fastafile header
        """
        return os.system('''sed 's/^>/>'''+species+'''_/g' '''+fasta+''' > '''+assembly)

    def sizefilter(self,prog_dir,assembly):
        self.prog_dir=prog_dir
        self.assembly=assembly
        """
        Function allowing to remove scaffolds with length less than 300 bp
        """
        os.system(prog_dir+'''/INTERCHANGE-V.1.0/Tools/seqlength.py -f'''+assembly+''' | sort -k1,1nr | awk -F " " '{if ($1>=300){print $2}}' > '''+assembly+'''.IDsup300''')
        os.system(f"{prog_dir}/INTERCHANGE-V.1.0/Tools/getSeq.py -f {assembly} -l {assembly}.IDsup300 -o {assembly}_sup300.fa")
        return

    def te_annotation(self,prog_dir,diamond_path,assembly,species,output_Repbase,threads,db_path):
        self.prog_dir=prog_dir
        self.diamond_path=diamond_path
        self.assembly=assembly
        self.species=species
        self.output_Repbase=output_Repbase
        self.threads=threads
        self.db_path=db_path

        """
        Function allowing to blast all scaffolds against Repbase database and keep best hits
        """
        os.system(f"{diamond_path}/diamond blastx -q {assembly}_sup300.fa -d {db_path}/Databases/repbase/repbase19.06_aaSeq_cleaned_TE.no.Gypsy-8_PX -o {output_Repbase}/{species}_sup300_Repbase.blx6 -e 1e-06 -f 6 -p {threads}")
        os.system(f"sort -k12,12nr {output_Repbase}/{species}_sup300_Repbase.blx6 | sort -u -k1,1 > {output_Repbase}/{species}_sup300_Repbase.blx6.best")
        return

    def mcr_annotation(self,prog_dir,blast_path,assembly,species,output_mcr,threads,db_path):
        self.prog_dir=prog_dir
        self.blast_path=blast_path
        self.assembly=assembly
        self.species=species
        self.output_mcr=output_mcr
        self.threads=threads
        self.db_path=db_path

        os.system(f"{blast_path}/blastn -task blastn -db {db_path}/Databases/MCR_2021/MCR_database_GetOrganelleDB_mitochondrion#1#2_plastid#1#2#3#4_TIGR_rrna.fasta -query {assembly}_sup300.fa -outfmt 6 -evalue 1e-20 -num_threads {threads} -out {output_mcr}/{species}_scaffolds_DBMCR.bln6")
        os.system(f"sort -k12,12nr {output_mcr}/{species}_scaffolds_DBMCR.bln6 | sort -u -k1,1 > {output_mcr}/{species}_scaffolds_DBMCR.bln6.best")
        return

    def cdd_delta(self,prog_dir,diamond_path,assembly,output_cdd,threads,db_path):
        self.prog_dir=prog_dir
        self.diamond_path=diamond_path
        self.assembly=assembly
        self.output_cdd=output_cdd
        self.threads=threads
        self.db_path=db_path

        os.system(f"{diamond_path}/diamond blastx -q {assembly}_sup300.fa -d {db_path}/Databases/CDD_delta/cdd_delta_TEfree -o {output_cdd}/{species}_scaffolds_cddDelta.blx6 -e 1e-06 -f 6 -p {threads}")
        os.system(f"sort -k12,12nr {output_cdd}/{species}_scaffolds_cddDelta.blx6 | sort -u -k1,1 > {output_cdd}/{species}_scaffolds_cddDelta.blx6.best")
        return

    def otherDB(self,prog_dir,blast_path,assembly,species,output_otherDB,db_path):
        self.assembly=assembly
        self.blast_path=blast_path
        self.species=species
        self.output_otherDB=output_otherDB
        self.prog_dir=prog_dir
        self.threads=threads
        self.db_path=db_path

        os.system(f"{blast_path}/blastn -task blastn -db {db_path}/Databases/others/last.database.mito.cholor.e.coli -query {assembly}_sup300.fa -outfmt 6 -evalue 1e-20 -num_threads {threads} -out {output_otherDB}/{species}_scaffolds_otherDB.bln6")
        os.system(f"sort -k12,12nr {output_otherDB}/{species}_scaffolds_otherDB.bln6 | sort -u -k1,1 > {output_otherDB}/{species}_scaffolds_otherDB.bln6.best")
        return

    fastaheader=classmethod(fastaheader)
    sizefilter=classmethod(sizefilter)
    te_annotation=classmethod(te_annotation)
    mcr_annotation=classmethod(mcr_annotation)
    cdd_delta=classmethod(cdd_delta)
    otherDB=classmethod(otherDB)

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

species=param[0]
fasta=param[1]
assembly=param[2]
output_Repbase=param[3]
output_mcr=param[4]
output_cdd=param[5]
output_otherDB=param[6]
threads=param[7]
prog_dir=param[8]
diamond_path=param[9]
blast_path=param[10]
db_path=param[11]

annotation=Annotation()
annotation.fastaheader(species,fasta,assembly)
annotation.sizefilter(prog_dir,assembly)
annotation.te_annotation(prog_dir,diamond_path,assembly,species,output_Repbase,threads,db_path)
annotation.mcr_annotation(prog_dir,blast_path,assembly,species,output_mcr,threads,db_path)
annotation.cdd_delta(prog_dir,diamond_path,assembly,output_cdd,threads,db_path)
annotation.otherDB(prog_dir,blast_path,assembly,species,output_otherDB,db_path)


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
