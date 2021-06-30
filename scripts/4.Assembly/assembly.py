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

class Assembler():

    def cleanTR(self,prinseq_path,kmerreads,good_reads):
        self.prinseq_path=prinseq_path
        self.kmerreads=kmerreads
        self.good_reads=good_reads
        """
        Function allowing to remove low complexity sequences in your data.
        """
        os.system(f"{prinseq_path}/prinseq-lite.pl -fasta {kmerreads} -out_format 1 -lc_method dust -lc_threshold 10 -out_good {good_reads}")

    def convert_pairedend(self,prog_dir,good_reads,reads):
        self.prog_dir=prog_dir
        self.kmerreads=good_reads
        self.reads=reads
        """
        Get the paired end for the assemby
        """
        os.system('''grep ">" '''+good_reads+'''.fasta | sed 's.>..g' | cut -f 1 -d '_' | sort -u | awk '{print $0"""_1"}' > '''+good_reads+'''_1.id''')
        os.system('''grep ">" '''+good_reads+'''.fasta | sed 's.>..g' | cut -f 1 -d '_' | sort -u | awk '{print $0"""_2"}' > '''+good_reads+'''_2.id''')
        os.system(f"{prog_dir}/Pipeline_HT/Tools/getSeq.py -f {reads} -l {good_reads}_1.id -o {good_reads}_1.fa")
        os.system(f"{prog_dir}/Pipeline_HT/Tools/getSeq.py -f {reads} -l {good_reads}_2.id -o {good_reads}_2.fa")
        return

    def spades(self,spades_path,threads,output,fa1,fa2,k):
        self.threads=threads
        self.spades_path=spades_path
        self.output=output
        self.fa1=fa1
        self.fa2=fa2
        self.k=k
        """
        Function to aplly the assembly with SPAdes software
        """
        return os.system(f"{spades_path}/spades.py -t {threads} -o {output} -1 {fa1} -2 {fa2} --only-assembler -k {k}")

    cleanTR=classmethod(cleanTR)
    convert_pairedend=classmethod(convert_pairedend)
    spades=classmethod(spades)

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

prog_dir=param[0]
kmerreads=param[1]
good_reads=param[2]
reads=param[3]
threads=param[4]
output=param[5]
fa1=param[6]
fa2=param[7]
k=param[8]
prinseq_path=param[9]
spades_path=param[10]

assembly=Assembler()
assembly.cleanTR(prinseq_path,kmerreads,good_reads)
assembly.convert_pairedend(prog_dir,good_reads,reads)
assembly.spades(spades_path,threads,output,fa1,fa2,k)

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################


print(f"Job done: {fa1} and {fa2} have been successfully assembled")

end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
