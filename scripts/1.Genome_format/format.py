#!/usr/bin/python3
# -*-coding:Utf-8 -*


import os, sys
import subprocess
import time, datetime
import pandas as pd
from itertools import permutations
import argparse
import shutil

start=time.time()

class Format():

    def uncompress(self,thread,file):
        self.thread=thread
        self.file=file
        return os.system(f"unpigz -p {thread} {file}")

    def compress(self,thread,file):
        self.thread=thread
        self.file=file
        return os.system(f"pigz -p {thread} {file}")

    def header(self,fq1,fq2,outfq1,outfq2):
        self.fq1=fq1
        self.fq2=fq2
        self.outfq1=outfq1
        self.outfq2=outfq2
        os.system('''awk -F " " '{if ($1~/^@/){print $1"_1"};if ($1!~/^@/){print $0}}' '''+fq1+''' > '''+outfq1)
        os.system('''awk -F " " '{if ($1~/^@/){print $1"_2"};if ($1!~/^@/){print $0}}' '''+fq2+''' > '''+outfq2)
        return

    def fqtofa(self,outfq1,outfq2,fa):
        self.outfq1=outfq1
        self.outfq2=outfq2
        self.fa=fa
        os.system('''sed -n '1~4s/^@/>/p;2~4p' '''+outfq1+''' > '''+fa+'''_reads_1.fa''')
        os.system('''sed -n '1~4s/^@/>/p;2~4p' '''+outfq2+''' > '''+fa+'''_reads_2.fa''')
        os.system(f"cat {fa}_reads_1.fa {fa}_reads_2.fa > {fa}_reads.fa")
        return

    uncompress=classmethod(uncompress)
    header=classmethod(header)
    compress=classmethod(compress)
    fqtofa=classmethod(fqtofa)


parser=argparse.ArgumentParser()
parser.add_argument('-p',dest='parameters',type=str,required=True, help='Define parameters file')
args=parser.parse_args()
paramfile=args.parameters

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

with open(paramfile, 'r') as param:
    param= param.read().split("\n")


fq1=param[0]
fq2=param[1]
outfq1=param[2]
outfq2=param[3]
fa=param[4]
thread=param[5]


reformat=Format()

file=f'{fq1}.gz'
reformat.uncompress(thread,file)
file=f'{fq2}.gz'
reformat.uncompress(thread,file)


reformat.header(fq1,fq2,outfq1,outfq2)
reformat.fqtofa(outfq1,outfq2,fa)
file=fq1
reformat.compress(thread,file)
file=fq2
reformat.compress(thread,file)
os.remove(outfq1)
os.remove(outfq2)
os.remove(f"{fa}_reads_1.fa")
os.remove(f"{fa}_reads_2.fa")

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

print(f"Job done : {fa}_reads.fa has been successfully generated")
end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
