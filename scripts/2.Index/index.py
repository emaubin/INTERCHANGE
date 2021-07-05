#!/usr/bin/python3
# -*-coding:Utf-8 -*


import os, sys
import subprocess
import time, datetime
import pandas as pd
from itertools import permutations
import argparse
import glob

class Index():


    def suffixerator(self,gt_path,fa,suffixindex):
        self.gt_path=gt_path
        self.fa=fa
        self.suffixindex=suffixindex
        return os.system(f'{gt_path}/gt suffixerator -db {fa} -indexname {suffixindex}_suffixindex -tis -suf -lcp -des -ssp -sds -dna -parts 6')


    def mkindex(self,gt_path,suffixindex,kmersize,indexmer):
        self.gt_path=gt_path
        self.suffixindex=suffixindex
        self.kmersize=kmersize
        self.indexmer=indexmer
        return os.system(f'{gt_path}/gt tallymer mkindex -esa {suffixindex}_suffixindex -mersize {kmersize} -indexname {indexmer}_{kmersize}mers_index -minocc 1 -counts -pl')

    suffixerator=classmethod(suffixerator)
    mkindex=classmethod(mkindex)


start=time.time()


parser=argparse.ArgumentParser()
parser.add_argument('-p',dest='parameters',type=str,required=True, help='Define parameters file')
args=parser.parse_args()
paramfile=args.parameters

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

# Open paramters file
with open(paramfile, 'r') as param:
    param= param.read().split("\n")

fa=param[0]
suffixindex=param[1]
indexmer=param[2]
kmersize=param[3]
gt_path=param[4]

# Use the functions to creat indexes
index=Index()
index.suffixerator(gt_path,fa,suffixindex)
index.mkindex(gt_path,suffixindex,kmersize,indexmer)

fileList = glob.glob(f"{suffixindex}_suffixindex.*")
#Iterate over the list of fileList and remove each file
for filePath in fileList:
    try:
        os.remove(filePath)
        print("This file has been successfully deleted : ", filePath)
    except:
        print("Error while deleting file : ", filePath)


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

print(f"Job done: k-mers index has been successfully generated for the fasta file {fa}")

end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
