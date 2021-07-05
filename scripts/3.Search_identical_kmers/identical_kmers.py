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

class KmerSearch():

    def search(self,gt_path,indexmer,query,output):
        self.gt_path=gt_path
        self.indexmer=indexmer
        self.query=query
        self.output=output
        return os.system(f'{gt_path}/gt tallymer search -tyr {indexmer} -q {query} -output qseqnum qpos counts sequence > {output}.tallymer')

    def coord(self,output,kmersize,kmercoord):
        self.output=output
        self.kmersize=kmersize
        self.kmercoord=kmercoord
        return os.system('''sed 's/+//g' '''+output+'''.tallymer | awk 'OFS="\t"{print $1,$2,$2+'''+kmersize+'''}'  | mergeBed -i -  | awk 'OFS="\t"{print $0,$3-$2}' > '''+kmercoord)

    def getSeq(self,query,output_filter,kmerreads):
        self.query=query
        self.output_filter=output_filter
        self.kmerreads=kmerreads
        return os.system(f'{prog_dir}/INTERCHANGE-V.1.0/Tools/getSeq.py -f {query} -l {output_filter}_ID_kmers_list.txt -o {kmerreads}')

    def filter(self,query,kmercoord,threshold,readsize,output_filter):
        self.query=query
        self.kmercoord=kmercoord
        self.threshold=threshold
        self.readsize=readsize
        self.output_filter=output_filter
        dico_Nread_kmersize={}
        list_num_read=[]
        with open(query,'r') as filin:
            list_ID=[]
            for line in filin:
                if '>' in line:
                    id=line[1:-1]
                    list_ID.append(id)
        filin.close()
        with open(kmercoord,'r') as filin:
            for line in filin:
                line=line.split('\t')
                Nread=line[0]
                kmerSize=int(line[3])
                if Nread not in dico_Nread_kmersize:
                    dico_Nread_kmersize[Nread]=[int(kmerSize)]
                else:
                    dico_Nread_kmersize[Nread].append(int(kmerSize))
        filin.close()
        for key, value in dico_Nread_kmersize.items():
            kmerLength=sum(value)
            kmerCov=(kmerLength/int(readsize))*100
            if str(kmerCov)>=threshold:
                list_num_read.append(int(key))
        list_IDkmer=[]
        counter1=0
        counter2=0
        while counter1<len(list_num_read):
            numread=list_num_read[counter2]
            IDkmer=list_ID[numread]
            list_IDkmer.append(IDkmer)
            counter1+=1
            counter2+=1
        with open(f'{output_filter}_ID_kmers_list.txt','w') as filout:
            for item in list_IDkmer:
                filout.write('%s\n' %(item))
        filout.close()


    search=classmethod(search)
    coord=classmethod(coord)
    getSeq=classmethod(getSeq)
    filter=classmethod(filter)


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

indexmer=param[0]
query=param[1]
output=param[2]
kmersize=param[3]
kmercoord=param[4]
threshold=param[5]
readsize=param[6]
output_filter=param[7]
kmerreads=param[8]
prog_dir=param[9]
gt_path=param[10]

search=KmerSearch()
search.search(gt_path,indexmer,query,output)
search.coord(output,kmersize,kmercoord)
search.filter(query,kmercoord,threshold,readsize,output_filter)
search.getSeq(query,output_filter,kmerreads)

try:
    os.remove(f"{output}.tallymer")
    print(f"Tallymer file has been successfully deleted : {output}.tallymer")

except:
    print(f"Error while deleting file : {output}.tallymer")

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

print(f"")

end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
