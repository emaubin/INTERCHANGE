#!/usr/bin/python3
# -*-coding:Utf-8 -*


import os, sys
import os.path
import subprocess
import time, datetime
import pandas as pd
from itertools import permutations
import argparse
import shutil

start=time.time()

class OrthoGenes():

    def mapping(self,minimap2_path,tmini,target,fq1,fq2,out_mapped,samtools_path,mapped_id):
        self.minimap2_path=minimap2_path
        self.tmini=tmini
        self.target=target
        self.fq1=fq1
        self.fq2=fq2
        self.out_mapped=out_mapped
        self.mapped_id=mapped_id
        self.samtools_path=samtools_path
        os.system(f"{minimap2_path}/minimap2 -t {tmini} -ax sr {target} {fq1} {fq2} | samtools view -b -F 4 -@ 2 > {out_mapped}")
        os.system(samtools_path+'''/samtools view '''+out_mapped+''' | awk '{print $1"\t"$3}' > '''+mapped_id)
        return

    def get_reads(self,prog_dir,mapped_id,reads_busco_id1,reads_busco_id2,reads_sp,reads_busco1,reads_busco2):
        self.mapped_id=mapped_id
        self.reads_busco_id1=reads_busco_id1
        self.reads_busco_id2=reads_busco_id2
        self.reads_sp=reads_sp
        self.reads_busco1=reads_busco1
        self.reads_busco2=reads_busco2
        self.prog_dir=prog_dir

        dict_ID={}
        with open(mapped_id, "r") as filin:
            for line in filin:
                line=line[:-1]
                line=line.split("\t")
                id=line[0]
                splt=line[1].split("#")
                gene=splt[1]
                if gene not in dict_ID:
                    dict_ID[gene]=[id]
                else:
                    dict_ID[gene].append(id)
            filin.close()

        for k, v in dict_ID.items():
            print(k)
            with open(reads_busco_id1, 'a') as filout:
                for val in v:
                    filout.write(f'{val}_1\n')
                filout.close()
            with open(reads_busco_id2, 'a') as filout:
                for val in v:
                    filout.write(f'{val}_2\n')
                filout.close()

        os.system(f"{prog_dir}/INTERCHANGE-V.1.0/Tools/getSeq.py -f {reads_sp} -l {reads_busco_id1} -o {reads_busco1}")
        os.system(f"{prog_dir}/INTERCHANGE-V.1.0/Tools/getSeq.py -f {reads_sp} -l {reads_busco_id2} -o {reads_busco2}")

        return

    def gene_assembly(self,spades_path,tspades,output_Assembly,reads_busco1,reads_busco2,k,new_scaffolfds,species):
        self.spades_path=spades_path
        self.tspades=tspades
        self.output_Assembly=output_Assembly
        self.reads_busco1=reads_busco1
        self.reads_busco2=reads_busco2
        self.k=k
        self.new_scaffolfds=new_scaffolfds

        os.system(f"{spades_path}/spades.py -t {tspades} -o {output_Assembly} -1 {reads_busco1} -2 {reads_busco2} --only-assembler -k {k}")
        scaffolds_busco=f"{output_Assembly}/scaffolds.fasta"
        if os.path.isfile(scaffolds_busco):
            os.system('''sed 's/>/>'''+species+'''_/g' '''+scaffolds_busco+''' > '''+new_scaffolfds)
        else:
            print ("File not exist")


    def blastn(self,db,query,output_blast,tblast,blast_path,prog_dir):
        self.db=db
        self.query=query
        self.output_blast=output_blast
        self.tblast=tblast
        self.blast_path=blast_path
        self.prog_dir=prog_dir

        os.system(f"{blast_path}/blastn -task blastn -query {query} -db {db} -outfmt 6 -evalue 1e-20 -out {output_blast}.bln6 -num_threads {tblast}")
        os.system(f"sort -k12,12nr {output_blast}.bln6 > {output_blast}.bestsort")
        besthit_dict={}
        with open(f"{output_blast}.bestsort", 'r') as filin:
            for line in filin:
                line = line[:-1]
                line=line.split("\t")
                id_scfd=line[0]
                splt=line[1].split("#")
                id_gene=splt[1]

                if id_gene not in besthit_dict:
                    besthit_dict[id_gene]=id_scfd
            filin.close()
        besthit_id_file=open(f"{output_blast}.id" ,'w')
        for k, v in besthit_dict.items():
            besthit_id_file.write(f"{v}\n")
        besthit_id_file.close()
        os.system(f"{prog_dir}/INTERCHANGE-V.1.0/Tools/getSeq.py -f {query} -l {output_blast}.id -o {output_blast}_vf.fa")
        return

    mapping=classmethod(mapping)
    get_reads=classmethod(get_reads)
    gene_assembly=classmethod(gene_assembly)

parser=argparse.ArgumentParser()
parser.add_argument('-p',dest='parameters',type=str,required=True, help='Define parameters file')
args=parser.parse_args()
paramfile=args.parameters


###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

with open(paramfile, 'r') as param:
    param= param.read().split("\n")

minimap2_path = param[0]
tmini= param[1]
target=db= param[2]
fq1= param[3]
fq2= param[4]
out_mapped= param[5]
samtools_path= param[6]
mapped_id= param[7]
reads_busco_id1= param[8]
reads_busco_id2= param[9]
reads_sp= param[10]
reads_busco1= param[11]
reads_busco2= param[12]
spades_path= param[13]
tspades= param[14]
output_Assembly= param[15]
k= param[16]
new_scaffolfds= query = param[17]
species= param[18]
prog_dir= param[19]
output_blast=param[20]
blast_path=param[21]
tblast=param[22]

busco_find=OrthoGenes()
#busco_find.mapping(minimap2_path,tmini,target,fq1,fq2,out_mapped,samtools_path,mapped_id)
#busco_find.get_reads(prog_dir,mapped_id,reads_busco_id1,reads_busco_id2,reads_sp,reads_busco1,reads_busco2)
#busco_find.gene_assembly(spades_path,tspades,output_Assembly,reads_busco1,reads_busco2,k,new_scaffolfds,species)
busco_find.blastn(db,query,output_blast,tblast,blast_path,prog_dir)

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

print(f"Job done : {output_Assembly} : Busco gene have been successfully assembled")
end=time.time()
time_sec=end-start
time_hour=str(datetime.timedelta(seconds=time_sec))
print(f"Job execution time : {time_hour}")
