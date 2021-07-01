#!/usr/bin/python3
# -*-coding:Utf-8 -*
from Bio import SeqIO
import sys
import os
import argparse


parser=argparse.ArgumentParser()

parser.add_argument('-f',dest='input_file',type=str,required=True,help='Define the fasta file to analyse')

args=parser.parse_args()
fastafile=args.input_file

for seq_record in SeqIO.parse(str(fastafile),"fasta"):
	output_line = '%i\t%s' %(len(seq_record), seq_record.id)
	print(output_line)

