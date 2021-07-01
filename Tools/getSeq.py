#!/usr/bin/python3
# -*-coding:Utf-8 -*

import os
import sys

from Bio import SeqIO
import argparse


parser=argparse.ArgumentParser()

parser.add_argument('-f',dest='fasta_file',type=str,required=True, help='Define the fasta file')
parser.add_argument('-l',dest='IDs_list',type=str,required=True, help='Define the list of IDs')
parser.add_argument('-o',dest='output_file',type=str,required=True, help='Define the name of output file')

args=parser.parse_args()

fasta_file=args.fasta_file
IDs_list=args.IDs_list
output_file=args.output_file

IDs=set()
with open(IDs_list) as file:
	for line in file:
		line=line.strip()
		if line != "":
			IDs.add(line)

fasta_seq=SeqIO.parse(open(fasta_file),'fasta')
with open(output_file, 'w') as filout:
	for seq in fasta_seq:
		if seq.id in IDs:
			SeqIO.write([seq],filout,"fasta")
