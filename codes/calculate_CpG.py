### This script counts A, T, G, C, CpG , TpG, and CpA sites of whole sequences, as well as calculates expected CpG, expected TpG, expected CpA, CpGo/e, TpGo/e CpAo/e, and then eta=(TpGo/e + CpAo/e)/(2 * CpGo/e)### 
### Usage: python script.py input.fasta > output.txt ###

#!/usr/bin/python

from __future__ import division
from Bio import SeqIO
import sys

print('ID Length Acount Ccount Gcount Tcount CGo_e TGo_e CAo_e eta\n')
fasta = sys.argv[1]

f = open(fasta, "rU")
for record in SeqIO.parse(f, "fasta"):
	sq = str(record.seq).upper()
	Acount = 0
	Ccount = 0
	Gcount = 0
	Tcount = 0
	CGexp = 0
	TGexp = 0
	CAexp = 0
	CGcount = 0
	TGcount = 0
	CAcount = 0
	CGo_e = 0
	TGo_e = 0
	CAo_e = 0
	eta = 0
	for i, nt in enumerate(sq):
		try:
			if sq[i] == "A":
				Acount += 1
			if sq[i] == "C":
				Ccount += 1
			if sq[i] == "G":
				Gcount += 1
			if sq[i] == "T":
				Tcount += 1     
			if sq[i] == "C" and sq[i+1] == "G":
				CGcount += 1
			if sq[i] == "T" and sq[i+1] == "G":
				TGcount += 1
			if sq[i] == "C" and sq[i+1] == "A":
				CAcount += 1
		except:
			continue
	CGexp = (Ccount*Gcount)/len(record.seq)
	TGexp = (Tcount*Gcount)/len(record.seq)
	CAexp = (Ccount*Acount)/len(record.seq)
	try:
		CGo_e = CGcount/CGexp
	except ZeroDivisionError: 
		CGo_e = float('Inf')
	try:
		TGo_e = TGcount/TGexp
	except ZeroDivisionError: 
		TGo_e = float('Inf')
	try:
		CAo_e = CAcount/CAexp
	except ZeroDivisionError: 
		CAo_e = float('Inf')
	try:
		eta = (TGo_e+CAo_e)/(CGo_e+CGo_e)
	except ZeroDivisionError: 
		eta= "NA"
	print record.id, len(record.seq), Acount, Ccount, Gcount, Tcount, CGo_e, TGo_e, CAo_e, eta
f.close()
