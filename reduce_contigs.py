#!/usr/bin/python -tt

from Bio import SeqIO
from mean_read import *


def big_small_contigs(filein,fileout1,fileout2):
	handle = open(fileout1,'w')
	handle2 = open(fileout2,'w')
	for records in SeqIO.parse(filein,"fasta"):
		temp = str(records.description)
		temp2 = temp.split('_')
		length = int(temp2[3])
		if length < 500:
			handle.write(">"+str(records.description)+"\n")
			handle.write(str(records.seq)+"\n")
		else:
			handle2.write(">"+str(records.description)+"\n")
			handle2.write(str(records.seq)+"\n")
			
			
#big_small_contigs("spades_out/contigs.fasta",'petit.txt','grand.txt')

def blast_small_big(petit,grand,grand_db,blastout,fileout):
	os.system("makeblastdb -in %s -dbtype 'nucl' -out %s"%(grand, grand_db))
	os.system("blastn -max_target_seqs 2 -db %s -query %s -outfmt 5 -out %s -evalue 0.001"%(grand_db, petit, blastout))
	xml(blastout,fileout)
	
#blast_small_big('petit.txt','grand.txt','grand_db','blast_out_small_big.xml','sortie_petit_grand_blast_apres_xml.txt')