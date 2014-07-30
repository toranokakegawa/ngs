#!/usr/bin/python -tt
import numpy as np
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO


#cfg = ConfigParser.ConfigParser()
#cfg.read('fichierdeconfig.cfg')
#rawfile = cfg.get('Filename', 'cle2')
#sizefile = cfg.get('Filename', 'cle42')
def size(f_in, f_out):
	os.system("cat %s | awk '{if(NR%%4==2) print length($1)}' > %s"%(f_in,f_out))
	sum=0
	count_lines=0
	with open(f_out,'r') as file:
		for line in file:
			sum=sum+int(line)
			count_lines=count_lines+1
	return sum/count_lines
	
	
	
def xml(filein,fileout):
	handle = open(fileout, "w")
	handle2 = open(filein,"r")
	for blast_record in NCBIXML.parse(handle2):
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				handle.write('****Alignment****\n')
				handle.write('identified genotype: '+alignment.title+'\n')
				handle.write('length: '+str(alignment.length)+'\n')
				handle.write('alignment length :'+str(hsp.align_length)+'\n')
				handle.write('score:'+str(hsp.score)+'\n')
				a = float(hsp.identities)
				b = float(hsp.align_length)
				handle.write('number of identities:'+str(a)+'\n')
				handle.write('identity percentage:'+str(round(((a/b)*100),3))+'%\n')
				handle.write('e-value:'+str(hsp.expect)+'\n')
				'''handle.write(hsp.query[:]+'\n')
				handle.write(hsp.match[:]+'\n')
				handle.write(hsp.sbjct[:]+'\n')'''
				handle.write(blast_record.query+'\n\n')
				
def xml2(filein,fileout):
	handle = open(fileout, "w")
	handle2 = open(filein,"r")
	list_hits = list()
	list_hit_al = list()
	list_hit_id = list()
	list_length_hits = list()
	list_coverage = list()
	list_hit_pid = list()
	list_hit_hsp = list()
	for blast_record in NCBIXML.parse(handle2):
		test_stop = False
		list_hits = list()
		list_hit_al = list()
		list_hit_id = list()
		list_length_hits = list()
		list_coverage = list()
		list_hit_pid = list()
		list_hit_hsp = list()
		for alignment in blast_record.alignments:
			length_hpv = float(alignment.length)
			count_hsp = 0
			list_al = list()
			list_id = list()
			for hsp in alignment.hsps:
				id = float(hsp.identities)
				al = float(hsp.align_length)
				#q_cover = (al/blast_record.query_length)*100
				#list_hits.append(alignment.title)
				list_al.append(al)
				list_id.append(id)
				#list_coverage.append(q_cover)
				#list_length_hits.append(alignment.length)
				count_hsp = count_hsp + 1
			#print list_al
			max_id = max(list_id)
			max_al = max(list_al)
			pid = (max_id/max_al)*100
			q_cover = (max_al/blast_record.query_length	)*100
			if (pid > 90):
				handle.write('****Alignment****\n')
				handle.write('identified genotype: '+alignment.title+'\n')
				handle.write('hpv sequence length: '+str(length_hpv)+'\n')
				handle.write('alignment length :'+str(max_al)+'\n')
				handle.write('numbers of identities :'+str(max_id)+'\n')
				handle.write('identity percentage:'+str(pid)+'%\n')
				handle.write('query coverage :'+str(q_cover)+'\n')
				handle.write('number of hsp :'+str(count_hsp)+'\n')
				handle.write('contig length :'+str(blast_record.query_length)+'\n')
				handle.write(blast_record.query+'\n\n')
				test_stop = True
				break
			else:
				list_hits.append(alignment.title)
				list_hit_id.append(max_id)
				list_hit_al.append(max_al)
				list_hit_pid.append(pid)
				list_coverage.append(q_cover)
				list_length_hits.append(length_hpv)
				list_hit_hsp.append(count_hsp)
		#if test_stop:
		#	break
		if len(list_hits)!=4:
			for i in range(len(list_hits),4):
				list_hits.append('empty')
		if len(list_hit_al)!=4:
			for i in range(len(list_hit_al),4):
				list_hit_al.append('empty')
		if len(list_hit_id)!=4:
			for i in range(len(list_hit_id),4):
				list_hit_id.append('empty')
		if len(list_length_hits)!=4:
			for i in range(len(list_length_hits),4):
				list_length_hits.append('empty')
		if len(list_coverage)!=4:
			for i in range(len(list_coverage),4):
				list_coverage.append('empty')
		if len(list_hit_pid)!=4:
			for i in range(len(list_hit_pid),4):
				list_hit_pid.append('empty')
		if len(list_hit_hsp)!=4:
			for i in range(len(list_hit_hsp),4):
				list_hit_hsp.append('empty')
		#print list_hits
		if not test_stop:
			if (list_coverage[0]<80):
				handle.write('****Alignment****\n')
				handle.write('identified genotype: new HPV genotype\n')
			else:
				handle.write('****Alignment****\n')
				handle.write('identified genotype: new HPV variant\n')
			handle.write('--------------------\n')
			handle.write('FIRST BLAST HIT: '+str(list_hits[0])+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[0])+'\n')
			handle.write('alignment length :'+str(list_hit_al[0])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[0])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[0])+'%\n')
			handle.write('number of hsp :'+str(list_hit_hsp[0])+'\n')
			handle.write('query coverage :'+str(list_coverage[0])+'\n')
			handle.write(blast_record.query+'\n')
			handle.write('--------------------\n')
			handle.write('SECOND BLAST HIT : '+str(list_hits[1])+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[1])+'\n')
			handle.write('alignment length :'+str(list_hit_al[1])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[1])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[1])+'%\n')
			handle.write('number of hsp :'+str(list_hit_hsp[1])+'\n')
			handle.write('query coverage :'+str(list_coverage[1])+'\n')
			handle.write('contig length :'+str(blast_record.query_length)+'\n')
			handle.write(blast_record.query+'\n\n')
			
			
		'''if not test_stop:
			if (list_coverage[0]<80):
				handle.write('****Alignment****\n')
				handle.write('identified genotype: new HPV genotype\n')
				handle.write('--------------------\n')
				handle.write('FIRST BLAST HIT: '+str(list_hits[0])+'\n')
				handle.write('hpv sequence length: '+str(list_length_hits[0])+'\n')
				handle.write('alignment length :'+str(list_hit_al[0])+'\n')
				handle.write('numbers of identities :'+str(list_hit_id[0])+'\n')
				handle.write('identity percentage:'+str(list_hit_pid[0])+'%\n')
				handle.write('number of hsp :'+str(list_hit_hsp[0])+'\n')
				handle.write('query coverage :'+str(list_coverage[0])+'\n')
				handle.write(blast_record.query+'\n')
				handle.write('--------------------\n')
				handle.write('SECOND BLAST HIT : '+str(list_hits[1])+'\n')
				handle.write('hpv sequence length: '+str(list_length_hits[1])+'\n')
				handle.write('alignment length :'+str(list_hit_al[1])+'\n')
				handle.write('numbers of identities :'+str(list_hit_id[1])+'\n')
				handle.write('identity percentage:'+str(list_hit_pid[1])+'%\n')
				handle.write('number of hsp :'+str(list_hit_hsp[1])+'\n')
				handle.write('query coverage :'+str(list_coverage[1])+'\n')
				handle.write('contig length :'+str(blast_record.query_length)+'\n')
				handle.write(blast_record.query+'\n\n')
				
			else:
				handle.write('****Alignment****\n')
				handle.write('identified genotype: new HPV variant\n')
				handle.write('--------------------\n')
				handle.write('FIRST BLAST HIT : '+str(list_hits[0])+'\n')
				handle.write('hpv sequence length: '+str(list_length_hits[0])+'\n')
				handle.write('alignment length :'+str(list_hit_al[0])+'\n')
				handle.write('numbers of identities :'+str(list_hit_id[0])+'\n')
				handle.write('identity percentage:'+str(list_hit_pid[0])+'%\n')
				handle.write('number of hsp :'+str(list_hit_hsp[0])+'\n')
				handle.write('query coverage :'+str(list_coverage[0])+'\n')
				handle.write('contig length :'+str(blast_record.query_length)+'\n')
				handle.write(blast_record.query+'\n')
				handle.write('--------------------\n')
				handle.write('SECOND BLAST HIT : '+str(list_hits[1])+'\n')
				handle.write('hpv sequence length: '+str(list_length_hits[1])+'\n')
				handle.write('alignment length :'+str(list_hit_al[1])+'\n')
				handle.write('numbers of identities :'+str(list_hit_id[1])+'\n')
				handle.write('identity percentage:'+str(list_hit_pid[1])+'%\n')
				handle.write('number of hsp :'+str(list_hit_hsp[1])+'\n')
				handle.write('query coverage :'+str(list_coverage[1])+'%\n')
				handle.write('contig length :'+str(blast_record.query_length)+'\n')
				handle.write(blast_record.query+'\n\n')'''

				
				
def type_hpv(typein):
	if ('type 2' in typein)	or ('type 6' in typein) or ('type 7' in typein) or ('type 10' in typein) or ('type 16' in typein) or ('type 18' in typein) or ('type 26' in typein) or ('type 32' in typein) or ('type 34' in typein) or ('type 53' in typein) or ('type 54' in typein) or ('type 61' in typein) or ('type 90' in typein):
		return 'alphapapillomavirus'
	if ('type 5' in typein)	or ('type 9' in typein) or ('type 49' in typein) or ('type 92' in typein) or ('type 96' in typein):
		return 'betapapillomavirus'
	if ('type 4' in typein)	or ('type 48' in typein) or ('type 50' in typein) or ('type 60' in typein) or ('type 88' in typein) or ('type 101' in typein) or ('type 109' in typein)	or ('type 112' in typein) or ('type 116' in typein) or ('type 121' in typein):
		return 'gammapapillomavirus'
	if ('type 1' in typein)	or ('type 63' in typein):
		return 'mupapillomavirus'
	if ('type 41' in typein):
		return 'nupapillomavirus'
	
				
				
#xml2('blast_spades_contigs_chosen_size.xml','taxo_spades_xml2.txt')		

def xmlBeta(filein,fileout):
	handle = open(fileout, "w")
	handle2 = open(filein,"r")
	list_hits = list()
	list_hit_al = list()
	list_hit_id = list()
	list_length_hits = list()
	list_coverage = list()
	list_hit_pid = list()
	list_hit_hsp = list()
	for blast_record in NCBIXML.parse(handle2):
		test_stop = False
		list_hits = list()
		list_hit_al = list()
		list_hit_id = list()
		list_length_hits = list()
		list_coverage = list()
		list_hit_pid = list()
		list_hit_hsp = list()
		for alignment in blast_record.alignments:
			length_hpv = float(alignment.length)
			count_hsp = 0
			list_al = list()
			list_id = list()
			for hsp in alignment.hsps:
				id = float(hsp.identities)
				al = float(hsp.align_length)
				#q_cover = (al/blast_record.query_length)*100
				#list_hits.append(alignment.title)
				list_al.append(al)
				list_id.append(id)
				#list_coverage.append(q_cover)
				#list_length_hits.append(alignment.length)
				count_hsp = count_hsp + 1
			#print list_al
			max_id = max(list_id)
			max_al = max(list_al)
			pid = (max_id/max_al)*100
			q_cover = (max_al/blast_record.query_length	)*100
			'''if (pid > 90):
				handle.write('****Alignment****\n')
				handle.write('identified genotype: '+alignment.title+'\n')
				handle.write('hpv sequence length: '+str(length_hpv)+'\n')
				handle.write('alignment length :'+str(max_al)+'\n')
				handle.write('numbers of identities :'+str(max_id)+'\n')
				handle.write('identity percentage:'+str(pid)+'%\n')
				handle.write('query coverage :'+str(q_cover)+'\n')
				handle.write('number of hsp :'+str(count_hsp)+'\n')
				handle.write('contig length :'+str(blast_record.query_length)+'\n')
				handle.write(blast_record.query+'\n\n')
				test_stop = True
				break
			else:'''
			list_hits.append(alignment.title)
			list_hit_id.append(max_id)
			list_hit_al.append(max_al)
			list_hit_pid.append(pid)
			list_coverage.append(q_cover)
			list_length_hits.append(length_hpv)
			list_hit_hsp.append(count_hsp)
		if len(list_hits)!=4:
			for i in range(len(list_hits),4):
				list_hits.append('empty')
		if len(list_hit_al)!=4:
			for i in range(len(list_hit_al),4):
				list_hit_al.append('empty')
		if len(list_hit_id)!=4:
			for i in range(len(list_hit_id),4):
				list_hit_id.append('empty')
		if len(list_length_hits)!=4:
			for i in range(len(list_length_hits),4):
				list_length_hits.append('empty')
		if len(list_coverage)!=4:
			for i in range(len(list_coverage),4):
				list_coverage.append('empty')
		if len(list_hit_pid)!=4:
			for i in range(len(list_hit_pid),4):
				list_hit_pid.append('empty')
		if len(list_hit_hsp)!=4:
			for i in range(len(list_hit_hsp),4):
				list_hit_hsp.append('empty')
		if (list_hit_pid[0] > 90):
			handle.write('\n\n****Alignment****\n')
			handle.write('identified genotype: '+str(list_hits[0])+'\n')
			handle.write('hpv genus: '+type_hpv(str(list_hits[0]))+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[0])+'\n')
			handle.write('alignment length :'+str(list_hit_al[0])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[0])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[0])+'%\n')
			handle.write('query coverage :'+str(list_coverage[0])+'\n')
			handle.write('number of hsp :'+str(list_hit_hsp[0])+'\n')
			handle.write('contig length :'+str(blast_record.query_length)+'\n')
			handle.write(blast_record.query+'\n')
			for i in range(1,len(list_hits)):
				list_hpv_blast_name = list()
				if list_hits[i] != 'empty':
					list_hpv_blast_name = str(list_hits[i]).split(' ')
					hpv_name = str(list_hpv_blast_name[1])
					for it in range(2,len(list_hpv_blast_name)):
						hpv_name = hpv_name+' '+str(list_hpv_blast_name[it])
					if not str(list_hits[0]) in hpv_name:
						handle.write('--------------------\n')
						handle.write('HIT NUMBER : '+str(i)+'\n')
						handle.write('hpv name : '+str(hpv_name)+'\n')
						handle.write('hpv genus: '+type_hpv(str(hpv_name))+'\n')
						handle.write('hpv sequence length: '+str(list_length_hits[i])+'\n')
						handle.write('alignment length :'+str(list_hit_al[i])+'\n')
						handle.write('numbers of identities :'+str(list_hit_id[i])+'\n')
						handle.write('identity percentage:'+str(list_hit_pid[i])+'%\n')
						handle.write('number of hsp :'+str(list_hit_hsp[i])+'\n')
						handle.write('query coverage :'+str(list_coverage[i])+'\n')
						handle.write('contig length :'+str(blast_record.query_length)+'\n')
						handle.write(blast_record.query+'\n')
		else:
			if (list_coverage[0]<80):
				handle.write('\n\n****Alignment****\n')
				handle.write('identified genotype: new HPV genotype\n')
			else:
				handle.write('\n\n****Alignment****\n')
				handle.write('identified genotype: new HPV variant\n')
			handle.write('--------------------\n')
			handle.write('FIRST BLAST HIT: '+str(list_hits[0])+'\n')
			handle.write('hpv genus: '+type_hpv(str(list_hits[0]))+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[0])+'\n')
			handle.write('alignment length :'+str(list_hit_al[0])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[0])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[0])+'%\n')
			handle.write('number of hsp :'+str(list_hit_hsp[0])+'\n')
			handle.write('query coverage :'+str(list_coverage[0])+'\n')
			handle.write('contig length :'+str(blast_record.query_length)+'\n')
			handle.write(blast_record.query+'\n')
			handle.write('--------------------\n')
			handle.write('SECOND BLAST HIT : '+str(list_hits[1])+'\n')
			handle.write('hpv genus: '+type_hpv(str(list_hits[1]))+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[1])+'\n')
			handle.write('alignment length :'+str(list_hit_al[1])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[1])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[1])+'%\n')
			handle.write('number of hsp :'+str(list_hit_hsp[1])+'\n')
			handle.write('query coverage :'+str(list_coverage[1])+'\n')
			handle.write('contig length :'+str(blast_record.query_length)+'\n')
			handle.write(blast_record.query+'\n')
			handle.write('--------------------\n')
			handle.write('THIRD BLAST HIT : '+str(list_hits[2])+'\n')
			handle.write('hpv genus: '+type_hpv(str(list_hits[2]))+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[2])+'\n')
			handle.write('alignment length :'+str(list_hit_al[2])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[2])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[2])+'%\n')
			handle.write('number of hsp :'+str(list_hit_hsp[2])+'\n')
			handle.write('query coverage :'+str(list_coverage[2])+'\n')
			handle.write('contig length :'+str(blast_record.query_length)+'\n')
			handle.write(blast_record.query+'\n')
			handle.write('--------------------\n')
			handle.write('FOURTH BLAST HIT : '+str(list_hits[3])+'\n')
			handle.write('hpv genus: '+type_hpv(str(list_hits[3]))+'\n')
			handle.write('hpv sequence length: '+str(list_length_hits[3])+'\n')
			handle.write('alignment length :'+str(list_hit_al[3])+'\n')
			handle.write('numbers of identities :'+str(list_hit_id[3])+'\n')
			handle.write('identity percentage:'+str(list_hit_pid[3])+'%\n')
			handle.write('number of hsp :'+str(list_hit_hsp[3])+'\n')
			handle.write('query coverage :'+str(list_coverage[3])+'\n')
			handle.write('contig length :'+str(blast_record.query_length)+'\n')
			handle.write(blast_record.query+'\n')
				
				
#xmlBeta('blast_spades_contigs_chosen_size.xml','xmlBeta/taxo_spades_xmlBeta.txt')				
				
'''#create a file containing all the viral genotype without duplicates
def coverage_taxo():
	handle=open('test_taxo_sans_doublon.txt','w')
	for blast_record in NCBIXML.parse(open('blast_spades_contigs_chosen_size.xml','r')):
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				handle.write(alignment.title+"\n")
	a=np.array([])
	with open('test_taxo_sans_doublon.txt','r') as file:
		for line in file:
			a = np.append(a,line)
	handle2 = open('test.txt','w')
	for i in a:
		if 'type' in i:
			b = i.split()
			b = b[1:]
			b = ' '.join(b)
			handle2.write(str(b)+"\n")
		else:
			handle2.write(i)
	with open('test.txt','r') as file:
		for line in file:
			a = np.append(a,line)
		b=set(a)
	handle3 = open('test2.txt','w')
	for z in b:
		handle3.write(z)'''
    	



#retrieve all the queries from the blast with their hits
'''def query():
	handle=open('query.txt','w')
	for blast_record in NCBIXML.parse(open('blast_spades_contigs_chosen_size.xml','r')):
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				handle.write(blast_record.query+" ")
				#handle.write(str(hsp.align_length)+" ")
				handle.write(alignment.title+"\n")'''
				


'''#create a file with all genotype (once) and the matching contigs with the coverage
def coverage():
	handle = open('profondeur.txt','w')
	cover_list = list()
	with open('test2.txt','r') as file:
		for line in file:
			with open('query.txt','r') as file2:		
				for line2 in file2:
					if line in line2:
						cover=line2.split()
						cover_list.append(cover[0])
				handle.write(">"+line+','.join(cover_list)+"\n")'''
	
	
			
#create a file containing all the viral genotype without duplicates
def coverage_taxo(doublon,temp1,temp2,fileout):
	handle=open(temp1,'w')
	for blast_record in NCBIXML.parse(open(doublon,'r')):
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				handle.write(alignment.title+"\n")
	a=np.array([])
	with open(temp1,'r') as file:
		for line in file:
			a = np.append(a,line)
	handle2 = open(temp2,'w')
	for i in a:
		if 'type' in i:
			b = i.split()
			b = b[1:]
			b = ' '.join(b)
			handle2.write(str(b)+"\n")
		else:
			handle2.write(i)
	'''with open(temp2,'r') as file:
		for line in file:
			a = np.append(a,line)'''
	handle4=open(temp2,'r')		
	a=handle4.readlines()
	b=set(a)
	handle3 = open(fileout,'w')
	for z in b:
		handle3.write(z)
	handle.close()
	handle2.close()
	handle3.close()	
	handle4.close()

#coverage_taxo('blast_spades_contigs_chosen_size.xml','temp1.txt','temp2.txt','sansdoublons.txt')
		
#retrieve all the queries from the blast with their hits
def query(filein,fileout):
	handle=open(fileout,'w')
	for blast_record in NCBIXML.parse(open(filein,'r')):
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				handle.write(blast_record.query+" ")
				handle.write(alignment.title+"\n")
	handle.close()
						
		
#create a file with all genotype (once) and the matching contigs with the coverage
def coverage(file1,file2,fileout):
	handle = open(fileout,'w')
	cover_list = list()
	with open(file1,'r') as file:
		for line in file:
			with open(file2,'r') as file1:		
				for line2 in file1:
					if line in line2:
						cover=line2.split()
						cover_list.append(cover[0])
				handle.write(">"+line+','.join(cover_list)+"\n")
	handle.close()
		
		
		
		
#create a file that holds only hpv with a type and the sorted contigs
def sort_contigs(filein,fileout):
	handle = open(fileout,'w')
	for records in SeqIO.parse(filein,'fasta'):
		if  'Human papillomavirus type' in records.description:
			temp = str(records.seq)
			temp2 = [x.strip() for x in temp.split(',')]
			temp = set(temp2)
			temp2 = sorted(temp)
			handle.write(">"+records.description+"\n")
			handle.write(str(temp2)+"\n")
	handle.close()
	
#sort_contigs('profondeur.txt','taxo_final.fasta')
		
		
def entire_genomes(filein,fileout,fileout2):
	handle = open(fileout,'w')
	liste=list()
	for records in SeqIO.parse(filein,'fasta'):
		temp = str(records.seq)
		temp2 = [x.strip() for x in temp.split(',')]
		count = 0
		for i in temp2:
			temp3 = i.split('_')
			length = int(temp3[3])
			if length > 6000:
				liste.append(i)
				count = count + 1
		if count > 1:
			handle.write(">"+records.description+"\n")
			handle.write(str(set(liste))+"\n")
	handle.close()
		
'''	handle3 = open(fileout,'r')
	handle2 = open(fileout2,'w')
	liste = handle3.readlines()
	ensemble = set(liste)
	for unique in ensemble:
		handle2.write(">"+unique)
	handle.close()
	handle2.close()
	handle3.close()'''
	
	
def hpv_phylo(bank,filein,fileout):
	handle = open(bank,'r')
	handle2 = open(filein,'r')
	handle3 = open(fileout,'w')
	liste = list()
	for records in SeqIO.parse(filein,'fasta'):
		liste.append(records.description)
	print liste
	for record in SeqIO.parse(bank,'fasta'):
		for item in liste:
			if item in record.description:
				handle3.write(">"+str(record.description)+"\n")
				handle3.write(str(record.seq)+"\n")
	handle.close()
	handle2.close()
	handle3.close()

def big_seq(filein,fileout):
	handle = open(fileout,'w')
	for record in SeqIO.parse(filein,'fasta'):
		if len(record.seq)>6000:
			handle.write(">"+record.description+"\n")
			handle.write(str(record.seq)+"\n")
	handle.close()
	
#big_seq('hpv_phylo.txt','hpv_grande_seq_id.txt')	

#entire_genomes('spades_genotype_coverage_no_duplicates.fasta','genome_entier.txt','genome_entier2.txt')
#hpv_phylo('viral_nuc.fasta','genome_entier.txt','hpv_phylo.txt')
		
		