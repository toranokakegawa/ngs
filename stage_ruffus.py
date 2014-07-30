#!/usr/bin/python -tt


import os
import re
import sys
import ConfigParser
from Bio import SeqIO
from Bio import GenBank
from xlrd import open_workbook
from switch import switch
from numberError import *
import pexpect
from ruffus import *
import shutil
#import pylab
import numpy
from mean_read import *
#from plot import *

cfg = ConfigParser.ConfigParser()
cfg.read('fichierdeconfig.cfg')
job1_input1 = cfg.get('Filename', 'cle3')#ICTV file
job1_output1 = cfg.get('Filename', 'cle4')#nucleotide output file
job1_output2 = cfg.get('Filename', 'cle4(0)')#protein output file 
job1_input2 = cfg.get('Filename', 'cle5')#genbank input file carrying the nucleotide sequences (used to create the data bank)
job1_input3 = cfg.get('Filename', 'cle5(0)')#genbank input file carrying the protein sequences (used to create the data bank)

job2_input1 = cfg.get('Filename', 'cle2')#phred function input
job2_output1 = cfg.get('Filename', 'cle2(0)')#phred function output/trim function input
job2_output2 = cfg.get('Filename', 'cle2(1)')#trim function output
job2_output3 = cfg.get('Filename', 'cle2(13)')#text file containing reads size


job3_output1 = cfg.get('Filename', 'cle9')#mapping output file (mapped reads)
job3_output2 = cfg.get('Filename', 'cle8')#mapping output file (unmapped reads)
job3_input1 = cfg.get('Filename', 'cle2(1)')#bowtie input file 
job3_input2 = cfg.get('Filename', 'cle2(12)')#bowtie input file converted in fasta
#job3_output1 = cfg.get('Filename', 'cle10')#mapping output file (index that will be created during the mapping)
#job3_output4 = cfg.get('Filename', 'cle12')#conversion output (fastq to fasta)/match
#job3_output5 = cfg.get('Filename', 'cle11')#conversion output (fastq to fasta)/unmatch

job4_input1 = cfg.get('Filename', 'cle4')#nucleotide input file, before quotes removal
job4_input2 = cfg.get('Filename', 'cle4(1)')#nucleotide input file after quotes removal
job4_input3 = cfg.get('Filename', 'cle4(0)')#protein input file, before quotes removal
job4_input4 = cfg.get('Filename', 'cle4(3)')#protein input file, after quotes removal
job4_output1 = cfg.get('Filename', 'cle13')#nucleotide database created (dbnuc output)
job4_output2 = cfg.get('Filename', 'cle33')#protein database created (dbprot output)
job4_key1 = cfg.get('data_input', 'cle3')#evalue of the blast
job4_output3 = cfg.get('Filename', 'cle15')#nuclotide blast output
job4_key2 = cfg.get('data_input', 'cle4')#blastn format output
job4_output4 = cfg.get('Filename', 'cle17')#awk function output that give the aligned reads with exact position
job4_input5 = cfg.get('Filename', 'cle11(2)')#blastn input
job4_output5 = cfg.get('Filename', 'cle30')#awk function output that gives the unaligned reads
job4_output6 = cfg.get('Filename', 'cle31')#blastx output
job4_output7 = cfg.get('Filename', 'cle32')#awk function output after blastx 
job4_output8 = cfg.get('Filename', 'cle34')#blastn output used to select reads regarding their size (>50)
job4_output9 = cfg.get('Filename', 'cle35')#blastx output used to select reads regarding their size (>50)
job4_output10 = cfg.get('Filename', 'cle20')#final blast output/assembly input


job5_output_dir = cfg.get('Filename', 'cle21')#spades output
#job5_output_subdir = cfg.get('Filename','cle21(1)')#spades subdir K-int
job5_output1 = cfg.get('Filename', 'cle6')#spades contigs directly stored in the spades output directory



job6_key1 = cfg.get('data_input','cle5')#ray kmer size
job6_output_dir = cfg.get('Filename', 'cle22')#ray output directory
job6_output1 = cfg.get('Filename', 'cle22(1)')#ray contigs 

job7_input1 = cfg.get('Filename', 'cle4(2)')#contigs assembled by spades stored in K-int file
job7_output1 = cfg.get('Filename', 'cle36')#mapping output (singletons) after spades
job7_output2 = cfg.get('Filename', 'cle37')#mapping output (singletons) after ray
job7_output3 = cfg.get('Filename', 'cle38')#the whole set of viral variants (singleton and blast positive) (spades)
job7_output4 = cfg.get('Filename', 'cle39')#the whole set of viral variants (singleton and blast positive) (ray)
job7_key1 = cfg.get('data_input', 'cle7')#blast taxo output format (xml)

job8_output1 = cfg.get('Filename', 'cle40')#spades taxonomic file (blast) 
job8_output2 = cfg.get('Filename', 'cle41')#ray taxonomic file (blast)
job8_output3 = cfg.get('Filename', 'cle42')#file containing the blast output after xml parsing[spades]
job8_output4 = cfg.get('Filename', 'cle43')#file containing the blast output after xml parsing[ray]

job9_output1 = cfg.get('Filename', 'cle44')#[sapdes]file containing the contigs with relevant sizes
job9_output2 = cfg.get('Filename', 'cle45')##[ray]file containing the contigs with relevant sizes
job9_output3 = cfg.get('Filename', 'cle46')#[spades]blast output (xml) from the chosen (by size) contigs
job9_output4 = cfg.get('Filename', 'cle47')#[ray]blast output (xml) from the chosen (by size) contigs
job9_output5 = cfg.get('Filename', 'cle48')#[spades]taxonomic designation from the chosen (by size) contigs
job9_output6 = cfg.get('Filename', 'cle49')#[ray]taxonomic designation from the chosen (by size) contigs
job9_output7 = cfg.get('Filename', 'cle60')#mapped reads against the selected contigs[spades]
job9_output8 = cfg.get('Filename', 'cle61')#unmapped reads against the selected contigs[spades]
job9_output9 = cfg.get('Filename', 'cle62')#mapped reads against the selected contigs[ray]
job9_output10 = cfg.get('Filename', 'cle63')#unmapped reads against the selected contigs[ray]

job10_input1 = cfg.get('Filename', 'cle84')#contigs size 100[spades]
job10_input2 = cfg.get('Filename', 'cle85')#contigs size 100-250[spades]
job10_input3 = cfg.get('Filename', 'cle86')#contigs size 250-500[spades]
job10_input4 = cfg.get('Filename', 'cle87')#contigs size 500-1000[spades]
job10_input5 = cfg.get('Filename', 'cle88')#contigs size 1000[spades]
job10_input6 = cfg.get('Filename', 'cle89')#contigs size 100[ray]
job10_input7 = cfg.get('Filename', 'cle90')#contigs size 100-250[ray]
job10_input8 = cfg.get('Filename', 'cle91')#contigs size 250-500[ray]
job10_input9 = cfg.get('Filename', 'cle92')#contigs size 500-1000[ray]
job10_input10 = cfg.get('Filename', 'cle93')#contigs size 1000[ray]
job10_input11 = job5_output_dir
job10_input12 = job5_output1
job10_input13 = job6_output_dir
job10_input14 = job6_output1
job10_output1 = cfg.get('Filename', 'cle64')
job10_output2 = cfg.get('Filename', 'cle69')
job10_output3 = cfg.get('Filename', 'cle65')
job10_output4 = cfg.get('Filename', 'cle70')
job10_output5 = cfg.get('Filename', 'cle66')
job10_output6 = cfg.get('Filename', 'cle71')
job10_output7 = cfg.get('Filename', 'cle67')
job10_output8 = cfg.get('Filename', 'cle72')
job10_output9 = cfg.get('Filename', 'cle68')
job10_output10 = cfg.get('Filename', 'cle73')
job10_output11 = cfg.get('Filename', 'cle74')
job10_output12 = cfg.get('Filename', 'cle79')
job10_output13 = cfg.get('Filename', 'cle75')
job10_output14 = cfg.get('Filename', 'cle80')
job10_output15 = cfg.get('Filename', 'cle76')
job10_output16 = cfg.get('Filename', 'cle81')
job10_output17 = cfg.get('Filename', 'cle77')
job10_output18 = cfg.get('Filename', 'cle82')
job10_output19 = cfg.get('Filename', 'cle78')
job10_output20 = cfg.get('Filename', 'cle83')
job10_output21 = cfg.get('Filename', 'cle94')
job10_output22 = cfg.get('Filename', 'cle95')
job10_output23 = cfg.get('Filename', 'cle96')
job10_output24 = cfg.get('Filename', 'cle97')
job10_output25 = cfg.get('Filename', 'cle98')
job10_output26 = cfg.get('Filename', 'cle99')
job10_output27 = cfg.get('Filename', 'cle100')
job10_output28 = cfg.get('Filename', 'cle101')
job10_output29 = cfg.get('Filename', 'cle102')
job10_output30 = cfg.get('Filename', 'cle103')
job10_output31 = cfg.get('Filename', 'cle104')
job10_output32 = cfg.get('Filename', 'cle105')
job10_output33 = cfg.get('Filename', 'cle106')
job10_output34 = cfg.get('Filename', 'cle107')
job10_output35 = cfg.get('Filename', 'cle108')
job10_output36 = cfg.get('Filename', 'cle109')
job10_output37 = cfg.get('Filename', 'cle110')
job10_output38 = cfg.get('Filename', 'cle111')
job10_output39 = cfg.get('Filename', 'cle112')
job10_output40 = cfg.get('Filename', 'cle113')
job10_output41 = cfg.get('Filename', 'cle114')





'''job10_input1 = job9_output3
job10_input2 = job9_output4
job10_output1 = cfg.get('Filename', 'cle50')#temporary file for the removal of duplicates in viral genotypes
job10_output2 = cfg.get('Filename', 'cle51')#temporary file for the removal of duplicates in viral genotypes
job10_output3 = cfg.get('Filename', 'cle52')#file containing the viral genotypes without duplicates after spades assembling
job10_output4 = cfg.get('Filename', 'cle53')#file containing the viral genotypes without duplicates after ray assembling

job11_input1 = job9_output3
job11_input2 = job9_output4
job11_output1 = cfg.get('Filename', 'cle54')#file with all the queries and the matching hits [spades]
job11_output2 = cfg.get('Filename', 'cle55')#file with all the queries and the matching hits [ray]

job12_input1 = job10_output3
job12_input2 = job10_output4
job12_input3 = job11_output1
job12_input4 = job11_output2
job12_output1 = cfg.get('Filename', 'cle56')#file with the viral genotypes and all the contigs [spades]
job12_output2 = cfg.get('Filename', 'cle57')#file with the viral genotypes and all the contigs [ray]


job13_input1 = job12_output1
job13_input2 = job12_output2
job13_output1 = cfg.get('Filename', 'cle58')#file with no duplicates in viral genotypes and no duplicates in the contigs [spades]
job13_output2 = cfg.get('Filename', 'cle59')#file with no duplicates in viral genotypes and no duplicates in the contigs [ray]'''



def check_file_exists(job1_input1,job1_input3,job1_input2,job1_output2,job1_output1):
	if not os.path.exists(job1_output1):
		return True, "Missing file %s" % job1_output1
	else:
		return False, "File %s exists" % job1_output1
	if not os.path.exists(job1_output2):
		return True, "Missing file %s" % job1_output2
	else:
		return False, "File %s exists" % job1_output2

def config():
	pexpect.run("python config.py")
	if(os.path.isfile("fichierdeconfig.cfg")):
		print("your configuration file is ready\n")
	else:
		print("the configuration file 'fichierdeconfig.cfg' doesn't exist\n")
		
		
def directory():
	name = raw_input("give the name of the directory in which you want your output files to be stored :  \n").lower()
	pexpect.run("mkdir %s"%(name))
	if(os.path.isdir(name)):
		handle = open(name+"/config2.py", "w")
		handle.write("#!/usr/bin/python -tt\nimport ConfigParser\n\n")
		handle.write("cfg = ConfigParser.ConfigParser()\ncfg.add_section('directory')\n\n")
		handle.write("S = 'directory'\n")
		handle.write("cfg.set(S, 'cle1', '"+name+"')\n\n")
		handle.write("cfg.write(open('fichierdeconfig2.cfg','w'))")
		handle.close()
		pexpect.run("python /Users/judicael/Documents/algo/%s/config2.py"%(name))
		if(os.path.isfile("fichierdeconfig2.cfg")):
			print ("the directory is ready for your output files\n")
		else:
			print("a problem occured and your directory is not ready\n")
	else:
		print("an error occured your directory couldn't be created\n")

			
def select2(question):
	exit = True
	while(exit):
		try:
			choice =  int(raw_input(question))#"would like to work in :\n 1: nucleic \n 2: nucleic/protein \n 3: protein \n (select by typing'1','2' or '3')  \n"
			if (choice!=1):
						if (choice!=2):
							raise NumberError
		except ValueError:
			print("you didn't give a number, try again \n ")
		except NumberError:
			print("try again")
		else:
			exit = False
			break
	return choice

def select(question):
	exit = True
	while(exit):
		try:
			choice =  int(raw_input(question))#"would like to work in :\n 1: nucleic \n 2: nucleic/protein \n 3: protein \n (select by typing'1','2' or '3')  \n"
			if (choice!=1):
						if (choice!=2):
							if (choice!=3):
								raise NumberError
		except ValueError:
			print("you didn't give a number, try again \n ")
		except NumberError:
			print("try again")
		else:
			exit = False
			break
	return choice
	
def select5(question):
	exit = True
	while(exit):
		try:
			choice =  int(raw_input(question))
			if (choice!=1):
				if (choice!=2):
					if (choice!=3):
						if (choice!=4):
							if (choice!=5):
								raise NumberError2
		except ValueError:
			print("you didn't give a number, try again \n ")
		except NumberError2:
			print("try again")
		else:
			exit = False
			break
	return choice
	
def bank_prot(filein,fileout,key):
	handle = open(fileout,'w')
	for record in SeqIO.parse(filein,'fasta'):
		if key in record.description:
			handle.write(">"+str(seq_record.description)+"\n")
			handle.write(str(seq_record.seq)+"\n")
		count = count + 1

@check_if_uptodate(check_file_exists)	
#@originate(job1_input1,job1_input3,job1_input2,job1_output2,job1_output1)	
@files(job1_input1,job1_input3,job1_input2,job1_output2,job1_output1)
#@follows(mkdir("pipeline_output"))
def bankChoice(file,file4,file3,file2,file1):
	'''cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	file = cfg.get('Filename', 'cle3')#ICTV file
	file1 = cfg.get('Filename', 'cle4')#nucleotide output file
	file2 = cfg.get('Filename', 'cle4(0)')#protein output file 
	file3 = cfg.get('Filename', 'cle5')#genbank input file carrying the nucleotide sequences (used to create the data bank)
	file4 = cfg.get('Filename', 'cle5(0)')#genbank input file carrying the protein sequences (used to create the data bank)'''
	bank_type = select("what kind of molecule would you like to work on? \n 1: nucleotide (press '1' then 'enter') \n 2: nucleotide and protein \n 3: protein (press '3' then 'enter')\n")
	for case in switch(bank_type):
		if case(1):
			print("the nucleotide sequences are about to get retrieve from the viral genbank file\n")
			data_bank2(file3, file1, file)
			break;
		if case(2):
			print("the nucleotide sequences are about to get retrieved from the "+file3+" file\n")
			data_bank2(file3, file1, file)
			print("now the protein sequences are also about to get retrieved from the "+file4+" file\n")
			data_bank2(file4, file2, file)
			break;
		if case(3):
			print("now the protein sequences are also about to get retrieved from the "+file4+" file\n")
			data_bank2(file4, file2, file)
			break;
	
def data_bank2(input, output, ictv):
	test1 = True
	count1 = 0
	count2 = 0
	choice = select("in order to create your data bank, you can use:\n 1: a key word (press '1' then 'enter') \n 2: the taxonomy in compliance with the classification from the International Committee on Viruses (press '2' then 'enter') \n 3: the accession or the gi number (press '3' then 'enter')  \n")
	choice = str(choice)
	choice = choice.replace(" ","")
	choice = int(choice)
	for case in switch(choice):
		if case(1):
			while(True):
				key_word = raw_input("give a key word in order to retrieve all of the sequences related: ").lower()
				#key_word = key_word.replace(" ","")
				if key_word:
					break
				else:
					print("you didn't type anything")
			output_handle = open(output,"w")
			print("your data bank is being created, it may take a few minutes \n")
			for seq_record in SeqIO.parse(input, "genbank"):
				if key_word in seq_record.annotations["organism"].lower() or key_word in seq_record.annotations["source"].lower():
					output_handle.write(">"+seq_record.annotations["organism"])
					output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
					count1 = count1 + 1
			if (count1==0):
				print("no organism related to "+key_word)
				raise TaskError
			else:
				print(str(count1)+" organisms have been added in "+output+" from "+input)
			output_handle.close()
			return 1
			break;
		if case(2):
			book = open_workbook(ictv)
			sheet = book.sheet_by_index(1)
			output_handle = open(output,"w")
			taxon = raw_input("according to the ICTV classification your research should use the Order (type 'order' then press enter) , the Family (type 'family' then press enter), the Subfamily (type 'subfamily' then press enter),the Genus (type 'genus' then press enter) or the Specie (type 'specie' then press enter):   ").lower()
			taxon = taxon.replace(" ","")
			test = True
			while(test):
				if(taxon=="order" or taxon=="family" or taxon=="subfamily" or taxon=="genus" or taxon=="specie"):
					break
				taxon = raw_input("eihter the taxon you have typed doesn't match the ICTV classification or you have typed wrong, you must choose by typing 'order', 'family', 'subfamily', 'genus' or 'specie':  ").lower()
				taxon = taxon.replace(" ","")
				if (taxon=="order" or taxon=="family" or taxon=="subfamily" or taxon=="genus" or taxon=="specie"):
					test = False
				else:
					test = True		
			key_word = raw_input("give the taxonomic name related to the taxon you have chosen: ").lower()
			key_word = key_word.replace(" ","")
			print("your data bank is being created, it may take a few minutes \n")
			for case in switch(taxon):
				if case('order'):
					for row_index in range(sheet.nrows):
						if key_word in sheet.cell(row_index,0).value.lower():
							count2 = count2 + 1
							for seq_record in SeqIO.parse(input, "genbank"):
								if key_word.lower() in map(lambda x:x.lower(),seq_record.annotations["taxonomy"]):
									output_handle.write(">"+seq_record.annotations["organism"])
									output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
									count1 = count1 + 1
							if count1==0:
								print ("no organism related to "+key_word+" has been found in "+input)
							else:
								print(str(count1)+" organisms have been added in "+output+" from "+input)
						if count2==1:
							break			
					if count2==0:
						print("the taxonomic name, you have given, doesn't belong to the taxon 'order' from ICTV viral classification")
						raise TaskError
					break;
				if case('family'):
					for row_index in range(sheet.nrows):
						if key_word in sheet.cell(row_index,1).value.lower():
							count2 = count2 + 1
							for seq_record in SeqIO.parse(input, "genbank"):
								if key_word in map(lambda x:x.lower(),seq_record.annotations["taxonomy"]):
									output_handle.write(">"+seq_record.annotations["organism"])
									output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
									count1 = count1 + 1
							if count1==0:
								print ("no organism related to "+key_word+" has been found in "+input)
							else:
								print(str(count1)+" organisms have been added in "+output+" from "+input)
						if count2==1:
							break				
					if count2==0:
						print("the taxonomic name, you have given, doesn't belong to the taxon 'family' from ICTV viral classification")
						raise TaskError
					break;
				if case('subfamily'):
					for row_index in range(sheet.nrows):
						if key_word in sheet.cell(row_index,2).value.lower():
							count2 = count2 + 1
							for seq_record in SeqIO.parse(input, "genbank"):
								if key_word in map(lambda x:x.lower(),seq_record.annotations["taxonomy"]):
									output_handle.write(">"+seq_record.annotations["organism"])
									output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
									count1 = count1 + 1
							if count1==0:
								print ("no organism related to "+key_word+" has been found in "+input)
							else:
								print(str(count1)+" organisms have been added in "+output+" from "+input)
						if count2==1:
							break					
					if count2==0:
						print("the taxonomic name, you have given, doesn't belong to the taxon 'subfamily' from ICTV viral classification")
						raise TaskError
					break;
				if case('genus'):
					for row_index in range(sheet.nrows):
						if key_word in sheet.cell(row_index,3).value.lower():
							count2 = count2 + 1
							for seq_record in SeqIO.parse(input, "genbank"):
								if key_word in map(lambda x:x.lower(),seq_record.annotations["taxonomy"]):
									output_handle.write(">"+seq_record.annotations["organism"])
									output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
									count1 = count1 + 1
							if count1==0:
								print ("no organism related to "+key_word+" has been found in "+input)
							else:
								print(str(count1)+" organisms have been added in "+output+" from "+input)
						if count2==1:
							break				
					if count2==0:
						print("the taxonomic name, you have given, doesn't belong to the taxon 'genus' from ICTV viral classification")
						raise TaskError
					break;
				if case('specie'):
					for row_index in range(sheet.nrows):
						if key_word in sheet.cell(row_index,5).value.lower():
							count2 = count2 + 1
							for seq_record in SeqIO.parse(input, "genbank"):
								if key_word in map(lambda x:x.lower(),seq_record.annotations["taxonomy"]):
									output_handle.write(">"+seq_record.annotations["organism"])
									output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
									count1 = count1 + 1
							if count1==0:
								print ("no organism related to "+key_word+" has been found in "+input)
							else:
								print(str(count1)+" organisms have been added in "+output+" from "+input)
						if count2==1:
							break				
					if count2==0:
						print("the taxonomic name, you have given, doesn't belong to the taxon 'specie' from ICTV viral classification")
						raise TaskError
					break;
				'''if case():
					print("eihter the taxon you tape doesn't match the ICTV classification or you have typed wrong, you must choose by typing 'order', 'family', 'subfamily', 'genus' or 'specie'")
					break;'''

			return 2
			break;
		if case(3):
			key_word = raw_input("give a accession or a gi number in order to retrieve all of the sequences related to them: ").upper()
			print("your data bank is being created, it may take a few minutes \n")		
			output_handle = open(output,"w")
			for seq_record in SeqIO.parse(file3, "genbank"):
				if key_word in seq_record.annotations["accessions"] or key_word==seq_record.annotations["gi"]:
					output_handle.write(">"+seq_record.annotations["organism"])
					output_handle.write("\n"+repr(seq_record.seq.tostring())+"\n")
					count1 = count1 + 1
			if (count1==0):
				print("no organism found under the accession/gi number "+key_word)
				raise TaskError
			output.close()
			return 3
			break;


	
	
def exist(count,file_out):
	size = os.path.getsize(file_out)
	if(isinstance(size, (int, long, float, complex))):	
		if (size > 0):
			print("Saved "+str(count)+" reads in the file "+file_out)
		else:
			print("an error occured, the '"+file_out+"' file is empty")
			raise TaskError
	else:
		print("an error occured, it seems the '"+file_out+"' file was not created")
		raise TaskError



	
def trim_reads(file,file_out):
	choice = "yes"
	while(choice == "yes"):
		#pexpect.run("fastqc", timeout=3600)
		os.system("fastqc")
		'''cfg = ConfigParser.ConfigParser()
		cfg.read('fichierdeconfig.cfg')
		file = cfg.get('Filename', 'cle2(0)')
		file_out = cfg.get('Filename', 'cle2(1)')'''
		exit = True
		while(exit):
			try:
				trim_size = int(raw_input('how many nucleotide(s) you want to strip off the end of the reads?'))
				if(trim_size < 0):
					raise ValueNegative
			except ValueError:
				print("you didn't give an integer, try again \n ")
			except ValueNegative:
				print("give an other number\n")
			else:
				exit = False
				break
		if(trim_size == 0):
			shutil.copy2(file,file_out)
			#trimmed_primer_reads = (rec for rec in \
			#							SeqIO.parse(file, "fastq"))
		else:
			trimmed_primer_reads = (rec[:-trim_size] for rec in \
										SeqIO.parse(file, "fastq"))
			count = SeqIO.write(trimmed_primer_reads, file_out, "fastq")
			size = os.path.getsize(file_out)
			if(isinstance(size, (int, long, float, complex))):	
				if (size > 0):
					print("Saved "+str(count)+" reads in the file "+file_out)
				else:
					print("an error occured, the '"+file_out+"' file is empty")
			else:
				print("an error occured, it seems the '"+file_out+"' file was not created")
				return False
				break
		#pexpect.run("fastqc", timeout=3600)
		os.system("fastqc")
		while(True):
			choice2 = raw_input("Are you satisfied with your read quality? enter 'yes' if you have finished and 'no' otherwise: ")
			if(choice2 == "yes"):
				print("now you are done with the trimming step")
				return True
				loop = False
			elif(choice2 !="no"):
				print("you didn't answer by 'yes' or 'no' \n \n try again ")
			elif(choice2 == "no"):
				choice = "yes"
				print("Since you are not satisfied  with the quality you have achieved, this step is about to start from the begining, using fastqc you will reassess the reads quality from "+file+" file and remove as many nucleotides as necessary\n")
				break
	
	

'''def trim_primers(records, primer):
	len_primer = len(primer)
	for record in records:
		if record.seq.startswith(primer):
			yield record[len_primer:]
		else:
			yield record'''

def trim_primers(records, primer):
    len_primer = len(primer) #cache this for later
    for record in records:
        if record.seq.startswith(primer):
            yield record[len_primer:]
        else:
            yield record

def primer_off():
	cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	file = cfg.get('Filename', 'cle2(0)')
	file_out = cfg.get('Filename', 'cle2(1)')
	primer = cfg.get('data_input', 'cle2')
	original_reads = SeqIO.parse(file, "fastq")
	trimmed_reads = trim_primers(original_reads, primer)
	count = SeqIO.write(trimmed_reads, file_out, "fastq")
	exist(count,file_out)
	
	

'''def phred(file,file_out):
	cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	file = cfg.get('Filename', 'cle2')
	file_out= cfg.get('Filename', 'cle2(0)')
	records = (rec for rec in \
		SeqIO.parse(file, "fastq") \
		if min(rec.letter_annotations[ "phred_quality"]) >= 20)
	count = SeqIO.write(records, file_out, "fastq")
	print( "%i reads with a phred score above 20" % count)'''
	
def phred2(file,file_out):
	records1 = (rec for rec in \
			SeqIO.parse(file, "fastq") \
			if min(rec.letter_annotations[ "phred_quality"]) >= 20)
	records2 = (rec for rec in \
			SeqIO.parse(file, "fastq") \
			if min(rec.letter_annotations[ "phred_quality"]) < 20)
	count1 = SeqIO.write(records1, file_out, "fastq")
	print( "%i reads with a phred score above 20" % count1)
	count2 = SeqIO.write(records2, 'test_phred_badq.fastq', "fastq")
	print( "%i reads with a phred score below 20" % count2)
	return count1,count2


@files(job2_input1,job2_output1,job2_output2,job2_output3)
@follows(bankChoice)
def fastqc(phred_in,phred_out,trim_out,sizefile):
	mean=size(phred_in,sizefile)
	print("the average size of the reads is: "+str(mean)+"\n")
	print("reads with a low phred quality are being removed...")
	count1,count2=phred2(phred_in,phred_out)
	tot=count1+count2
	if(count1 < tot*25/100):
		choice=select2("the read process led to a loss of data: more than 25% of the reads have been removed regarding their phred score.\n You have to choose if you want to go on further analysis from:\n 1\ the original file \n 2\ the one without low quality reads\n")
		if (choice==1):		
			print("We need to see the data quality, FASTQC is about to be opened, you will have to click on 'open' and browse to give it your file containing the raw data.\n After seeing the reads quality, close FASTQC so you can trim the reads in order to improve the reads quality")
			shutil.copy2(phred_in,phred_out)
			trim_reads(phred_out,trim_out)
		else:
			print("We need to see the data quality, FASTQC is about to be opened, you will have to click on 'open' and browse to give it your file containing the raw data.\n After seeing the reads quality, close FASTQC so you can trim the reads in order to improve the reads quality")
			trim_reads(phred_out,trim_out)
	else:
		print("We need to see the data quality, FASTQC is about to be opened, you will have to click on 'open' and browse to give it your file containing the raw data.\n After seeing the reads quality, close FASTQC so you can trim the reads in order to improve the reads quality")
		trim_reads(phred_out,trim_out)
	
	
	



def fastqToFasta(file_name, file_out):
	count = 0
	handle = open(file_out,"w")
	for seq_record in SeqIO.parse(file_name, "fastq"):
		 handle.write(">"+str(seq_record.description)+"\n")
		 handle.write(str(seq_record.seq)+"\n")
		 count = count + 1
	handle.close()
	return count







@files(job3_input1,job3_input2,job3_output1,job3_output2)
#@follows(fastqc)
def bowtie2(rawfile,rawfile2,mapfile,unmapfile):
	cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	#mapfile = cfg.get('Filename', 'cle9')
	#unmapfile = cfg.get('Filename', 'cle8')
	#rawfile = cfg.get('Filename', 'cle2(1)')#bowtie input file 
	index = cfg.get('Filename', 'cle10')
	#match = cfg.get('Filename', 'cle12')
	#unmatch = cfg.get('Filename', 'cle11')
	print("for avoiding all the problems related to fastq format the "+rawfile+" file is being converted to fasta format before the mapping")
	fastqToFasta(rawfile,rawfile2)
	print("the mapping is in process")
	#pexpect.run("bowtie2 -f --very-sensitive-local --un %s --al %s -x  %s -U %s"%(unmapfile, mapfile, index, rawfile2))#figure the option --no-unal out
	os.system("bowtie2 -f --very-sensitive-local --un %s --al %s -x  %s -U %s"%(unmapfile, mapfile, index, rawfile2))#figure the option --no-unal out
	#count1 = fastqToFasta(unmapfile, unmatch)
	count1 = pexpect.run('count.sh %s'%(unmapfile))
	print(str(count1)+" unmapped reads were saved in "+unmapfile+" \n")
	#count2 = fastqToFasta(mapfile, match)
	count2 = pexpect.run('count.sh %s'%(mapfile))
	print(str(count2)+" mapped reads were saved in "+mapfile+" \n")
	'''choice = False
	handle = open(mapfile,"r")
	file = handle.read()
	handle2 = open(unmapfile,"r")
	file2 = handle2.read()
	loop = True
	while (loop):
		choice = raw_input("would you like to display the reads that mapped? (type 'yes' or 'no'):  ")
		if( choice == "yes" ):
			print(file)
			break
		elif ( choice =="no" ):
			break
		else:
			print("sorry your answer was not clear, you have to type 'yes' or 'no'")
		if( loop=='yes' or loop=='no'):
			loop = False
		else:
			loop = True
	loop2 = True
	while (loop2):
		choice = raw_input("would you like to display the reads that didn't mapped? (type 'yes' or 'no'):  ")
		if( choice == "yes" ):
			print(file2)
			break
		elif ( choice =="no" ):
			break
		else:
			print("sorry your answer was not clear, you have to type 'yes' or 'no'")
		if( loop=='yes' or loop=='no'):
			loop = False
		else:
			loop = True
	handle.close()
	handle2.close()'''
	
	
	
			
			
			
			
def awk_blast_out(before_blast,blastout,outawk):
	exist = os.path.isfile(blastout)
	size = os.path.getsize(blastout)
	if(exist):
		if(size>0):
			print("the blast algorithm has stop running without errors")
			print("the exact positions (at the beginning and at the end of the reads) that have been aligned are being calculated\n")
			#p = pexpect.run("awk -f select_target_blast.awk %s %s"%(before_blast, blastout))
			#print p
			os.system("awk -f select_target_blast.awk %s %s > %s"%(before_blast, blastout,outawk))
			#blast_position = open(outawk, "w")
			#blast_position.write(p)
			print("you can find the entire sequence of the aligned reads from the first aligned nucleotide to the last one in "+outawk+" file\n")
		else:
			print("No reads from the input file "+before_blast+" have been aligned. The output file "+blastout+" is empty\n")
			#raise TaskError
	else:
		print("the blast query didn't produce any result file")
		raise TaskError


def makedbprot(input,output):
	print("the files required for the blastx database are being created")
	#o = pexpect.run("makeblastdb -in %s -dbtype 'prot' -out %s"%(input, output))
	#print o
	os.system("makeblastdb -in %s -dbtype 'prot' -out %s"%(input, output))
	size1 = os.path.getsize("%s.phr"%(output))
	size2 = os.path.getsize("%s.pin"%(output))
	size3 = os.path.getsize("%s.psq"%(output))
	list = [size1,size2,size3]
	return list

def makedbnuc(input, output):
	print("the files required for the blastn database are being created")
	#o = pexpect.run("makeblastdb -in %s -dbtype 'nucl' -out %s"%(input, output))
	#print o
	os.system("makeblastdb -in %s -dbtype 'nucl' -out %s"%(input, output))
	size1 = os.path.getsize("%s.nhr"%(output))
	size2 = os.path.getsize("%s.nin"%(output))
	size3 = os.path.getsize("%s.nsq"%(output))
	list = [size1,size2,size3]
	return list
	
	
	
def dbnuc(input, output):
	if  os.path.exists("%s.nhr"%(output)):
		if os.path.exists("%s.nin"%(output)):
			if os.path.exists("%s.nsq"%(output)):
				print("the three files '%s.nhr', '%s.nin' and '%s.nsq' required for the blastn database already exists\n"%(output, output, output))
				return 1
			else:
				print("the file %s.nsq is missing"%(output))
				list = makedbnuc(input, output)
				return list
		else:
			print("the files '%s.nin' is missing\n"%(output))
			list = makedbnuc(input, output)
			return list
	else:
		print("the files '%s.nhr' is missing\n"%(output))
		list = makedbnuc(input, output)
		return list
		
def dbprot(input, output):
	if  os.path.exists("%s.phr"%(output)):
		if os.path.exists("%s.pin"%(output)):
			if os.path.exists("%s.psq"%(output)):
				print("the three files '%s.phr', '%s.pin' and '%s.psq' required for the blastn database already exists\n"%(output, output, output))
				return 1
			else:
				print("the file %s.psq is missing"%(output))
				list = makedbnuc(input, output)
				return list
		else:
			print("the files '%s.pin' is missing\n"%(output))
			list = makedbnuc(input, output)
			return list
	else:
		print("the files '%s.phr' is missing\n"%(output))
		list = makedbnuc(input, output)
		return list	
			
			
			
#select blast output reads with a length above 20 (the kmer size)
def select_blastout(filein,fileout):
	#handle = open(fileout,"w")
	#p = pexpect.run("awk '//{if ($3>20) {print $0}}' %s"%(filein))
	#handle.write(p)
	#handle.close()
	os.system("awk '//{if ($3>20) {print $0}}' %s > %s"%(filein,fileout))
	
def fusion(filein1,filein2,fileout):
	filenames = [filein1,filein2]
	with open(fileout, 'w') as outfile:
		for fname in filenames:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)
			
def removeQuote(filein,fileout):
	handle = open(fileout,"w")		
	for seq_record in SeqIO.parse(filein, "fasta"):
		handle.write(">"+str(seq_record.description)+"\n")
		handle.write(str(seq_record.seq).replace("'","")+"\n")
	handle.close()
	
@files(job4_input1,job4_input2,job4_input3,job4_input4,job4_output1,job4_output2,job4_key1,job4_output3,job4_key2,job4_output4,job4_input5,job4_output5,job4_output6,job4_output7,job4_output8,job4_output9,job4_output10)
@follows(bowtie2)
def blast(file1,file11,file120,file121,file2,file22,evalue,output,format,outawk,unmatch,neg,blastx_output,outawk_blastx,pre_awk_blastn,pre_awk_blastx,final):
	'''cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	file1 = cfg.get('Filename', 'cle4')#nucleotide input file, before quotes removal
	file11 = cfg.get('Filename', 'cle4(1)')#nucleotide input file after quotes removal
	file120 = cfg.get('Filename', 'cle4(0)')#protein input file, before quotes removal
	file121 = cfg.get('Filename', 'cle4(3)')#protein input file, after quotes removal
	file2 = cfg.get('Filename', 'cle13')#nucleotide database created (dbnuc output)
	file22 = cfg.get('Filename', 'cle33')#protein database created (dbprot output)
	evalue = cfg.get('data_input', 'cle3')#evalue of the blast
	output = cfg.get('Filename', 'cle15')#nuclotide blast output
	format = cfg.get('data_input', 'cle4')#blastn format output
	outawk = cfg.get('Filename', 'cle17')#awk function output that give the aligned reads with exact posi
	unmatch = cfg.get('Filename', 'cle11(1)')#blastn input
	neg = cfg.get('Filename', 'cle30')#awk function output that gives the unaligned reads
	blastx_output = cfg.get('Filename', 'cle31')#blastx output
	outawk_blastx = cfg.get('Filename', 'cle32')#awk function output after blastx 
	pre_awk_blastn = cfg.get('Filename', 'cle34')#blastn output used to select reads regarding their size (>50)
	pre_awk_blastx =cfg.get('Filename', 'cle35')#blastx output used to select reads regarding their size (>50)'''
	#choice of the type of molecule you'd like to work on (nucleic or protein)
	choice = select("would like to work in :\n 1: nucleic \n 2: nucleic/protein \n 3: protein \n (select by typing'1','2' or '3')  \n")
	for case in switch(choice):	
			if case(1):
				removeQuote(file1,file11)
				size = dbnuc(file11,file2)
				exist = True
				if not (size==1):
					for item in size:
						if(item==0):
							exist = False
				if(exist):
					print("the blastn is in process \n")
					#pexpect.run("blastn -db %s -query %s -outfmt %s -out %s -evalue %s"%(file2, unmatch, format, output, evalue))
					os.system("blastn -max_target_seqs 1 -db %s -query %s -outfmt %s -out %s -evalue %s"%(file2, unmatch, format, output, evalue))
					#select_blastout(output,pre_awk_blastn)
					#awk_blast_out(unmatch,pre_awk_blastn,outawk)
					awk_blast_out(unmatch,output,outawk)
					shutil.copy2(outawk,final)
				else:
					print("sorry an error occured, the blastn database can't be used")
					raise TaskError
				break;
			if case(2):
				removeQuote(file1,file11)
				removeQuote(file120,file121)
				size = dbnuc(file11,file2)
				exist = True
				if not (size == 1):
					for item in size:
						if(item==0):
							exist = False
				if(exist):
					print("the blastn is in process \n")
					#pexpect.run("blastn -db %s -query %s -outfmt %s -out %s -evalue %s"%(file2, unmatch, format, output, evalue))
					os.system("blastn -db %s -query %s -outfmt %s -out %s -evalue %s"%(file2, unmatch, format, output, evalue))
					select_blastout(output,pre_awk_blastn)
					awk_blast_out(unmatch,pre_awk_blastn,outawk)
					#p = pexpect.run("awk -f selection_read_neg.awk %s %s"%(unmatch,output))
					#handle = open(neg, "w")
					#handle.write(p)
					os.system("awk -f selection_read_neg.awk %s %s > %s"%(unmatch,output,neg))
					#remove the first 2 lines
					lines = open(neg).readlines()
					open(neg, 'w').writelines(lines[2:])
					size = dbprot(file121, file22)
					exist = True
					if not (size == 1):
						for item in size:
							if(item==0):
								exist = False
					if(exist):
						print("the blastx is in process \n")
						#p = pexpect.run("blastx -db %s -query %s -outfmt %s -out %s -evalue %s"%(file22, neg, format, blastx_output, evalue))
						#print p
						os.system("blastx -db %s -query %s -outfmt %s -out %s -evalue %s"%(file22, neg, format, blastx_output, evalue))
						select_blastout(blastx_output,pre_awk_blastx)
						awk_blast_out(neg,pre_awk_blastx,outawk_blastx)
						#il faut concatener les deux fichiers a la fin (outawk et outawk_blastn)
						fusion(outawk,outawk_blastn,final)
					else:
						print("sorry an error occured, the blastx database can't be used")
						raise TaskError
				else:
					print("sorry an error occured, the blastn database can't be used")
					raise TaskError
				break;
			if case(3):
				removeQuote(file120,file121)
				size = dbprot(file121, file22)
				exist = True
				if not (size == 1):
					for item in size:
						if(item==0):
							exist = False
				if(exist):
					print("the blastx is in process \n")
					#p = pexpect.run("blastx -db %s -query %s -outfmt %s -out %s -evalue %s"%(file22, unmatch, format, blastx_output, evalue))
					os.system("blastx -db %s -query %s -outfmt %s -out %s -evalue %s"%(file22, unmatch, format, blastx_output, evalue))
					select_blastout(blastx_output,pre_awk_blastx)
					awk_blast_out(unmatch,output,outawk_blastx)
					shutil.copy2(outawk_blastx,final)
				else:
					print("sorry an error occured, the blastx database can't be used")
					raise TaskError
				break;

'''#@follows(bowtie2, mkdir("output/results/here"))
def blast():
	cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	file1 = cfg.get('Filename', 'cle4')
	file11 = cfg.get('Filename', 'cle4(1)')
	file2 = cfg.get('Filename', 'cle13')
	evalue = cfg.get('data_input', 'cle3')
	output = cfg.get('Filename', 'cle15')
	format = cfg.get('data_input', 'cle4')
	outawk = cfg.get('Filename', 'cle17')
	unmatch = cfg.get('Filename', 'cle11(1)')
	handle = open(file11,"w")
	#choice of the type of
	choice = select("would like to work in :\n 1: nucleic \n 2: nucleic/protein \n 3: protein \n (select by typing'1','2' or '3')  \n")
	#removal of the quote from the file generated by the data_bank2() function
	for seq_record in SeqIO.parse(file1, "fasta"):
		handle.write(">"+str(seq_record.description)+"\n")
		handle.write(str(seq_record.seq).replace("'","")+"\n\n")
	handle.close()    
	#creation of the database you can blast against 
	print("your data bank is being formated so it can be used as a blast database in a query \n ")
	o = pexpect.run("makeblastdb -in %s -dbtype 'nucl' -out %s"%(file11, file2))
	print o
	size = os.path.getsize("%s.nhr"%(file2))
	size2 = os.path.getsize("%s.nin"%(file2))
	size3 = os.path.getsize("%s.nsq"%(file2))
	count = 0
	if(isinstance(size, (int, long, float, complex))):
		print("your database is ready to be used in your blast query  \n")
		print("the blast is in process \n")
		pexpect.run("blastn -db %s -query %s -outfmt %s -out %s -evalue %s"%(file2, unmatch, format, output, evalue))
		exist = os.path.isfile(output)
		size = os.path.getsize(output)
		if(exist):
			if(size>0):
				print("the blast algorithm has stop running without errors")
				print("the exact positions (at the beginning and at the end of the reads) that have been aligned are being calculated\n")
				p = pexpect.run("awk -f select_target_blast.awk %s %s"%(unmatch, output))
				print p
				blast_position = open(outawk, "w")
				blast_position.write(p)
				print("you can find the entire sequence of the aligned reads from the first aligned nucleotide to the last one in "+outawk+" file\n")
			else:
				print("nothing have been aligned since the output file "+output+" is empty\n")
				raise TaskError		
		else:
			print("the blast query didn't produce any result file")
			raise TaskError
	else:
		print("sorry an error occured, the database can't be used")'''
		
		
def kmer(inputfile):
	while(exit):
		try:
			size = int(raw_input("in order to reduce the gap between the reads size, you can give a minimal size for the aligned reads :  \n"))
			if(size>50):
				raise Kmer
			if(size<20):
				raise Kmer
		except ValueError:
			print("you didn't give a number, try again \n ")
		except Kmer:
			print("the size is between 20 and 50, try again")
		else:
			exit = False
			break
	

@files(job4_output10,job5_output_dir,job5_output1)
#@follows(blast)
def spades(filein,fileout,contigs):
	'''cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	filein = cfg.get('Filename', 'cle17')#spades input
	fileout = cfg.get('Filename', 'cle21')#spades output'''
	print("the spades reads assembly is starting...\n")
	#p = pexpect.run("spades.py --only-assembler -o %s -s %s"%(fileout, filein))
	#print p
	os.system("spades.py --only-assembler -o %s -s %s"%(fileout, filein))
	count = len(os.listdir('spades_out'))
	if(count>2):
		print("the current task is done")
	else:
		print("Oups, an error seems to have occured")
		raise TaskError
	try:
		contig = os.path.join(fileout,contigs)
		check_contig=os.path.getsize(contig)
		if(check_contig == 0):
			contigs_path = os.path.join(fileout, 'K21/final_contigs.fasta')
			open(fileout+'/contigs.fasta','w')
			shutil.copy2(contigs_path,fileout+'/contigs.fasta')
	except OSError:
		print("there is no 'contigs.fasta' file in "+fileout+" directory,therefore we have to use the one in K21 subdirectory")
		contigs_path = os.path.join(fileout, 'K21/final_contigs.fasta')
		open(fileout+'/contigs.fasta','w')
		shutil.copy2(contigs_path,fileout+'/contigs.fasta')
			
			
#@follows(spades, "assembly step", mkdir("output/results/here"))
def contig_plot():
	cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	file0 = cfg.get('Filename', 'cle21')
	file1 = cfg.get('Filename', 'cle4(2)')
	file = 	"/"+file0+"/"+file1
	sizes = [len(rec) for rec in SeqIO.parse(file, "fasta")]
	print("there are "+str(len(sizes))+" contigs, with a size varying from "+str(min(sizes))+" to "+str(max(sizes))+" \n the distribution of the sequence lengths is about to pop up")
	pylab.hist(sizes, bins=20)
	pylab.title("%i contigs\nLengths %i to %i" \
				% (len(sizes),min(sizes),max(sizes)))
	pylab.xlabel("Sequence length (bp)")
	pylab.ylabel("Count")
	pylab.show()

@files(job4_output10,job6_output_dir,job6_key1,job6_output1)
@follows(spades)
def ray(filein,fileout,kmer,contig_file):
	'''cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	filein = cfg.get('Filename', 'cle20')
	fileout = cfg.get('Filename', 'cle22')
	kmer = cfg.get('data_input','cle5')'''
	print("the ray reads assembly is starting...\n")
	#p = pexpect.run("mpiexec Ray -k %s -o %s -s %s"%(kmer, fileout, filein))
	#print p
	os.system("mpiexec Ray -k %s -o %s -s %s"%(kmer, fileout, filein))
	try:
		contigs_path = os.path.join(fileout, contig_file)
		exist = os.path.getsize(contigs_path)
		if (exist == 0):
			print(" Error : there no contigs yielded by ray!!!")
			raise TaskError
	except OSError:
		print("there is no 'contigs.fasta' file in "+fileout+" directory,therefore we have to use the one in K55 subdirectory")



	

'''def singleton():
	cfg = ConfigParser.ConfigParser()
	cfg.read('fichierdeconfig.cfg')
	output_blast = cfg.get('Filename', 'cle17')
	map = cfg.get('Filename', 'cle26')
	unmap = cfg.get('Filename', 'cle27')
	spades_contig = cfg.get('Filename', 'cle28')
	index = cfg.get('Filename', 'cle29')
	file0 = cfg.get('Filename', 'cle21')
	file1 = cfg.get('Filename', 'cle4(2)')
	output_assembly = os.path.join(file0, file1)
	shutil.copy2(output_assembly,spades_contig)
	size = os.path.getsize(spades_contig)
	bam = os.path.splitext(spades_contig)[0]+".bam"
	sorted = os.path.splitext(spades_contig)[0]+".sorted"
	mappedbam = os.path.splitext(spades_contig)[0]+".mapped.bam"
	unmappedbam = os.path.splitext(spades_contig)[0]+".unmapped.bam"
	#sai = os.path.splitext(spades_contig)[0]+".sai"
	sam = os.path.splitext(spades_contig)[0]+".sam"
	choice = select2("BWA OR BOWTIE2? :\n 1: BWA \n 2: BOWTIE2 \n")
	for case in switch(choice):	
		if case(1):
			#bwa mem
			pexpect.run("bwa index %s"%(output_blast))
			handle2 = open(sam, 'w')
			p = pexpect.run("bwa mem %s %s"%(output_blast, spades_contig))
			handle2.write(p)
			handle2.close()
			pexpect.run("samtools view -h -b -S -o %s %s"%(bam,sam))
			pexpect.run("samtools sort %s %s"%(bam, sorted))
			sortedbam = sorted+".bam"
			print sortedbam
			p = pexpect.run("samtools view -c  -F 4 %s"%(sortedbam))
			print("unmapped read number: "+p+"\n")
			p = pexpect.run("samtools view -c -f 4 %s"%(sortedbam))
			print p
			p = pexpect.run("samtools view -c -b -F 4 -o %s %s"%(mappedbam, bam))
			print p
			p = pexpect.run("samtools view -b -f4 -o %s %s"%(unmappedbam, bam))
			print p
			break;
		if case(2):
			#avec bowtie
			if(size > 0):
				pexpect.run("bowtie2-build %s %s"%(spades_contig, index))
				p = pexpect.run("bowtie2 -f --un unmap.fasta --al map.fasta -x  %s -U %s"%(index, output_blast))
				print p
			else:
				print("sorry an error occured during the "+spades_contig+" copying")'''
				
@follows(ray)				
@files(job5_output_dir,job5_output1,job7_output1,job6_output_dir,job6_output1,job7_output2,job4_output10,job7_output3,job7_output4)				
def singleton(spades_out,contig_spades,singleton_out_spades,ray_dir,ray_contig,singleton_out_ray,reference,entire_viral_spades,entire_viral_ray):
	file_in_spades = os.path.join(spades_out,contig_spades)
	file_ray = os.path.join(ray_dir,ray_contig)
	print("the index from the assembled contigs is being built\n")
	o = os.system("bowtie2-build %s 'assembled_index'"%(reference))
	if(o==0):
		print("the mapping with spades output is in process")
		#pexpect.run("bowtie2 -f --very-sensitive-local --un %s -x 'assembled_index' -U %s"%(singleton_out_spades,file_in_spades))
		os.system("bowtie2 -f --very-sensitive-local --un %s -x 'assembled_index' -U %s"%(singleton_out_spades,file_in_spades))
		print("the mapping with ray output is in process")
		#pexpect.run("bowtie2 -f --very-sensitive-local --un %s -x 'assembled_index' -U %s"%(singleton_out_ray,file_ray))
		os.system("bowtie2 -f --very-sensitive-local --un %s -x 'assembled_index' -U %s"%(singleton_out_ray,file_ray))
		fusion(file_in_spades,singleton_out_spades,entire_viral_spades)
		fusion(file_ray,singleton_out_ray,entire_viral_ray)
	else:
		print("an error occured, the index couldn't be created")
		raise TaskError
		
		
		
#@follows(singleton)		
@files(job7_output3,job7_output4,job4_output1,job7_key1,job8_output1,job8_output2,job8_output3,job8_output4)
#follows(contig_plot, mkdir("output/results/here"))
def taxo(total_spades,total_ray,blast_db,blast_format,blast_spades_output,blast_ray_output,hpv,hpv2):
	#cfg = ConfigParser.ConfigParser()
	#cfg.read('fichierdeconfig.cfg')
	#file2 = cfg.get('Filename', 'cle13')
	#evalue = cfg.get('data_input', 'cle3')
	#output = cfg.get('Filename', 'cle23')
	#format = cfg.get('Filename', 'cle6')
	#file0 = cfg.get('Filename', 'cle21')
	#file1 = cfg.get('Filename', 'cle4(2)')
	#file = 	"/"+file0+"/"+file1
	print("the taxon of the virus is being determined with blast: \n")
	print("-at first with the spades contigs and the potential singletons\n")
	#pexpect.run("blastn -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, total_spades, blast_format, blast_spades_output))	
	os.system("blastn -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, total_spades, blast_format, blast_spades_output))
	xml(blast_spades_output,hpv)
	print("\n -secondly with the ray contigs and the potenetial singletons\n")
	#pexpect.run("blastn -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, total_ray, blast_format, blast_ray_output))
	os.system("blastn -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, total_ray, blast_format, blast_ray_output))
	xml(blast_ray_output,hpv2)
	choice = False
	handle = open(blast_spades_output,"r")
	file = handle.read()
	handle2 = open(blast_ray_output,"r")
	file2 = handle2.read()
	loop = True
	while (loop):
		choice = raw_input("[spades assembly] Would you like to see the file containing the viral types present in the samples (type 'yes' or 'no'):  ")
		if( choice == "yes" ):
			print(file)
			break
		elif ( choice =="no" ):
			break
		else:
			print("sorry your answer was not clear, you have to type 'yes' or 'no'")
		if( loop=='yes' or loop=='no'):
			loop = False
		else:
			loop = True
	loop2 = True
	while (loop2):
		choice = raw_input("[ray assembly] Would you like to see the file containing the viral types present in the samples (type 'yes' or 'no'):  ")
		if( choice == "yes" ):
			print(file2)
			break
		elif ( choice =="no" ):
			break
		else:
			print("sorry your answer was not clear, you have to type 'yes' or 'no'")
		if( loop=='yes' or loop=='no'):
			loop = False
		else:
			loop = True
	handle.close()
	handle2.close()
	
	
#cut-off size is chosen for the taxonomic designation
def relevant_contig(filein,size_chosen,blast_db,blast_format,blast_contig_size_chosen,taxo_contig_size_chosen):
	choice = select5("What size of contigs is relevant (the cut off value) ?\n 1: <100 (press '1' then 'enter') \n 2: >100 (press '2' then 'enter') \n 3: >250 (press '3' then 'enter') \n 4: >500 (press '4' then 'enter') \n 5: >1000 (press '5' then 'enter')  \n")
	choice = str(choice)
	choice = choice.replace(" ","")
	choice = int(choice)
	count = 0
	for case in switch(choice):
		if case(1):
			handle=open(size_chosen,'w')
			for record in SeqIO.parse(open(filein, "r"), "fasta"):
				if(len(record.seq) < 100):
					handle.write(">"+str(record.description)+"\n")
		 			handle.write(str(record.seq)+"\n")
		 			count = count + 1
			handle.close()
			print ("there is "+count+" contigs with a size lower than 100")
			os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_choosen, blast_format, blast_contig_size_chosen))
			xml2(blast_contig_size_chosen,taxo_contig_size_chosen)
			break;
		if case(2):
			handle=open(size_chosen,'w')
			for record in SeqIO.parse(open(filein, "r"), "fasta"):
				if(len(record.seq) >= 100):
					handle.write(">"+str(record.description)+"\n")
		 			handle.write(str(record.seq)+"\n")
		 			count = count + 1
			handle.close()
			print ("there is "+str(count)+" contigs with a size longer than 100")
			os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_chosen, blast_format, blast_contig_size_chosen))
			xml2(blast_contig_size_chosen,taxo_contig_size_chosen)
			break;
		if case(3):
			handle=open(size_chosen,'w')
			for record in SeqIO.parse(open(filein, "r"), "fasta"):
				if(len(record.seq) >= 250):
					handle.write(">"+str(record.description)+"\n")
		 			handle.write(str(record.seq)+"\n")
		 			count = count + 1
			handle.close()
			print ("there is "+str(count)+" contigs with a size longer than 250")
			os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_chosen, blast_format, blast_contig_size_chosen))
			xml2(blast_contig_size_chosen,taxo_contig_size_chosen)
			break;
		if case(4):
			handle=open(size_chosen,'w')
			for record in SeqIO.parse(open(filein, "r"), "fasta"):
				if(len(record.seq) >= 500):
					handle.write(">"+str(record.description)+"\n")
		 			handle.write(str(record.seq)+"\n")
	 				count = count + 1
			handle.close()
			print ("there is "+str(count)+" contigs with a size longer than 500")
			os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_chosen, blast_format, blast_contig_size_chosen))
			xml2(blast_contig_size_chosen,taxo_contig_size_chosen)
			break;
		if case(5):
			handle=open(size_chosen,'w')
			for record in SeqIO.parse(open(filein, "r"), "fasta"):
				if(len(record.seq) >= 1000):
					handle.write(">"+str(record.description)+"\n")
		 			handle.write(str(record.seq)+"\n")
		 			count = count + 1
			handle.close()
			print ("there is "+str(count)+" contigs with a size longer than 1000")
			os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_chosen, blast_format, blast_contig_size_chosen))
			xml2(blast_contig_size_chosen,taxo_contig_size_chosen)
			break;




	

@files(job5_output_dir,job5_output1,job6_output_dir,job6_output1,job9_output1,job9_output2,job4_output1,job7_key1,job9_output3,job9_output4,job9_output5,job9_output6,job9_output7,job9_output8,job9_output9,job9_output10,job4_output10)

#@follows(taxo)
def plot_contig(spades_out,contig_spades,ray_dir,ray_contig,size_chosen_spades,size_chosen_ray,blast_db,blast_format,blast_contig_size_chosen_spades,blast_contig_size_chosen_ray,taxo_contig_size_chosen_spades,taxo_contig_size_chosen_ray,map,unmap,map2,unmap2,blastout):
	file_in_spades = os.path.join(spades_out,contig_spades)
	file_ray = os.path.join(ray_dir,ray_contig)
	cS1 = cS2 = cS3 = cS4 = cS5 = cS6 = cS7 = cS8 = cS9 = cS10 = 0
	for record in SeqIO.parse(open(file_in_spades, "r"), "fasta"):
		if(len(record.seq) < 100):
			cS1 = cS1 + 1
		else:
			if(len(record.seq) < 250):
				cS2 = cS2 + 1
			else:
				if(len(record.seq) < 500):
					cS3 = cS3 + 1
				else:
					if(len(record.seq) < 1000):
						cS4 = cS4 + 1
					else:
						cS5 = cS5 + 1
	#print("cs1: "+cS1+", cs2: "+cS2+", cs3:"+cS3+", cs4:"+cS4+", cs5:"+cS5)
	
	for record in SeqIO.parse(open(file_ray, "r"), "fasta"):
		if(len(record.seq) < 100):
			cS6 = cS6 + 1
		else:
			if(len(record.seq) < 250):
				cS7 = cS7 + 1
			else:
				if(len(record.seq) < 500):
					cS8 = cS8 + 1
				else:
					if(len(record.seq) < 1000):
						cS9 = cS9 + 1
					else:
						cS10 = cS10 + 1
		
	plot_histo(cS1,cS2,cS3,cS4,cS5,cS6,cS7,cS8,cS9,cS10)
	file_in_spades = os.path.join(spades_out,contig_spades)
	file_ray = os.path.join(ray_dir,ray_contig)
	print("The chosen contigs of spades are being identified")
	relevant_contig(file_in_spades,size_chosen_spades,blast_db,blast_format,blast_contig_size_chosen_spades,taxo_contig_size_chosen_spades)
	print("The chosen contigs of ray are being identified")
	relevant_contig(file_ray,size_chosen_ray,blast_db,blast_format,blast_contig_size_chosen_ray,taxo_contig_size_chosen_ray)
	
	print("The viral reads are about to be mapped against the selected contigs from spades\n")
	print("The index from the selected contigs is being built\n")
	o1 = os.system("bowtie2-build %s 'spades_selected_contigs_index'"%(size_chosen_spades))
	if(o1==0):
		print("the mapping is starting")
		os.system("bowtie2 -f --very-sensitive-local --al %s --un %s -x 'spades_selected_contigs_index' -U %s"%(map,unmap,blastout))
	else:
		print("an error occured, the index couldn't be created")
		raise TaskError
	print("The viral reads are about to be mapped against the selected contigs from ray\n")
	print("The index from the selected contigs is being built\n")
	o2 = os.system("bowtie2-build %s 'ray_selected_contigs_index'"%(size_chosen_ray))
	if(o2==0):
		print("the mapping is starting")
		os.system("bowtie2 -f --very-sensitive-local --al %s --un %s -x 'ray_selected_contigs_index' -U %s"%(map2,unmap2,blastout))
	else:
		print("an error occured, the index couldn't be created")
		raise TaskError
	countMap = pexpect.run("count.sh %s"%(map))
	countRead = pexpect.run("count.sh %s"%(blastout))
	countMapPercent = int(countMap)/int(countRead)
	print (str(countMapPercent)+'% of the viral reads have mapped against the selected contigs from spades\n')
	countMap2 = pexpect.run("count.sh %s"%(map2))
	countRead2 = pexpect.run("count.sh %s"%(blastout))
	countMapPercent2 = int(countMap2)/int(countRead2)
	print (str(countMapPercent2)+'% of the viral reads have mapped against the selected contigs from ray')
	


def read_map_contig(size_chosen_spades,size_chosen_ray,blastout,index1,map,unmap,map2,unmap2,index2):
	print("The viral reads are about to be mapped against the selected contigs from spades\n")
	print("The index from the selected contigs is being built\n")
	o1 = os.system("bowtie2-build %s %s"%(size_chosen_spades,index1))
	if(o1==0):
		print("the mapping is starting")
		os.system("bowtie2 -f --very-sensitive-local --al %s --un %s -x %s -U %s"%(map,unmap,index1,blastout))
	else:
		print("an error occured, the index couldn't be created")
		raise TaskError
	print("The viral reads are about to be mapped against the selected contigs from ray\n")
	print("The index from the selected contigs is being built\n")
	o2 = os.system("bowtie2-build %s %s"%(size_chosen_ray,index2))
	if(o2==0):
		print("the mapping is starting")
		os.system("bowtie2 -f --very-sensitive-local --al %s --un %s -x %s -U %s"%(map2,unmap2,index2,blastout))
	else:
		print("an error occured, the index couldn't be created")
		raise TaskError
	countMap = pexpect.run("count.sh %s"%(map))
	countRead = pexpect.run("count.sh %s"%(blastout))
	countMapPercent = int(countMap)/int(countRead)
	print (str(countMapPercent)+'% of the viral reads have mapped against the selected contigs from spades\n')
	countMap2 = pexpect.run("count.sh %s"%(map2))
	countRead2 = pexpect.run("count.sh %s"%(blastout))
	countMapPercent2 = int(countMap2)/int(countRead2)
	print (str(countMapPercent2)+'% of the viral reads have mapped against the selected contigs from ray\n')



def size_category(filein,size_category1,size_category2,size_category3,size_category4,size_category5,blast_db,blast_format,blast_contig_size_chosen,taxo_contig_size_chosen1,taxo_contig_size_chosen2,taxo_contig_size_chosen3,taxo_contig_size_chosen4,taxo_contig_size_chosen5):
	print("\nthe assembled contigs are considered throught 5 category regarding their sizes\n")
	handle1=open(size_category1,'w')
	handle2=open(size_category2,'w')
	handle3=open(size_category3,'w')
	handle4=open(size_category4,'w')
	handle5=open(size_category5,'w')
	count1=count2=count3=count4=count5=0
	for record in SeqIO.parse(open(filein, "r"), "fasta"):
		if(len(record.seq) < 100):
			handle1.write(">"+str(record.description)+"\n")
			handle1.write(str(record.seq)+"\n")
			count1 = count1 + 1
		else:
			if(len(record.seq) < 250):
				handle2.write(">"+str(record.description)+"\n")
				handle2.write(str(record.seq)+"\n")
				count2 = count2 + 1
			else:
				if(len(record.seq) < 500):
					handle3.write(">"+str(record.description)+"\n")
					handle3.write(str(record.seq)+"\n")
					count3 = count3 + 1
				else:
					if(len(record.seq) < 1000):
						handle4.write(">"+str(record.description)+"\n")
						handle4.write(str(record.seq)+"\n")
						count4 = count4 + 1
					else:
						handle5.write(">"+str(record.description)+"\n")
						handle5.write(str(record.seq)+"\n")
						count5 = count5 + 1
	handle1.close()
	handle2.close()
	handle3.close()
	handle4.close()
	handle5.close()
	print ("there is "+str(count1)+" contigs with a size lower than 100")
	if os.path.getsize(size_category1)>0:
		print("the contigs with a size smaller than 100 are being identified (blast)")
		os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_category1, blast_format, blast_contig_size_chosen))
	if os.path.exists(blast_contig_size_chosen):
		xmlBeta(blast_contig_size_chosen,taxo_contig_size_chosen1)
	print ("there is "+str(count2)+" contigs with a size belonging to [100;250[")
	if os.path.getsize(size_category2)>0:
		print("the contigs with a size belonging to [100;250[ are being identified (blast)")
		os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_category2, blast_format, blast_contig_size_chosen))
	if os.path.exists(blast_contig_size_chosen):
		xmlBeta(blast_contig_size_chosen,taxo_contig_size_chosen2)
	print ("there is "+str(count3)+" contigs with a size belonging to [250;500[")
	if os.path.getsize(size_category3)>0:
		print("the contigs with a size belonging to [250;500[ are being identified (blast)")
		os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_category3, blast_format, blast_contig_size_chosen))
	if os.path.exists(blast_contig_size_chosen):
		xmlBeta(blast_contig_size_chosen,taxo_contig_size_chosen3)
	print ("there is "+str(count4)+" contigs with a size belonging to [500;1000[")
	if os.path.getsize(size_category4)>0:
		print("the contigs with a size belonging to [500;1000[ are being identified (blast)")
		os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_category4, blast_format, blast_contig_size_chosen))
	if os.path.exists(blast_contig_size_chosen):
		xmlBeta(blast_contig_size_chosen,taxo_contig_size_chosen4)
	print ("there is "+str(count5)+" contigs with a size over 1000")
	if os.path.getsize(size_category5)>0:
		print("the contigs with a size over 1000 are being identified (blast)")
		os.system("blastn -max_target_seqs 4 -db %s -query %s -outfmt %s -out %s -evalue 0.0001"%(blast_db, size_category5, blast_format, blast_contig_size_chosen))
	if os.path.exists(blast_contig_size_chosen):
		xmlBeta(blast_contig_size_chosen,taxo_contig_size_chosen5)




@files(job10_input1,job10_input2,job10_input3,job10_input4,job10_input5,job10_input6,job10_input7,job10_input8,job10_input9,job10_input10,job10_input11,job10_input12,job10_input13,job10_input14,job10_output1,job10_output2,job10_output3,job10_output4,job10_output5,job10_output6,job10_output7,job10_output8,job10_output9,job10_output10,job10_output11,job10_output12,job10_output13,job10_output14,job10_output15,job10_output16,job10_output17,job10_output18,job10_output19,job10_output20,job4_output1,job4_output10,job7_key1,job10_output21,job10_output22,job10_output23,job10_output24,job10_output25,job10_output26,job10_output27,job10_output28,job10_output29,job10_output30,job10_output31,    job10_output32,job10_output33,job10_output34,job10_output35,job10_output36,job10_output37,job10_output38,job10_output39,job10_output40,job10_output41)
def category(size_spades1,size_spades2,size_spades3,size_spades4,size_spades5,size_ray1,size_ray2,size_ray3,size_ray4,size_ray5,spades_out,contig_spades,ray_dir,ray_contig,map_spades1,unmap_spades1,map_ray1,unmap_ray1,map_spades2,unmap_spades2,map_ray2,unmap_ray2,map_spades3,unmap_spades3,map_ray3,unmap_ray3,map_spades4,unmap_spades4,map_ray4,unmap_ray4,map_spades5,unmap_spades5,map_ray5,unmap_ray5,blast_db,blastout,blast_format,blast_contig_size_chosen,taxo_contig_size_chosen_spades1,taxo_contig_size_chosen_spades2,taxo_contig_size_chosen_spades3,taxo_contig_size_chosen_spades4,taxo_contig_size_chosen_spades5,taxo_contig_size_chosen_ray1,taxo_contig_size_chosen_ray2,taxo_contig_size_chosen_ray3,taxo_contig_size_chosen_ray4,taxo_contig_size_chosen_ray5,    spades_index_100,spades_index_100_250,spades_index_250_500,spades_index_500_1000,spades_index_1000,ray_index_100,ray_index_100_250,ray_index_250_500,ray_index_500_1000,ray_index_1000):
	file_in_spades = os.path.join(spades_out,contig_spades)
	file_ray = os.path.join(ray_dir,ray_contig)
	size_category(file_in_spades,size_spades1,size_spades2,size_spades3,size_spades4,size_spades5,blast_db,blast_format,blast_contig_size_chosen,taxo_contig_size_chosen_spades1,taxo_contig_size_chosen_spades2,taxo_contig_size_chosen_spades3,taxo_contig_size_chosen_spades4,taxo_contig_size_chosen_spades5)
	size_category(file_ray,size_ray1,size_ray2,size_ray3,size_ray4,size_ray5,blast_db,blast_format,blast_contig_size_chosen,taxo_contig_size_chosen_ray1,taxo_contig_size_chosen_ray2,taxo_contig_size_chosen_ray3,taxo_contig_size_chosen_ray4,taxo_contig_size_chosen_ray5)
	if os.path.getsize(size_spades1)>0:
		if os.path.getsize(size_ray1)>0:
			read_map_contig(size_spades1,size_ray1,blastout,spades_index_100,map_spades1,unmap_spades1,map_ray1,unmap_ray2,ray_index_100)
	else:
		print("can't create the index")
	if os.path.getsize(size_spades2)>0:
		if os.path.getsize(size_ray2)>0:
			read_map_contig(size_spades2,size_ray2,blastout,spades_index_100_250,map_spades2,unmap_spades2,map_ray2,unmap_ray2,ray_index_100_250)
	else:
		print("can't create the index")
	if os.path.getsize(size_spades3)>0:
			if os.path.getsize(size_ray3)>0:
				read_map_contig(size_spades3,size_ray3,blastout,spades_index_250_500,map_spades3,unmap_spades3,map_ray3,unmap_ray3,ray_index_250_500)
	else:
		print("can't create the index")
	if os.path.getsize(size_spades4)>0:
		if os.path.getsize(size_ray4)>0:
			read_map_contig(size_spades4,size_ray4,blastout,spades_index_500_1000,map_spades4,unmap_spades4,map_ray4,unmap_ray4,ray_index_500_1000)
	else:
		print("can't create the index")
	if os.path.getsize(size_spades5)>0:
		if os.path.getsize(size_spades5)>0:
			read_map_contig(size_spades5,size_ray5,blastout,spades_index_1000,map_spades5,unmap_spades5,map_ray5,unmap_ray5,ray_index_1000)
	else:
		print("can't create the index")
		
		
		
		
'''@files(job10_input1,job10_input2,job10_output1,job10_output2,job10_output3,job10_output4)
@follows(plot_contig)
def no_duplicates(spades_doublon,ray_doublon,temp1,temp2,spades_set,ray_set):
	print("the file containing all the viral genotype (after spades assembling) is being created...")
	coverage_taxo(spades_doublon,temp1,temp2,spades_set)
	print("the file containing all the viral genotype (after ray assembling) is being created...")
	coverage_taxo(ray_doublon,temp1,temp2,ray_set)
	
	
	
@files(job11_input1,job11_input2,job11_output1,job11_output2)
@follows(no_duplicates)
def retrieve_query_with_hits(spades_in,ray_in,spades_out,ray_out):
	print("the queries and the matching hits are being retrieved from the alignment made from the spades contig file")
	query(spades_in,spades_out)
	print("the queries and the matching hits are being retrieved from the alignment made from the ray contig file")
	query(ray_in,ray_out)
	
	
@files(job12_input1,job12_input2,job12_input3,job12_input4,job12_output1,job12_output2)
@follows(retrieve_query_with_hits)	
def final_with_duplicates(spades_no_duplicates,ray_no_duplicates,spades_query,ray_query,spades_coverage,ray_coverage):
	print("[spades]the file with viral genotypes and all the matching contigs and their converages is being created...")
	coverage(spades_no_duplicates,spades_query,spades_coverage)
	print("[ray]the file with viral genotypes and all the matching contigs and their converages is being created...")
	coverage(ray_no_duplicates,ray_query,ray_coverage)
	
@files(job13_input1,job13_input2,job13_output1,job13_output2)
@follows(final_with_duplicates)	
def final_without_duplicates(spades_coverage,ray_coverage,spades_final,ray_final):
	print("duplicates are being removed from the final taxonomic file [spades]...")
	sort_contigs(spades_coverage,spades_final)	
	print("duplicates are being removed from the final taxonomic file [ray]...")
	sort_contigs(ray_coverage,ray_final)'''


@files(job3_output1,job8_output1,job2_input1)
#@follows(taxo)		
def pie(map,virus,others):
	print("the distribution of the organisms is about to printed in a pie chart \n")
	count_host = int(pexpect.run("count.sh %s"%(map)))
	count_virus = int(pexpect.run("count.sh %s"%(virus)))
	#count_bacteria = pexpect.run("count.sh %s"%bact)
	count_tot = int(pexpect.run("count_fastq.sh %s"%(others)))
	count_others = count_tot-(count_host+count_virus)
	#list_in=(map,virus,bact,others)
	list_in=(map,virus,others)	
	#list_out=(count_host,count_virus, count_bacteria, count_others)
	list_out=(count_host,count_virus, count_others)
	'''iterator = list_out.__iter__()
	for file in list_in:
		liste_out[i] = pexpect.run("count.sh %s"%file)
		i++
		
	for it in list_in:
		print it
	#traitement des donnees pour affichage dans le graphe'''
	# make a square figure and axes
	figure(1, figsize=(6,6))
	ax = axes([0.1, 0.1, 0.8, 0.8])
	# The slices will be ordered and plotted counter-clockwise.
	#labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
	#labels = 'host', 'virus', 'bacteria', 'others'
	labels = 'host', 'virus', 'others'
	fracs = [count_host, count_virus, count_others]
	#explode=(0, 0.05, 0, 0)
	explode=(0.05, 0.05, 0.075, 0.075)
	
	pie(fracs, explode=explode, labels=labels,
	                autopct='%1.1f%%', shadow=True, startangle=90)
	                # The default startangle is 0, which would start
	                # the Frogs slice on the x-axis.  With startangle=90,
	                # everything is rotated counter-clockwise by 90 degrees,
	                # so the plotting starts on the positive y-axis.
	
	title('distributtion of organisms within the sample', bbox={'facecolor':'0.8', 'pad':5})
	
	show()
	
	
	
	
	
def main():
	'''#bankChoice(job1_input1,job1_input3,job1_input2,job1_output2,job1_output1)
	#blast(job4_input1,job4_input2,job4_input3,job4_input4,job4_output1,job4_output2,job4_key1,job4_output3,job4_key2,job4_output4,job4_input5,job4_output5,job4_output6,job4_output7,job4_output8,job4_output9,job4_output10)
	#bowtie2(job3_input1,job3_input2,job3_output1,job3_output2)
	pipeline_printout(sys.stdout, [category])
	pipeline_run([category])
	#bankChoice()'''
	category(job10_input1,job10_input2,job10_input3,job10_input4,job10_input5,job10_input6,job10_input7,job10_input8,job10_input9,job10_input10,job10_input11,job10_input12,job10_input13,job10_input14,job10_output1,job10_output2,job10_output3,job10_output4,job10_output5,job10_output6,job10_output7,job10_output8,job10_output9,job10_output10,job10_output11,job10_output12,job10_output13,job10_output14,job10_output15,job10_output16,job10_output17,job10_output18,job10_output19,job10_output20,job4_output1,job4_output10,job7_key1,job10_output21,job10_output22,job10_output23,job10_output24,job10_output25,job10_output26,job10_output27,job10_output28,job10_output29,job10_output30,job10_output31,    job10_output32,job10_output33,job10_output34,job10_output35,job10_output36,job10_output37,job10_output38,job10_output39,job10_output40,job10_output41)
	
if __name__ == '__main__':
	main()