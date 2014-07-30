
from Bio import SeqIO



def fastaToFasta():
	count = 0
	file_name=raw_input("donnez input\n")
	file_out=raw_input("donner output\n")
	handle = open(file_out,"w")
	for seq_record in SeqIO.parse(file_name, "fasta"):
		 handle.write(">"+str(seq_record.description)+"\n")
		 handle.write(str(seq_record.seq+"\n"))
		 count = count + 1
	handle.close()
	return count


def main():
	fastaToFasta()
	
if __name__ == '__main__':
	main()