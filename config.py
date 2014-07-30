#!/usr/bin/python -tt
import ConfigParser
import os

cfg = ConfigParser.ConfigParser()
cfg.add_section('data_input')
cfg.add_section('Filename')

#cfg2 = ConfigParser.ConfigParser()
#cfg2.read('fichierdeconfig2.cfg')
#directory = cfg2.get('directory', 'cle1')

S = 'data_input'
cfg.set(S, 'cle1', 'Human_papillomavirus')#key word indicating the organism that is looked for
cfg.set(S, 'cle2', '')#the primer sequence
cfg.set(S, 'cle3', '0.0001')#cut off value for blast query
cfg.set(S, 'cle4', '6')#format of the blast output file
cfg.set(S, 'cle5', '31')#kmer size
cfg.set(S, 'cle6', '6')
cfg.set(S, 'cle7', '5')#taxo blast ouptut format


S = 'Filename'
cfg.set(S, 'cle1', 'gbvrl_nucleic.fst')#file used for the previous function creating the data bank
cfg.set(S, 'cle2', '394_S4_L001_TOTAL.fastq')#raw data (read) in fastq format
cfg.set(S, 'cle2(0)', 'quality_seq.fastq')#data without high quality reads
cfg.set(S, 'cle2(1)', 'trimmed.fastq')#trimmed reads
cfg.set(S, 'cle2(12)', 'trimmed.fasta')# trimmed reads in fasta format
cfg.set(S, 'cle2(13)', 'size.txt')#text file containing reads size
cfg.set(S, 'cle3', 'ICTV.xls')#ICTV viral classification
cfg.set(S, 'cle4', 'bank_nucl.fasta')#fasta file containing nucleotide sequences of the organism choosed from viral genbank
cfg.set(S, 'cle4(0)', 'bank_prot.fasta')#fasta file containing protein sequences of the organism choosed from viral genbank
cfg.set(S, 'cle4(2)', 'final_contigs.fasta')#test file for the contigs plot
cfg.set(S, 'cle4(1)', 'viral_nuc.fasta')#file without quote that will be used as a blast database
cfg.set(S, 'cle4(3)', 'viral_prot.fasta')#file without quote that will be used as a blast database
cfg.set(S, 'cle5', 'gbvrl.seq')#genbank file carrying the nucleotide sequences (used to create the data bank)
cfg.set(S, 'cle5(0)', 'gbvrl_prot.txt')#genbank file carrying the protein sequences (used to create the data bank)
cfg.set(S, 'cle6', 'contigs.fasta')
cfg.set(S, 'cle7', 'hg19')#reference genome
cfg.set(S, 'cle8', 'unmatch.fasta')#unmatched reads
cfg.set(S, 'cle9', 'match.fasta')#matched reads
cfg.set(S, 'cle10', 'hg19_index')#reference index
cfg.set(S, 'cle11', 'unmatch.fasta')#converted file after mapping
cfg.set(S, 'cle11(1)', 'IonXpress_oreillon.fa')#test file with hiseq reads before blast
cfg.set(S, 'cle11(2)', 'trimmed.fasta')#trimmed reads after fastq to fasta conversion
cfg.set(S, 'cle12', 'match.fasta')#converted file after mapping
cfg.set(S, 'cle13', 'viral_db')#the database (used in the blast query) made from the data bank
cfg.set(S, 'cle15', 'blast_output.fasta')#aligned reads from blast results (positions)
cfg.set(S, 'cle17', 'awk_output.fasta')
cfg.set(S, 'cle20','final_blast_output.fasta')#final blast output/assemmbly input
cfg.set(S, 'cle21', 'spades_out')#spades output directory
cfg.set(S, 'cle21(1)', 'K55')#spades output subdirectory
cfg.set(S, 'cle22', 'ray_out')
cfg.set(S, 'cle22(1)', 'Contigs.fasta')#ray contigs
cfg.set(S, 'cle23', 'taxo_out.fasta')
#cfg.set(S, 'cle24', os.path.join(directory,'testbonrepertoire.txt'))
cfg.set(S, 'cle25', 'test_bwa')#test sur bwa aln pour le pb de version
cfg.set(S, 'cle26', 'map.fasta')#bowtie file to find the used reads in the assembly
cfg.set(S, 'cle27', 'unmap.fasta')#bowtie file to find the unused reads in the assembly
cfg.set(S, 'cle28', 'contigs.fasta')#spades contig file copied in the algo directory
cfg.set(S, 'cle29', 'assembly_index')#bowtie index file made from spades contig
cfg.set(S, 'cle30', 'blastn_neg.fasta')##reads that were not aligned by blast
cfg.set(S, 'cle31', 'blastx_output.fasta')
cfg.set(S, 'cle32', 'awk_out_blastx.fasta')
cfg.set(S, 'cle33', 'viral_blastx_db')
cfg.set(S, 'cle34', 'blastnout_50.fasta')
cfg.set(S, 'cle35', 'blastxout_50.fasta')
cfg.set(S, 'cle36', 'singleton_spades.fasta')#singletons found by mapping against positive blast after spades assembly
cfg.set(S, 'cle37', 'singleton_ray.fasta')#singletons found by mapping against positive blast after ray assembly
cfg.set(S, 'cle38', 'entire_viral_spades.fasta')#the whole set of viral variants (singleton and blast positive) (spades)
cfg.set(S, 'cle39', 'entire_viral_ray.fasta')#the whole set of viral variants (singleton and blast positive) (ray)
cfg.set(S, 'cle40', 'blast_final_taxo_spades.xml')#spades taxonomic file (blast)
cfg.set(S, 'cle41', 'blast_final_taxo_ray.xml')#ray taxonomic file (blast)
cfg.set(S, 'cle42', 'hpv_genotype_spades.txt')#file containing the blast output after xml parsing[spades]
cfg.set(S, 'cle43', 'hpv_genotype_ray.txt')#file containing the blast output after xml parsing[ray]
cfg.set(S, 'cle44', 'spades_contigs_chosen_size.fasta')#[sapdes]file containing the contigs with relevant sizes
cfg.set(S, 'cle45', 'ray_contigs_chosen_size.fasta')#[ray]file containing the contigs with relevant sizes
cfg.set(S, 'cle46', 'blast_spades_contigs_chosen_size.xml')#[spades]blast output (xml) from the chosen (by size) contigs
cfg.set(S, 'cle47', 'blast_ray_contigs_chosen_size.xml')#[ray]blast output (xml) from the chosen (by size) contigs
cfg.set(S, 'cle48', 'taxo__spades_contigs_chosen_size.txt')#[spades]taxonomic designation from the chosen (by size) contigs
cfg.set(S, 'cle49', 'taxo_ray_contigs_chosen_size.txt')#[ray]taxonomic designation from the chosen (by size) contigs
cfg.set(S, 'cle50', 'temporary1.txt')#temporary file for the removal of duplicates in viral genotypes 
cfg.set(S, 'cle51', 'temporary2.txt')#temporary file for the removal of duplicates in viral genotypes 
cfg.set(S, 'cle52', 'spades_taxo_set.txt')#file containing the viral genotypes without duplicates after spades assembling
cfg.set(S, 'cle53', 'ray_taxo_set.txt')#file containing the viral genotypes without duplicates after ray assembling
cfg.set(S, 'cle54', 'spades_queries.txt')#file with all the queries and the matching hits [spades]
cfg.set(S, 'cle55', 'ray_queries.txt')#file with all the queries and the matching hits [ray]
cfg.set(S, 'cle56', 'spades_genotype_coverage.fasta')#file with the viral genotypes and all the contigs [spades]
cfg.set(S, 'cle57', 'ray_genotype_coverage.fasta')#file with the viral genotypes and all the contigs[ray]
cfg.set(S, 'cle58', 'spades_genotype_coverage_no_duplicates.fasta')#file with no duplicates in viral genotypes and no duplicates in the contigs [spades]
cfg.set(S, 'cle59', 'ray_genotype_coverage_no_duplicates.fasta')#file with no duplicates in viral genotypes and no duplicates in the contigs [ray]
cfg.set(S, 'cle60', 'spades_read_contig_mapped.fasta')#mapped reads against the selected contigs[spades]
cfg.set(S, 'cle61', 'spades_read_contig_unmapped.fasta')#unmapped reads against the selected contigs[spades]
cfg.set(S, 'cle62', 'ray_read_contig_mapped.fasta')#mapped reads against the selected contigs[ray]
cfg.set(S, 'cle63', 'ray_read_contig_unmapped.fasta')#unmapped reads against the selected contigs[ray]
cfg.set(S, 'cle64', 'spades_mapped_100.fasta')#mapped reads against the selected contigs[spades]
cfg.set(S, 'cle65', 'spades_mapped_100_250.fasta')#mapped reads against the selected contigs[spades]
cfg.set(S, 'cle66', 'spades_mapped_250_500.fasta')#mapped reads against the selected contigs[spades]
cfg.set(S, 'cle67', 'spades_mapped_500_1000.fasta')#mapped reads against the selected contigs[spades]
cfg.set(S, 'cle68', 'spades_mapped_1000.fasta')#mapped reads against the selected contigs[spades]
cfg.set(S, 'cle69', 'spades_unmapped_100].fasta')#unmapped reads against the selected contigs[spades]
cfg.set(S, 'cle70', 'spades_unmapped_100_250.fasta')#unmapped reads against the selected contigs[spades]
cfg.set(S, 'cle71', 'spades_unmapped_250_500.fasta')#unmapped reads against the selected contigs[spades]
cfg.set(S, 'cle72', 'spades_unmapped_500_1000.fasta')#unmapped reads against the selected contigs[spades]
cfg.set(S, 'cle73', 'spades_unmapped_1000.fasta')#unmapped reads against the selected contigs[spades]
cfg.set(S, 'cle74', 'ray_mapped_100.fasta')#mapped reads against the selected contigs[ray]
cfg.set(S, 'cle75', 'ray_mapped_100_250.fasta')#mapped reads against the selected contigs[ray]
cfg.set(S, 'cle76', 'ray_mapped_250_500.fasta')#mapped reads against the selected contigs[ray]
cfg.set(S, 'cle77', 'ray_mapped_500_1000.fasta')#mapped reads against the selected contigs[ray]
cfg.set(S, 'cle78', 'ray_mapped_1000.fasta')#mapped reads against the selected contigs[ray]
cfg.set(S, 'cle79', 'ray_unmapped_100.fasta')#unmapped reads against the selected contigs[ray]
cfg.set(S, 'cle80', 'ray_unmapped_100_250.fasta')#unmapped reads against the selected contigs[ray]
cfg.set(S, 'cle81', 'ray_unmapped_250_500.fasta')#unmapped reads against the selected contigs[ray]
cfg.set(S, 'cle82', 'ray_unmapped_500_1000.fasta')#unmapped reads against the selected contigs[ray]
cfg.set(S, 'cle83', 'ray_unmapped_1000.fasta')#unmapped reads against the selected contigs[ray]
cfg.set(S, 'cle84', 'spades_contig_100.fasta')#contigs size 100[spades]
cfg.set(S, 'cle85', 'spades_contig_100_250.fasta')#contigs size 100-250[spades]
cfg.set(S, 'cle86', 'spades_contig_250_500.fasta')#contigs size 250-500[spades]
cfg.set(S, 'cle87', 'spades_contig_500_1000.fasta')#contigs size 500-1000[spades]
cfg.set(S, 'cle88', 'spades_contig_1000.fasta')#contigs size 1000[spades]
cfg.set(S, 'cle89', 'ray_contig_100.fasta')#contigs size 100[ray]
cfg.set(S, 'cle90', 'ray_contig_100_250.fasta')#contigs size 100-250[ray]
cfg.set(S, 'cle91', 'ray_contig_250_500.fasta')#contigs size 250-500[ray]
cfg.set(S, 'cle92', 'ray_contig_500_1000.fasta')#contigs size 500-1000[ray]
cfg.set(S, 'cle93', 'ray_contig_1000.fasta')#contigs size 1000[ray]
cfg.set(S, 'cle94', 'out_multi_taxo.xml')#contigs size 1000[ray]
cfg.set(S, 'cle95', 'taxo_contig_spades_100.txt')#contigs size 100[spades]
cfg.set(S, 'cle96', 'taxo_contig_spades_100_250.txt')#contigs size 100-250[spades]
cfg.set(S, 'cle97', 'taxo_contig_spades_250_500.txt')#contigs size 250-500[spades]
cfg.set(S, 'cle98', 'taxo_contig_spades_500_1000.txt')#contigs size 500-1000[spades]
cfg.set(S, 'cle99', 'taxo_contig_spades_1000.txt')#contigs size 1000[spades]
cfg.set(S, 'cle100', 'taxo_contig_ray_100.txt')#contigs size 100[ray]
cfg.set(S, 'cle101', 'taxo_contig_ray_100_250.txt')#contigs size 100-250[ray]
cfg.set(S, 'cle102', 'taxo_contig_ray_250_500.txt')#contigs size 250-500[ray]
cfg.set(S, 'cle103', 'taxo_contig_ray_500_1000.txt')#contigs size 500-1000[ray]
cfg.set(S, 'cle104', 'taxo_contig_ray_1000.txt')#contigs size 1000[ray]






cfg.write(open('fichierdeconfig.cfg','w'))