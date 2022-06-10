# -----------------------------------------
# Title: Reference_Analyzer_16s.py
# Author: Silver A. Wolf
# Last Modified: Fri, 10.06.2022
# Version: 0.1.5
# -----------------------------------------

import csv
import glob
import os
import numpy as np
import multiprocessing
from multiprocessing import Pool

def analyze_samples(read1):
	read2 = read1.split("_R1.fastq.gz")[0] + "_R2.fastq.gz"
	sample_name = read1.split("/")[1].split("_")[0]
	output_dir = "output/" + sample_name

	print("Analyzing: " + sample_name)

	os.system("join_paired_ends.py " +
			  "-m fastq-join " +
			  "-j 6 " +
			  "-p 8 " +
			  "-f " + read1 + " " +
			  "-r " + read2 + " " +
			  "-o " + output_dir
			 )

	os.system("seqkit seq " +
			  "-m 300 " +
			  "-M 470 " +
			  "-j 4 " +
			  output_dir + "/fastqjoin.join.fastq " +
			  "> " + output_dir + "/fastqjoin.filtered.fastq"
			 )

	os.system("split_libraries_fastq.py " +
			  "-r 3 " +
			  "-n 0 " +
			  "-p 0.5 " +
			  "-q 19 " +
			  "--barcode_type not-barcoded " +
			  "--sample_ids " + sample_name + " " +
			  "-i " + output_dir + "/fastqjoin.filtered.fastq " +
			  "-o " + output_dir
			 )

	os.system("./usearch11.0.667_i86linux32 " +
			  "-otutab " + output_dir + "/seqs.fna " +
			  "-otus databases/gg_13_8_otus/rep_set/99_otus.fasta " +
			  "-id 0.97 " +
			  "-strand both " +
			  "-sizeout " +
			  "-threads 4 " +
			  "-biomout " + output_dir + "/otus.biom " +
			  "-mapout " + output_dir + "/otus_map.txt " +
			  "-otutabout " + output_dir + "/otus_about.txt " +
			  "-dbmatched " + output_dir + "/otus_with_sizes.fa " +
			  "-notmatched " + output_dir + "/otus_notmatched.fa"
			 )

	os.system("./usearch11.0.667_i86linux32 " +
			  "-sortbysize " + output_dir + "/otus_with_sizes.fa " +
			  "-fastaout " + output_dir + "/otus_sorted.fa " +
			  "-minsize 0"
			 )

	os.system("./usearch11.0.667_i86linux32 " +
			  "-unoise3 " + output_dir + "/otus_sorted.fa " +
			  "-tabbedout " + output_dir + "/zotus.tsv " +
			  "-zotus " + output_dir + "/zotus.fa " +
			  "-minsize 0"
			 )

	#os.system("SequenceMatch " +
	#		  "train databases/current_Bacteria_unaligned.fa databases/rdp/rdp"
	#		 )

	os.system("SequenceMatch " +
			  "seqmatch databases/rdp/ " + output_dir + "/zotus.fa " +
			  "-k 1 " +
			  "-s 0.5 " +
			  "> " + output_dir + "/seq_results.txt"
			 )

	# Continuing analysis of the unknown sequences not found in the GreenGenes database
	os.system("./usearch11.0.667_i86linux32 " +
			  "-fastx_uniques " + output_dir + "/otus_notmatched.fa " +
			  "-fastaout " + output_dir + "/otus_notmatched_derep.fa " +
			  "-relabel " + sample_name + " " +
			  "-sizeout " +
			  "-strand both " +
			  "-threads 4"
			 )
	
	os.system("./usearch11.0.667_i86linux32 " +
			  "-unoise3 " + output_dir + "/otus_notmatched_derep.fa " +
			  "-tabbedout " + output_dir + "/zotus_notmatched.tsv " +
			  "-zotus " + output_dir + "/zotus_notmatched.fa " +
			  "-minsize 0"
			 )
	
	os.system("SequenceMatch " +
			  "seqmatch databases/rdp/ " + output_dir + "/zotus_notmatched.fa " +
			  "-k 1 " +
			  "-s 0.5 " +
			  "> " + output_dir + "/seq_results_notmatched.txt"
			 )

	print("Finished: " + sample_name)

def create_biom_table():
	dirs = os.listdir("output")
	tax_ids = []

	# For each sample
	# Open its seq_results file and unmapped file
	# Create list of taxids without any duplicates
	for d in dirs:
		rdp_results = "output/" + d + "/seq_results.txt"
		with open(rdp_results, "r") as infile:
			for line in infile:
				split = line.split("\t")
				if ((split[1] not in tax_ids) and (split[0] != "query name")):
					tax_ids.append(split[1])
		rdp_results_unmapped = "output/" + d + "/seq_results_notmatched.txt"
		with open(rdp_results_unmapped, "r") as infile:
			for line in infile:
				split = line.split("\t")
				if ((split[1] not in tax_ids) and (split[0] != "query name")):
					tax_ids.append(split[1])

	tab = open("output/final_otu_table.tsv", "w")
	tab.write("# Constructed from biom file\n")
	tab.write("#OTU ID" + "\t" + "\t".join(dirs) + "\t" + "taxonomy\n")

	# For each taxid in list
	for t in tax_ids:
		final_line = t
		# Go through each samples seq_results files
		for d in dirs:
			# Set abundance to 0
			abundance = 0
			
			# Scan through seq_results
			rdp_results = "output/" + d + "/seq_results.txt"
			with open(rdp_results, "r") as infile:
				for line in infile:
					split = line.split("\t")
					# If id in file -> retrieve abundance from tsv
					if split[1] == t:
						c = 0
						number = int(split[0].split("Zotu")[1])
						zotu_results = "output/" + d + "/zotus.tsv"
						with open(zotu_results, "r") as outfile:
							for l in outfile:
								s = l.split("\t")
								if "zotu" in s[2]:
									c = c + 1
									if number == c:
										abundance = abundance + int(s[0].split("size=")[1].split(";")[0])
										break
			
			# Scan through seq_results_unmapped
			rdp_results_unmapped = "output/" + d + "/seq_results_notmatched.txt"
			with open(rdp_results_unmapped, "r") as infile:
				for line in infile:
					split = line.split("\t")
					# If id in file -> retrieve abundance from tsv
					if split[1] == t:
						c = 0
						number = int(split[0].split("Zotu")[1])
						zotu_results = "output/" + d + "/zotus_notmatched.tsv"
						with open(zotu_results, "r") as outfile:
							for l in outfile:
								s = l.split("\t")
								if "zotu" in s[2]:
									c = c + 1
									if number == c:
										abundance = abundance + int(s[0].split("size=")[1].split(";")[0])
										break
			
			final_line = final_line + "\t" + str(abundance)

		# Once done with the abundancies, retrieve taxa information from reference fasta
		rdp_database = "databases/current_Bacteria_unaligned.fa"
		with open(rdp_database, "r") as db:
			domain = "k__"
			phylum = "p__"
			class_tax = "c__"
			order = "o__"
			family = "f__"
			genus = "g__"
			species = "s__"
			for entry in db:
				if t in entry:
					if ";domain" in entry:
						domain = domain + entry.split(";domain")[0].split(";")[-1].replace("\"", "")
					if ";phylum" in entry:
						phylum = phylum + entry.split(";phylum")[0].split(";")[-1].replace("\"", "")
					if ";class" in entry:
						class_tax = class_tax + entry.split(";class")[0].split(";")[-1].replace("\"", "")
					if ";order" in entry:
						order = order + entry.split(";order")[0].split(";")[-1].replace("\"", "")
					if ";family" in entry:
						family = family + entry.split(";family")[0].split(";")[-1].replace("\"", "")
					if ";genus" in entry:
						genus = genus + entry.split(";genus")[0].split(";")[-1].replace("\"", "")
					if ";species" in entry:
						species = species + entry.split(";species")[0].split(";")[-1].replace("\"", "")
					tax = domain + "; " + phylum + "; " + class_tax + "; " + order + "; " + family + "; " + genus + "; "+ species
					break

		# Write line with abundancies per sample plus taxa information into tsv file
		final_line = final_line + "\t" + tax + "\n"
		tab.write(final_line)

	tab.close()

	# Update phyla labels according to newest taxonomy
	update_phyla()

	# Export tsv to BIOM format for downstream analysis
	os.system("biom convert " +
			  "-i output/final_otu_table.tsv " +
			  "-o output/final_otu_table.biom " +
			  "--to-hdf5 " +
			  "--table-type=\"OTU table\" " +
			  "--process-obs-metadata taxonomy"
			 )

	os.system("filter_otus_from_otu_table.py " +
			  "-i output/final_otu_table.biom " +
			  "-o output/final_otu_table_clean.biom " +
			  "-n 4"
			 )
	
	os.system("biom convert " +
			  "-i output/final_otu_table_clean.biom " +
			  "-o output/final_otu_table_clean.tsv " +
			  "--to-tsv " +
			  "--header-key taxonomy"
			 )

def update_phyla():
	dict_tax = {}
	
	with open("databases/42_Bacterial_Phyla_2022.csv", mode = "r") as tax:
		full_tax = csv.reader(tax)
		dict_tax = {rows[0]:rows[1] for rows in full_tax}
	
	del dict_tax["Old"]
	
	for key in dict_tax:
		old_name = "p__" + key
		new_name = "p__" + dict_tax[key]
		os.system("sed -i -- \'s/" + old_name + "/" + new_name + "/g\' output/final_otu_table.tsv")

def helperfunction(list_reads):
	for read in list_reads:
		analyze_samples(read)

# Main
def main():
	seqs = glob.glob("sequences/*_R1.fastq.gz")
	
	# Parallelize Script
	amount_threads = 54
	list_reads = np.array_split(seqs, amount_threads)
	pool = Pool(amount_threads)
	results = pool.map(helperfunction, list_reads)
	
	create_biom_table()

if __name__ == "__main__":
	main()