
#!/usr/bin/python

##This is in-house script that was used for read mapping of vMAGs and MAGs using the parameters mentioned in steps 2.0 and 3.0

#PYTHON SCRIPT 
#written by: Richard Wolfe
#
#to run type: python sam_file_analize_reads.py -i <inputfile> 
#         or: ./ongest_sequence.py -i <inputfile> 
#
#   if error: /usr/bin/python^M: bad interpreter: No such file or directory
#      -there is a windows endl after shebang
#      -open in vi 
#         once in vi type:
#           :set ff=unix<return>
#           :x<return>
#
#
#
# goes through the sam file and analyze how the reads mapped
# 
#   Note:
#



import sys      #for exit command and maxint
import argparse #to get command line args 
                #needed to install argparse module because using python 2.6
                #and argparse comes with python 2.7
                #  sudo easy_install argparse
import os.path
import glob     #for * in file name
import toolbox

import count_multi_map_reads   #to make sam_stats table

		
def check_fastq_lines(R1,R2):
	#R1 and R2 are file names

	"""
	f1 = open(R1, "r")
	f2 = open(R2, "r")

	f1_lines = 0
	f2_lines = 0

	line = f1.readline()  
	while line:
		f1_lines += 1
		line = f1.readline() 

	line = f2.readline()  
	while line:
		f2_lines += 1
		line = f2.readline()
	"""
	f1_lines = 0
	f2_lines = 0

	with open(R1) as f:
        	for f1_lines, l in enumerate(f):
            		pass

	with open(R2) as f:
        	for f2_lines, l in enumerate(f):
            		pass
	
	#need to add 1 line?
	f1_lines += 1
	f2_lines += 1

	if f1_lines != f2_lines:
		print "Error....Files are not same length"
		print R1 + " lines = ", f1_lines
		print R2 + " lines = ", f2_lines
		sys.exit(1)


	#f1.close()
	#f2.close()

	#print "Files are same length"
	#print R1 + " lines = ", f1_lines
	#print R2 + " lines = ", f2_lines
	#sys.exit(1)

	return

def get_from_file(f):
	#
	# f = file handle
	# returns a list of the lines in the file

	l = []

	line = f.readline()
	while line:
		l.append(line.rstrip())   #removes the endline
		line = f.readline()

	return l



#create an argument parser object
#description will be printed when help is used
parser = argparse.ArgumentParser(description='A script to analyze how the reads mapped')

#add the available arguments -h and --help are aqdded by default
#if the input file does not exist then program will exit
#if output file does not exit it will be created
# args.input is the input file Note: cant write to this file because read only
# args.output is the output file
# args.m is the minimum seq length
parser.add_argument('-f', '--fasta_file',type=argparse.FileType('rU'), help='fasta file to map to ', required=True)

#parser.add_argument('-r', '--reads_files', help='comma seperated list of reads files to map(must be 2 and no spaces)', required=True)
#parser.add_argument('-n', '--names', help='comma seperated list of sample names for the reads(no spaces)', required=True)
parser.add_argument('-r', '--reads_files', type=argparse.FileType('rU'), help='text file of R1 and R2 reads files to map(must be 1 file per line)', required=True)
parser.add_argument('-n', '--names', type=argparse.FileType('rU'), help='text file of sample names for the reads(mut be 1 per line)', required=True)

parser.add_argument('-t', '--output',  help='Output file table name. Also, sam_stats_<filename>', required=True)

parser.add_argument('-m', '--max_mismatch', type=int, help='max number of mismatches') #required=True) 
parser.add_argument('--percent_read', type=int, help='percent of (read_length - mismatches)/read_length to keep EX: 90')  #optional 

parser.add_argument('-p', '--processors', type=int, help='number of processors', required=True)

parser.add_argument('--min_coverage', type=float, help='Minimum coverage cutoff Ex:250 or 250.00', required=True)
parser.add_argument('--min_percent_contig_coverage', type=float, help='Minimum percent of contig coverage cutoff Ex:80 or 80.00', required=True)

parser.add_argument('-u', '--multi_map', help='multi map T or F', default="F")
parser.add_argument('-d', '--delete_files', help='delete files made EX: bowtie and sam files T or F', default="T")

#get the args
args = parser.parse_args()

#can only have -m or --percent_read not both
if args.max_mismatch == None: #because int
	if args.percent_read == None: #because int
        	print "Error ... must have -m or --percent_read"
		sys.exit(1)

if not args.max_mismatch == None: #because int
	if not args.percent_read == None: #because int
        	print "Error ... can not have both -m and --percent_read"
		sys.exit(1)

#if using --percent_read set max_mismatch very high
if args.max_mismatch == None: #because int
	args.max_mismatch = 5000

#print args.max_mismatch
#sys.exit(0)

if args.max_mismatch < 0:
	print "Error ... -m must be greater than or equal to 0"
	sys.exit(1)

if args.processors < 1:
	print "Error ... -p must be greater than 0"
	sys.exit(1)

#see if output file exists so it is not overwritten
if os.path.exists(args.output):
	print "Error .... Output file already exists"
	sys.exit(1)
	
if args.multi_map != "T":
	if args.multi_map != "F":
		print "Error ... -u must be T or F"
    		sys.exit(1)

if args.delete_files != "T":
	if args.delete_files != "F":
		print "Error ... -d must be T or F"
    		sys.exit(1)



print "Script Started ... "


#reads = args.reads_files.split(",")
#names = args.names.split(",")

reads = get_from_file(args.reads_files)
names = get_from_file(args.names)

reads_str = ",".join(reads)
names_str = ",".join(names)

#print reads_str
#print names_str
#sys.exit()

#check if reads file is a pair
if len(reads) % 2 != 0:
	print "Error .... reads files not in sets of 2"
	sys.exit(1)


#make sure all files exist
for item in reads:
	if not os.path.exists(item):
		print "Error .... reads file does not exist...." + item
		sys.exit(1)

#check if each pair of reads has a name and not 2 names the same
if len(reads) / 2 != len(names):
	print "Error .... number of reads files divided by 2 not equal to number of names"
	sys.exit(1)
temp = []
for item in names:
	if item in temp:
		print "Error ... Same name occurs more than once"
		sys.exit(1)
	else:
		temp.append(item)
temp = []

#check to make sure reads files have same number of lines
i = 0
while i < len(reads):
	check_fastq_lines(reads[i],reads[i+1])

	i = i + 2


#map each pair of reads to the fasta file
# 0. This script checks to make sure no sam files are in the folder before proceeding
# 1. This runs bowtie2 with -N 1 -- alows mismatches in seed
#    bowtie2 is run on each set of reads mapping to the fasta file
# 2. Then runs script that makes a new sam files by marking reads that map with mismatches > specified  to not map
# 3. Also runs a script that reads all the new sam files and makes a report table named map_reads_fasta_output.txt
#      This table is not used
cmd = "python /ORG-Data/scripts/map_reads_fasta.py -f " + args.fasta_file.name + " -r " + reads_str
cmd = cmd + " -n " + names_str + " -t map_reads_fasta_output.txt " + " -m " + str(args.max_mismatch) 
cmd = cmd + " -p " + str(args.processors) + " -d F -u " + args.multi_map
if not args.percent_read == None: #because int
	cmd = cmd + " --percent_read " + str(args.percent_read)
print "Running command: " + cmd 	
toolbox.run_system(cmd)


for item in names:
	f_name = args.fasta_file.name.split("/")[-1]  #get the file name in case it is a path

	sam_file = "mismatches_" + str(args.max_mismatch) + "_" + item + "_to_" + f_name + ".sam"

	cmd = "python /ORG-Data/scripts/make_contig_cov_file_abundance.py -s " + args.fasta_file.name 
   	cmd = cmd + " -o " + item + "_make_contig_cov_abundance.txt" + " -b " + sam_file

	print "Running command: " + cmd 	
	toolbox.run_system(cmd)

for item in names:
	f_name = args.fasta_file.name.split("/")[-1]  #get the file name in case it is a path

	sam_file = "mismatches_" + str(args.max_mismatch) + "_" + item + "_to_" + f_name + ".sam"

	cmd = "python /ORG-Data/scripts/make_contig_cov_file.py -s " + args.fasta_file.name 
   	cmd = cmd + " -o " + item + "_make_contig_cov_output.txt" + " -b " + sam_file

	print "Running command: " + cmd 	
	toolbox.run_system(cmd)

for item in names:
	f_name = args.fasta_file.name.split("/")[-1]  #get the file name in case it is a path

	sam_file = "mismatches_" + str(args.max_mismatch) + "_" + item + "_to_" + f_name + ".sam"

	cmd = "python /ORG-Data/scripts/calc_percent_contig_cov.py -s " + sam_file 
   	cmd = cmd + " -o " + item + "_cal_percent_contig_cov.txt" + " -v " + str(args.max_mismatch) 
	cmd = cmd + " -m 1 -p 0"

	print "Running command: " + cmd 	
	toolbox.run_system(cmd)




#make a list of genes from the fasta file
gene_ids = []
results = []  #will be a list of lists .. the percent_reads_mapped for each sample

line = args.fasta_file.readline()
while line:
	if line.startswith(">"):
		id = line.split()[0]
		id = id[1:]  #remove >
		if id not in gene_ids:
			gene_ids.append(id)
		else:
			print "Error... duplicate id in fasta file = " + id
			sys.exit(1)
	line = args.fasta_file.readline()



#extract the info from the ouptut files
for item in names:
	percent_reads_mapped = []
	contig_cov = []
	percent_contig_cov = []
	for i in gene_ids:
		percent_reads_mapped.append("NA")
		contig_cov.append("NA")
		percent_contig_cov.append("NA")

	#get the % reads mapped for this sample
	f = open(item + "_make_contig_cov_abundance.txt", "r")
	line = f.readline() #header line
	line = f.readline()  
	while line:
		line = line.rstrip()
		p = line.split("\t")[3]
		id = line.split("\t")[0]

		index = gene_ids.index(id)
		percent_reads_mapped[index] = p

		line = f.readline()

	#make sure everything found
	for i in percent_reads_mapped:
		if i == "NA":
			print "ERROR in " + item + "_make_contig_cov_abundance.txt"
			sys.exit(1)

	#get the contig cov for this sample
	f = open(item + "_make_contig_cov_output.txt", "r")
	
	line = f.readline()  
	while line:
		line = line.rstrip()
		c = line.split("\t")[4]
		id = line.split("\t")[0]

		index = gene_ids.index(id)
		contig_cov[index] = c

		line = f.readline()

	#make sure everything found
	for i in contig_cov:
		if i == "NA":
			print "ERROR in " + item + "_make_contig_cov_output.txt"
			sys.exit(1)

	#get the % contig cov for this sample
	f = open(item + "_cal_percent_contig_cov.txt", "r")
	
	line = f.readline()  
	while line:
		line = line.rstrip()
		pc = line.split("\t")[3]
		id = line.split("\t")[0]

		index = gene_ids.index(id)
		percent_contig_cov[index] = pc

		line = f.readline()

	#make sure everything found
	for i in percent_contig_cov:
		if i == "NA":
			print "ERROR in " + item + "_cal_percent_contig_cov.txt"
			sys.exit(1)

	#see if this % reads mapped meets thresholds and save list
	#parser.add_argument('--min_coverage', type=float, help='Minimum coverage cutoff Ex:250 or 250.00', required=True)
	#parser.add_argument('--min_percent_contig_coverage', type=float, help='Minimum percent of contig coverage cutoff Ex:80 or 80.00', required=True)

	i = 0
	while i < len(gene_ids):
		c = float(contig_cov[i]) * 1.0
		p = float(percent_contig_cov[i]) * 100.0
		
		if c < args.min_coverage or p < args.min_percent_contig_coverage:
			percent_reads_mapped[i] = "0.0"

		i += 1

	results.append(percent_reads_mapped)

#print the output table
output = open(args.output, "w")
output.write("Gene_id\t")
i = 0
while i < len(names):
	output.write(names[i] + "\t")
	i += 1
output.write("\n")

i = 0
while i < len(gene_ids):
	output.write(gene_ids[i] + "\t")
	for item in results:
		output.write(item[i] + "\t")

	output.write("\n")

	i += 1

##################################################
###make sam_stats table for each of the samples

sam_files = []     #the sam files for each sample
sample_names = []  #the sample names (could have used names [])
sam_data = []      #will be a list of tuples with the data from 1 or more sam files

for item in names:
	sam_file_name = "mismatches_" + str(args.max_mismatch) + "_" + item + "_to_" + f_name + ".sam"
        sam_files.append(sam_file_name)
       	sample_names.append(item)

#read each sam file and put the data tuple into sam_data
i = 0
while i < len(sam_files):
	data_tuple = count_multi_map_reads.read_sam_file(sam_files[i])
	sam_data.append(data_tuple)
		
 	i += 1

#print the data to the output file
sam_data_table = open("sam_stats_" + args.output, "w")
count_multi_map_reads.write_table(sam_data_table, sam_data, sample_names)


############################
if args.delete_files == "T":
	for item in names:
			try:
				os.remove(item + "_make_contig_cov_abundance.txt")
				os.remove(item + "_make_contig_cov_output.txt")
				os.remove(item + "_cal_percent_contig_cov.txt")
			except OSError:
				print "Some files for " + item + " could not be removed"
	try:		
		#os.remove(glob('*.sam'))
		#os.remove(glob('*.bt2'))
		cmd = "rm *.sam"  
   		
		print "Running command: " + cmd 	
		toolbox.run_system(cmd)

		cmd = "rm *.bt2"  
   		
		print "Running command: " + cmd 	
		toolbox.run_system(cmd)
		
	except OSError:
		print "Some sam or bowtie files could not be removed"

print ""
print "Number of reads files = ", len(reads)
print "Number of names = ", len(names)
print ""
print "Number of ids in fasta file = ",len(gene_ids)


print "Script finished..."
sys.exit(0)
