#!/usr/bin/python

#PYTHON SCRIPT 
#written by: Richard Wolfe
#
#to run type: python single_copy_genes.py -i <inputfile> -o <outputfile> -m <min sequence length>
#         or: ./single_copy_genes.py attributes
#
#   if error: /usr/bin/python^M: bad interpreter: No such file or directory
#      -there is a windows endl after shebang
#      -open in vi 
#         once in vi type:
#           :set ff=unix<return>
#           :x<return>
#
#
#  runs scg on all of the metabat bins in the folder
#  it make a seperate folder for each bin    

#
# this script runs AMPHORA2 pipeline and finds single copy genes
# 
#   -input file would be a DNA fasta file for example Methanohalophilus_348_143_42_scaffold.fa which
#     was made by using pullseq to pull certain scaffolds from the scaffold.fa file
#   -output file  is the name you want to call the output file
#
#   Example command to run :
#      python single_copy_genes.py -i Methanohalophilus_348_143_42_scaffold.fa -o Methanohalophilus_348_143_42_phylotype_1e-10.result -e 1e-10 -t Archaea -p 20
#

import sys      #for exit command and maxint
import argparse #to get command line args 
                #needed to install argparse module because using python 2.6
                #and argparse comes with python 2.7
                #  sudo easy_install argparse
import os       #to run system commands
import datetime #to make timestamp
import glob     #for * in file name
from string import maketrans #to make translate table for string timetable

import toolbox


print "Script started ..."

def cat_SCG_results():
	#cat all the SCG files and then exit

	#make a list of all the SCG results files
	SCG_dir = glob.glob('SCG_*.metabat-bins-.*.fa/*.metabat-bins-.*.fa.result.SCG_table.txt')

	#get header from first file and write to out file
	f = open(SCG_dir[0], "rU")
	line = f.readline()
	args.output_file.write(line)
	f.close()

	for item in SCG_dir:
		print item
		f = open(item, "rU")
		line = f.readline()
		line = f.readline()
 		args.output_file.write(line)
		f.close()


	
	print "Script finished..."
	sys.exit(0)

#create an argument parser object
#description will be printed when help is used
parser = argparse.ArgumentParser(description='A script to run AMPHORA2')

#add the available arguments -h and --help are added by default
#if the input file does not exist then program will exit
#if output file does not exit it will be created
# args.input is the input file Note: cant write to this file because read only
# args.output is the output file
# args.m is the minimum seq length
#parser.add_argument('-i', '--input_file', type=argparse.FileType('r'), help='input file',required=True)
parser.add_argument('-o', '--output_file', type=argparse.FileType('w'), help='output file', required=True)
parser.add_argument('-e', '--evalue', help='evalue', required=True)
parser.add_argument('-t', '--type', help='type: Archaea, Bacteria, Mixed', required=True)
parser.add_argument('-p', '--processors', type=int, help='number of processors', required=True)
parser.add_argument('-s', '--seq_type', help='seq_type: DNA or PROT sequence', required=True)

#get the args
args = parser.parse_args()

#additional argument tests
if args.type != 'Archaea':
	if args.type != 'Bacteria':
		if args.type != 'Mixed':
			print "Error: argument -t not Archaea, Bacteria, or Mixed"
			sys.exit(0)
if args.seq_type != 'DNA':
	if args.seq_type != 'PROT':
		print "Error: argument -s not DNA or PROT"
		sys.exit(0)


if args.processors <= 0:
	print "Error: argument --processors <= 0"
	sys.exit(0)


#Test print the args
#print args


#cat_SCG_results() #temp remove when done



#check if files exist from a previous time
#if os.path.isfile('*.aln') or os.path.isfile('*.mask') or os.path.isfile('*.pep') or os.path.isfile('*.phylotype'):
#	print("Files exist from previous run. Please delete or move to a different folder")
#        sys.exit(0)

bin_files = glob.glob('*.metabat-bins-.*.fa')

#check if directory already exists
for item in bin_files:
	#print item
	if os.path.exists("SCG_" + item):
		print "directory already exists = " + item
		sys.exit(1)
print "The number of Metabat files = ", len(bin_files)

#make all the directories
for item in bin_files:
	os.mkdir("SCG_" + item)

	ret_dir = os.getcwd() #get the current working dir
	os.chdir("SCG_" + item)
	
	#run SCG on the fastafile
	print "Running SCG on " + item

	cmd = "python /ORG-Data/scripts/bin/Phylogeny_Protpipe/single_copy_genes.py -i " + ret_dir + "/" + item 
        cmd = cmd + " -o " + item + ".result -e " + args.evalue + " -t " + args.type + " -p " + str(args.processors) + " -s " + args.seq_type
	toolbox.run_system(cmd)

	os.chdir(ret_dir)  #return back to the original directory

#sys.exit(0)

cat_SCG_results()


#print "Script finished..."
