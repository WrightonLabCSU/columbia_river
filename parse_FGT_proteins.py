#!/usr/bin/python

#PYTHON SCRIPT 
#written by: Richard Wolfe
#
#to run type: python make_5_column_feature_file.py -s <scaffold.fa> -p <prodigal protein file> ...
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
#
#
#   Example command to run :
#      
#  cols[0] = Protein Name
#  cols[1] = Protein count
#  cols[2] = comma separated list of proteins
#
#  had to import libraries to try and read Excel file
#  sudo /usr/bin/pip install xlrd openpyxl
#      Successfully installed et-xmlfile-1.0.1 jdcal-1.4.1 openpyxl-3.0.5 xlrd-1.2.0




import sys      #for exit command and maxint
import argparse #to get command line args 
                #needed to install argparse module because using python 2.6
                #and argparse comes with python 2.7
                #  sudo easy_install argparse
import os       #to run system commands
#import pandas as pd  #to read excel file


print "\n\nScript started..."

#create an argument parser object
#description will be printed when help is used
parser = argparse.ArgumentParser(description='A script to make a genbank file')

#add the available arguments -h and --help are added by default
#if the input file does not exist then program will exit
#if output file does not exit it will be created
# args.input is the input file Note: cant write to this file because read only
# args.output is the output file
# args.m is the minimum seq length
parser.add_argument('-i', '--input_file', type=argparse.FileType('rU'), help='Input tab separated file',required=True)
parser.add_argument('-o', '--output_file', type=argparse.FileType('w'), help='Name for output file',required=True)
#parser.add_argument('-s', '--sheet_name',  help='Name of Excel sheet',required=True)


#get the args
args = parser.parse_args()

#additional argument tests


infile_lines = 0
outfile_lines = 0
total_protein_count = 0   #summ all the protein counts -- should = the number of lines to out file

#read first line it is a header line
line = args.input_file.readline()
infile_lines += 1

#write the header to the out file
args.output_file.write(line)
outfile_lines+= 1

print line

#read next line
line = args.input_file.readline()

#if the file is not empty keep reading one at a time
while line:
	infile_lines += 1

       	#protein_count is in the second column
	protein_count = int(line.split("\t")[1])  #converts to int

   	total_protein_count += protein_count

	if  protein_count == 1:
		args.output_file.write(line)
		outfile_lines+= 1

	elif protein_count > 1:
		line = line.strip() # remove endlines
		#These are in the 3rd column
		protein_list = line.split("\t")[2]

		
		#this column may have qoutes
		if protein_list.startswith('"') and protein_list.endswith('"'):
   			print protein_list
			protein_list = protein_list[1:-1]  #remove 1st and last char
			print protein_list

		#the p
		proteins = protein_list.split(',')  #the list is comma separated

		#make a new line
		cols = line.split("\t")
		i = 0
		while i < len(proteins):
			#cols = line.split("\t")
			cols[0] = proteins[i]
			cols[1] = "1"
			#cols[2] = proteins[i]

			new_line = '\t'.join(cols)

			args.output_file.write(new_line + "\n")
			outfile_lines+= 1

			i += 1

	else:
		#should not happen
		print "Error protein count = 0"
		sys.exit(1) #exit on error
		


        line = args.input_file.readline()


print ""
print "Lines read from input file = ", infile_lines 
print "Lines wrote to output file = ", outfile_lines

print "Total Protein Count (should = lines wrote to out file - 1)= ", total_protein_count



print ""
print "\n\nScript finished..."
