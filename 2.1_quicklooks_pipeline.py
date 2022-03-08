#!/bin/bash

#  /ORG-Data/scripts/quicklooks/QUICKLOOKS_PIPELINE.sh  $1 $2 $3 $4 $5 $6 $7
#  This script replaced SCAFFOLDS_TO_QUICKLOOKS.sh
#
# script written by Richard Wolfe
#
# pipeline from a fasta file of scaffolds to quicklooks
#
#    This script makes a folder quicklooks_length_timestamp
#    and all the results are put into this folder
#    This folder is made in the current directory so you need to cd ito the directory
#    where you want this folder made and then run this script
#
#    Ex: SCAFFOLDS_TO_QUICKLOOKS.sh  ../scaffold.fa ../../R1_All_trimmed.fastq ../../R2_All_trimmed.fastq 1000 F02_w
#
#    Note: scaffold.fa file can use any id it does not need to be scaffold_0 it 
#          could be gene1...g|any|thing (sp) more data on line..g_e_n_e_1....
#          it just needs to be a DNA file with unique ids
#          Also make sure no duplicate ids in fasta file, if duplicates ids
#               then there will be 2 genes and reads = 0 and cov = 0 in results
#          Also if * in id then there will be spaces in id in results file
#
#    Note: if running on a bin then you can use the original reads
#
#    $1 = path to scaffold.fa file Ex ../scaffold.fa
#         -this could be a degapped scaffold file such as 100_percent_scaffold_degap
#    $2 = path to forward reads Ex  ../../R1_All_trimmed.fastq or NA if NO_BOWTIE
#    $3 = path to reverse reads Ex: ../../R2_All_trimmed.fastq or NA if NO_BOWTIE
#    $4 = min length of cotigs  Ex:  1000
#    $5 = Job ID  Ex: F02_w 
#       -can be anything (no spaces) just added to front of scaffold as an id
#       -ex: F02_4_scafold_67 or F_0_@_scaffold_67
#    $6 = BOWTIE or NO_BOWTIE  
#         note if NO_BOWTIE is used then $2 and $3 can say NA because these files are not used
#    $7 = email or NO_EMAIL
#    $8 = NA if BOWTIE or if NO_BOWTIE and you dont want contig stats
#         or path to Contig_coverage_score_100_Final_paired.txt (file with stats) when NO_BOWTIE 
#    $9 = number of threads
#    removed $10 = evalue default should be 0.0001
#    removed $11 = SAM or NO_SAM  default should be NO_SAM  (optional)


#make a unique directory name
prefix="quicklooks_$4_$5_"
suffix=$(date +%s)   #get seconds from 1/1/1970
dir_name=$prefix$suffix
mkdir $dir_name
cd $dir_name

#the dependicies
#usearch_path=/opt/usearch7.0.1090/usearch
#usearch_db=/ORG-Data/Database/UniRef/uniref90_no_spaces.udb

echo "Made a working directory $dir_name"

start_dir=$(pwd)
logfile=$start_dir/quicklooks_script.log

#make a variable because I dont think $7 can be seen inside the function
email_add=$7

#start a logfile
echo "Quicklooks script started" > $logfile
echo "Quicklooks script started and a log file was made at $logfile"
echo "logfile location = $logfile" >> $logfile
echo " " >> $logfile

echo "Arg 1 scaffold.fa file= $1" >> $logfile
echo "Arg 2 forward reads file = $2" >> $logfile
echo "Arg 3 reverse reads file = $3" >> $logfile
echo "Arg 4 min length of contigs  = $4"  >> $logfile
echo "Arg 5 job id = $5" >> $logfile
echo "Arg 6 BOWTIE or NO_BOWTIE = $6" >> $logfile
echo "Arg 7 email or NO_EMIL = $email_add" >> $logfile
echo "Arg 8 Contig_coverage_score_100_Final_paired.txt file = $8" >> $logfile

echo "Arg 9 number of threads  = $9" >> $logfile
#echo "Arg 10 evalue for usearch = ${10}" >> $logfile
#echo "Arg 11 SAM or NO_SAM  evalue for usearch = ${11}" >> $logfile

email_exit()
{
#this is a function that will email the log and then exit
if [ "$email_add" == "NO_EMAIL" ]
then
   echo "Script finished and No email was sent" >> $logfile
   echo "Script finished and No email was sent"
else
   #send email to user logfile will be written to the email
   echo "Script finished and sending email to $email_add" >> $logfile
   echo "Script finished and sending email to $email_add"
   mail -vvv -s "Quicklooks Job Finished" $email_add < $logfile
fi

exit 0

}



echo ' '
echo 'QUICKLOOKS_PIPELINE.sh script'
echo ' '

#check to make sure at least 10 variables
#check if 9th variable is empty
if [ -z ${9} ]
then 
   echo "You did not provide 9 attributes"
   echo "You did not provide 9 attributes" >> $logfile
   #dont send email because no email address
   exit 1
fi

#set default e_value
#e_value="0.0001"
#if [  -z ${10} ]  #if $10 is not empty
#then 
#   echo "Default evalue used"
#else
#   e_value="${10}"
#fi

#echo "evalue = $e_value"


#because we cd into another directory we need to add ../ to the files
#unless they are the full path
scaffold_file=$1
forward_reads=$2
reverse_reads=$3
contig_cov_file=Contig_coverage_score_100_Final_paired.txt

#check if files are full path
if [[ $scaffold_file != /* ]]
then
   scaffold_file=../$scaffold_file
fi

if [[ $forward_reads != /* ]]
then
   forward_reads=../$forward_reads
fi

if [[ $reverse_reads != /* ]]
then
   reverse_reads=../$reverse_reads
fi

if [[ "${6}" == "NO_BOWTIE" ]] &&  [[ "${8}" != "NA" ]]
then
   contig_cov_file=$8

   if [[ $contig_cov_file != /* ]]
   then
      contig_cov_file=../$contig_cov_file
   fi

   if [ ! -f "$contig_cov_file" ]
   then
       echo "Error .. file $8 does not exist"
       echo "Error .. file $8 does not exist" >> $logfile
       email_exit
   fi

fi


#check if files exists
if [ ! -f "$scaffold_file" ]
then
    echo "Error .. file $1 does not exist"
    echo "Error .. file $1 does not exist" >> $logfile
    email_exit
fi

if [ "${6}" == "BOWTIE" ]
then
   if [ ! -f "$forward_reads" ]
   then
       echo "Error .. file $2 does not exist"
       echo "Error .. file $2 does not exist" >> $logfile
       email_exit
   fi

   if [ ! -f "$reverse_reads" ]
   then
       echo "Error .. file $3 does not exist"
       echo "Error .. file $3 does not exist" >> $logfile
       email_exit
   fi
fi


#check if  arg $6 is correct
#check the arg
if [ "${6}" == "NO_BOWTIE" ]
then
   echo "Bowtie will not be performed"
elif [ "${6}" == "BOWTIE" ]
then
   echo "Bowtie will be performed"
   echo "Bowtie will be performed" >> $logfile
else
   echo "Optional variable 6 must be NO_BOWTIE or BOWTIE"
   echo "Optional variable 6 must be NO_BOWTIE or BOWTIE" >> $logfile
   email_exit 
fi

#check if  arg $11 is correct
#check the arg
#if [ "${11}" == "NO_SAM" ]
#then
#   echo "sam file will be deleted"
#elif [ "${11}" == "SAM" ]
#then
#   echo "sam file  will not be deleted"
#   echo "sam file  will not be deleted" >> $logfile
#else
#   echo "sam file will be deleted"
#   echo "sam file will be deleted" >> $logfile
   #email_exit
#fi


echo 'Getting stats from scaffold.fa file'
contig_stats.pl -i $scaffold_file -o scaffold.fa

#if arg 6 is BOWTIE then run bowtie
if [ ${6} == "BOWTIE" ]
then
  #start in another process
  (
  #bowtie2 is used to map the reads back to the scaffolds giving us
  #   number of reads, scaffold length, GC content and coverage

  echo 'Building bowtie2 index'
  echo ' '
  bowtie2-build $scaffold_file scaffold_fa

  echo 'Running bowtie2 on the scaffolds'
  echo ' '
  #  -p is number of processors this was originally 33
  #
  bowtie2 --fast -p "$9" -x scaffold_fa -S All_mappedtoall_paired.sam -1 $forward_reads -2 $reverse_reads --un unmapped_paired.fq --al mapped_paired.fq > bowtie_log
  #Note output such as percent reads mapped is put into the file bowtie_log

  echo 'Extracting data from sam file'
  echo ' '
  
  #run script to extract info from sam file because above may include reads that do not map
  python /ORG-Data/scripts/make_contig_cov_file.py -s $scaffold_file -o $contig_cov_file -b All_mappedtoall_paired.sam

  #remove the 6 bowtie index files
  rm -f scaffold_fa.1.bt2
  rm -f scaffold_fa.2.bt2
  rm -f scaffold_fa.3.bt2
  rm -f scaffold_fa.4.bt2
  rm -f scaffold_fa.rev.1.bt2
  rm -f scaffold_fa.rev.2.bt2

  #remove files created by this script
  #rm -f R1R2_ALL_contig_reads_paired.txt
  #rm -f R1R2_ALL_contig_reads_paired.txt_NEW
  #rm -f scaffold_info_genomes.txt
  #rm -f Contig_coverage_merge_100_paired.txt

  ) &

fi
#end if not a 6th argument

#start another process
(
echo "Pulling sequences of length $4 from $scaffold_file"
echo ' '
pullseq.py -i $scaffold_file -m $4 -o contigs_$4.fa

echo 'Running prodigal to select the genes from the pulled scaffolds'
echo ' '
prodigal -i contigs_$4.fa -o contigs_$4.genes -a contigs_$4.genes.faa -d contigs_$4.genes.fna -p meta -m


#echo 'Running usearch on the prodigal protein sequences'
#echo ' '
#$usearch_path -ublast contigs_$4.genes.faa -db $usearch_db -maxhits 1 -evalue "$e_value" -blast6out ublast_output_2.b6 -threads "$9" 
) &

#wait for both processes to complete
wait


#if the file Contig_coverage_score_100_Final_paired.txt does not exist then make an empty file
#This file would not exist if running without Bowtie or there was an error with bowtie
if [[ "${6}" == "NO_BOWTIE" ]] && [[ "${8}" == "NA" ]]
then
  echo "" > Contig_coverage_score_100_Final_paired.txt
  echo "Made an empty file Contig_coverage_score_100_Final_paired.txt"
fi

#if error in bowtie
if [[ ! -f "Contig_coverage_score_100_Final_paired.txt" ]] && [[ "${6}" == "BOWTIE" ]]
then
  echo "" > Contig_coverage_score_100_Final_paired.txt
  echo "Made an empty file Contig_coverage_score_100_Final_paired.txt"
fi


#exit

#comment out for new quicklooks
#check if files exist
#if [ ! -f "ublast_output_2.b6" ]
#then
#   echo "Error .. file ublast_output_2.b6 does not exist"
#   exit 1
#fi
#
##Step -2 - Merging the Coverage and Annotations files together
#perl_quicklook1.pl ublast_output_2.b6 $contig_cov_file > 1.txt
#end comment out

#make and empty file
#echo "" > 1.txt
touch 1.txt

#add the genes that were not annotated
python /ORG-Data/scripts/quicklooks/add_missing_annotations.py -i 1.txt -g contigs_$4.genes.faa -o 1.txt.NEW -s $contig_cov_file


#sort col1 version nums in text and then col 2 num
sort -k 1,1V -k 2,2n < 1.txt.NEW > Quicklooks_$5.txt


##Step -3 Replacing all the *'s in file
#sed 's/\*/_/g' 1.txt.NEW >  1_2.txt
sed -i  's/\*/ /g' Quicklooks_$5.txt

#Add prefix to begin of each line
sed -i "s/^/$5_/g" Quicklooks_$5.txt

#insert the header line
#sed -i '1i contig\tgene\tReads\tlength\tGC\tcoverage\tsbjct_id\tProtein_name\tTaxon\tIdentity\tevalue\tbit_score' Quicklooks_$5.txt
sed -i '1i contig\tgene\tReads\tlength\tGC\tcoverage' Quicklooks_$5.txt
echo "created file Quicklooks_$5.txt"
echo "Analysis Completed Successfully!"

#commented out and removed sam file below
#if [ -z ${11} ]
#then 
#   rm -f All_mappedtoall_paired.sam
#elif [ "${11}" == "SAM" ]
#then
#   echo "sam file not deleted"
#else
#   rm -f All_mappedtoall_paired.sam
#fi

#remove all unneeded files
rm -f All_mappedtoall_paired.sam
rm -f 1.txt
rm -f 1.txt.NEW
rm -f bowtie_log
rm -f mapped_paired.fq
rm -f contigs_5000.fa
rm -f contigs_5000.genes    
rm -f contigs_5000.genes.faa
rm -f contigs_5000.genes.fna
rm -f ublast_output_2.b6
rm -f unmapped_paired.fq
rm -f Contig_coverage_score_100_Final_paired.txt

echo "Script completed"
echo " "

echo "Script completed" >> $logfile
email_exit

