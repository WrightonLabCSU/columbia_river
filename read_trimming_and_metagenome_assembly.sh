#!/bin/bash

#Takes in a list of untrimmed reads that will be trimmed. Ran on directory with all the raw read locations.
for element in $(<$1)
do
sickle pe -f "$element"_01.fastq -r "$element"_02.fastq -t sanger -o  R1_"$element"_trimmed.fastq -p R2_"$element"_trimmed.fastq -s "$element"_trimmed_out.fastq
done

#Takes in a tab-delim file where column A contains read paths for trimmed reads and column B is the sample ID.

while read a b
 do
cd $a
echo "Moved into $a"
fq2fa --merge --filter R1_"$b"_trimmed.fastq R2_"$b"_trimmed.fastq R1R2_"$b"_trimmed.fa
idba_ud -r R1R2_"$b"_trimmed.fa -o /idba_assembly --num_threads 20 #For IDBA-ud assemblies
/opt/SPAdes-3.13.0-Linux/bin/metaspades.py -t 10 -1 R1_"$element"_trimmed.fastq	-2 R2_"$element"_trimmed.fastq -m 300 -o /metaspades_assembly #For metaSPAdes assemblies on reads

done
exit
