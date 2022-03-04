#bin/bash!
#Genome Binning
#AMPHORA2

#GTDB-tk
gtdbtk classify_wf --extension fa --genome_dir /genome_dir --out_dir /gtdb_out

#Dereplication of genomes
/home2/opt/bin/dRep dereplicate_wf /dereplication_bins/ -g /genone_dir --skipCheckM #AMPHORA2 was used to inspect completion so checkM is skipped here.

#Annotation of genomes using DRAM
cd /genome_dir
source /opt/Miniconda2/miniconda2/bin/activate DRAM
DRAM.py annotate -i '*.fa' -o DRAM_out --min_contig_size 2500 --threads 35 &> log.txt

#Protein Trees made using ProtPipeliner available here: https://github.com/TheWrightonLab/Protpipeliner
muscle -in seqs_for_tree.faa -out seqs_for_tree_aligned.faa #align using MUSCLE
protpipeliner.py -i full_tree_seqs_aligned_trimmed.fasta -b 100 -m none -a T -t 24 #run Protpipeliner after manual trimming of genes to correct sizes.

#Mapping of metagenomic reads to MAGs to calculate relative abundance. -f is concatenated genomes file. -r is list of read locations. -n is names of each read set. Script map_reads_fasta_abundance.py can be found here: XXX.
python /ORG-Data/scripts/map_reads_fasta_abundance.py -f 102_genomes_concatenated.fasta -r reads_list.txt -n reads_names.txt -m 7 -u T -p 40 -t doutput_mapping_95_75_102_GenomesAbund_table.txt --min_coverage 3 --min_percent_contig_coverage 75
