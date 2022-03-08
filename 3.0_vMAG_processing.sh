#!/bin/bash

##Identify viruses from assembled metaG.
#VirSorter 1.0.3 was run on CyVerse Discovery Environment. Briefly, metagenomes were screened for DNA viral sequences using VirSorter v1.0.3 with the ViromeDB database option (65), retaining viral contigs ranked 1, 2, 4 or 5 where category 1-2 indicate high confidence predicted lytic viruses and 4-5 indicate high-confidence prophage sequences from VirSorter output.

#Clustering of VirSorter identified viral contigs
Cluster_genomes_5.1.pl -f unclustered_viruses_ge10kb.fasta -c 85 -i 95

#Annotation of vMAGs with DRAM-V
#First we remove "bad" characters from fasta of concatenated, clustered viruses ge10kb:
DRAM-v.py remove_bad_characters -i clustered_vMAGs_95-85.fasta -o clustered_vMAGs_95-85_nobadchars.fasta

#Then we remove bad characters from the concatenated VirSorter tab output files:
DRAM-v.py remove_bad_characters -v viral-affi-contigs-for-dramv.tab -o viral-affi-contigs-for-dramv_nobadchars.tab

#Finally we annotate and distill.
DRAM-v.py annotate -i clustered_vMAGs_95-85_nobadchars.fasta -v viral-affi-contigs-for-dramv_nobadchars.tab -o /DRAMv_output --min_contig_size 0 --threads 30 --min_contig_size 0 &> log.txt

DRAM-v.py distill -i DRAMv_output/annotations.tsv -o /DRAMv_distill

#Phyre2 analyses were done on the Phyre2 online servers (http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index) by uploading an amino acid sequence containing AMGs of interest.

##Viral Taxonomy
#vContact2 was run on CyVerse Discovery Environment. This was done on software version 0.9.8 with default settings. Files that were used to run this were deposited on Zenodo under DOI: 10.5281/zenodo.6310084

##Virus - Host Identification
#After both hosts and viruses were acquired, we then run VirHostMatcher on viral and microbial genomes to determine host predictions:
python /opt/VirHostMatcher/vhm.py -v virus_genomes_dir -b host_genomes_dir -o output_virhostmatcher_results

#Read mapping for relative abundances of vMAGs. Script map_reads_fasta_abundance.py is available in 2.2 above and was the same used for MAGs with different parameters:
map_reads_fasta_abundance.py -f clustered_vMAGs_95-85_nobadchars.fasta -r reads_list.txt -n reads_names.txt -m 15 -u T -p 90 -t output_abundances.txt --min_coverage 2 --min_percent_contig_coverage 75
