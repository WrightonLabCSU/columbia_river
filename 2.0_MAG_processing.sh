#Genome Binning 

##In order to get bins, we need to get coverage of our assembled scaffolds. We use Bowtie for this, and "$element" is a list of sample names.

bowtie2 --fast -p "$9" -x "$element"_scaffold.fa -S "$element"_all_mappedtoall_paired.sam -1 $forward_reads -2 $reverse_reads --un unmapped_paired.fq --al mapped_paired.fq > bowtie_log

##Next we can use this to generate our bins with metabat. This takes a list of sample names and is run in the directory where you have assembled scaffolds and the reads mapped back to the assemblies. runMetaBat.sh is the same as Metabat, but in addition it runs samtools view sam to bam conversion, and then samtools sort on that bam file for sorted bam.

runMetaBat.sh --verysensitive "$element"_scaffolds.fa "$element"_all_mappedtoall_paired.sam

#AMPHORA2 was run with in-house pipeline that uses AMPHROA2. Command is below and script details are in 2.1. "MAG_ID" is the ID of each MAG.

single_copy_genes.py  -i "MAG_ID"_scaffold.fa -o "MAG_ID"_scg.txt -e 1e-3 -t Mixed -p 6 -s DNA

#Then we run GTDB-tk to get taxonomy of our MAGs.

gtdbtk classify_wf --extension fa --genome_dir /genome_dir --out_dir /gtdb_out

#Then we dereplicate genomes with dRep.

/home2/opt/bin/dRep dereplicate_wf /dereplication_bins/ -g /genone_dir --skipCheckM #AMPHORA2 was used to inspect completion so checkM is skipped here.

#We can then annotate genomes using DRAM.

source /opt/Miniconda2/miniconda2/bin/activate DRAM
DRAM.py annotate -i '*.fa' -o DRAM_out --min_contig_size 2500 --threads 35 &> log.txt

#For Protein Trees, we used an in-house script called ProtPipeliner that is available here: https://github.com/TheWrightonLab/Protpipeliner.

#first align using MUSCLE
muscle -in seqs_for_tree.faa -out seqs_for_tree_aligned.faa
#Then run the pipeliner.
protpipeliner.py -i full_tree_seqs_aligned_trimmed.fasta -b 100 -m none -a T -t 24 #run Protpipeliner after manual trimming of genes to correct sizes.

#Mapping of metagenomic reads to MAGs to calculate relative abundance. -f is concatenated genomes file. -r is list of read locations. -n is names of each read set. Script map_reads_fasta_abundance.py can be found in 2.2.

python /ORG-Data/scripts/map_reads_fasta_abundance.py -f 102_genomes_concatenated.fasta -r reads_list.txt -n reads_names.txt -m 7 -u T -p 40 -t doutput_mapping_95_75_102_GenomesAbund_table.txt --min_coverage 3 --min_percent_contig_coverage 75
