# Human Phageprints
#### Individuals can be identified by a small number of globally distributed, temporally persistent resident phage families
<img src="https://github.com/gitamahm/phageprint/blob/main/phageprint.png" width="100%" height="100%">

Please see the following jupyter notebooks where much of the time series and machine learning analysis was performed:
- **`phageprint_machine_learning_modeling.ipynb`**
    - the purpose of this notebook is to analyze the HB1 time series dataset.
- **`phageprint_machine_learning_modeling_HA.ipynb`**
    - this notebook is very similar to the one above but is for the time series analysis of the HA phage family.

In addition to these notebooks, there are three custom python scripts that we wrote to complement the QIIME package:

- **`seqQualityFilters.py`**
    - performs sequence quality control
- **`processOtuTable.py`**
    - performs additional quality control and denoising of the OTU table
- **`createNetwork.py`**
    - creates a network where OTUs are connected to the samples that they are found in with weighted edges based on OTU relative abundance

1) join the forward and reverse reads using fastqJoin from qiime 

	join_paired_ends.py -f ./GM01-SS-GMAC_S1_L001_R1_001.fastq.gz  -r ./GM01-SS-GMAC_S1_L001_R2_001.fastq -p 0 -o ./fastqJoinedP0

2) Use your own script quality_filter.py in the ~/ngs folder to get rid of sequences with wrong barcode length, sequences lacking forward or reverse primers at the right positions.
	
	this needs to be run from normal terminal (rather than the qiime shell). Be in the same folder as the one containing the joined reads though. And run the quality control code using its absolute address. 
	here's the absolute path of this script
	/Users/octatig88/Desktop/all/phd/rsch/bioinfo/biopython-1.64/ngs/seqQualityFilters.py
	
3) In Macqiime (type macqiime in terminal): Validate the mapping file
	
	validate_mapping_file.py -m  mapping_file_forward_reads.txt -o validate_mapping_file_output
	
	Primer 1's I's need to be converted to Ns. 

4) demultiplex the reads based on split_libraries qiime script

	split_libraries_fastq.py -i ./trimmed64offsetM5.fastq -o ./demultiplexed/ -m ../../mapWeek2Marker5 -b ./barcode64offsetM5.fastq --barcode_type hamming_8 --max_barcode_errors 0.5 -q 29 -p .99
	
	apparently no error correction for hamming barcodes are available, used :  
	
	split_libraries_fastq.py -i ./trimmed64offsetM1.fastq -o ./demultiplexedM1q29p99/ -b ./barcode64offsetM1.fastq --barcode_type 8 --max_barcode_errors 0.5 -q 29 -p .99 -m ./mapWeek2Marker5 

5) cluster the reads into OTUssplit_libraries_fastq.py

	used to use the following but it's no longer working for reasons that I can't figure out: pick_otus.py -i ./seqs.fna -o ./otus_s97usearch/ -m usearch -s 0.97 -g 50  --word_length 64  --suppress_reference_chimera_detection
	so I used uclust instead (its default is 97% similarity as the definition of cluster)
	
	pick_otus.py -i ./seqs.fna -o ./otus_s98uclust/ -s 0.98 -m uclust
	
	then cd  otus_s97uclust/
	
6) Making an OTU table

	make_otu_table.py -i ./seqs_otus.txt -o ./otu_table.biom
	
7) Rarefying the OTU table to contain at least 4000 sequences

	single_rarefaction.py -i ./otu_table.biom -o ./otu_table_rarefied -d 4000

8) Converting OTU table from .biom format to the classic version which is human readable (optional step just for viewing the OTU table)
	
	biom convert -i ./otu_table_rarefied -o ./otu_table_rarefied.txt -b
	
	***qiime 2 command is 
	biom convert -i ./otu_table_rarefied -o ./otu_table_rarefied.txt --to-tsv
	
	
if you don't want to use your own code, you can run the following to filter out low-abundant OTUS. Didn't use my own code in step 9 last run. 

	filter_otus_from_otu_table.py -i ./otu_table_rarefied.biom -n 10 -o ./otu_table_rarefied_filtered_n10.biom
	
	this will filter any OTUs with less than 10 sequences associated with it across the entire table. 
	
9) Use your own code processOtuTable.py
	
	IMPORTANT : when using uclust, get rid of "denovo" in the rarefied OTU table before running the code
	see step 2 for instruction on how to run your own scripts
	
	here's the absolute path of this script
	/Users/octatig88/Desktop/all/phd/rsch/bioinfo/biopython-1.64/ngs/processOtuTable.py


Used the following to get filter out the seq.fna file to only contain sequences that belong to OTUs that passed a filter. So the seq_otus_filtered_n10.txt had to be made using a code I wrote (in ipython) to get rid of OTUs that hadn't passed the filter. 
filter_fasta.py -f ../seqs.fna -o ./seqs_fitered_n10.fna -m ./seqs_otus_filted_n10.txt 	
	
To do all demultiplexing for different sequencing runs

split_libraries_fastq.py -i AC2trimmed64offsetM1.fastq,AC3trimmed64offsetM1.fastq -b barcode64offsetM1AC2.fastq,barcode64offsetM1AC3.fastq -m mapAC2.txt,mapAC3.txt -o ./demultiplexed --barcode_type 8 --max_barcode_errors 0.5 -q 29 -p .99
split_libraries_fastq.py -i AC1trimmed64offsetM1.fastq,AC2trimmed64offsetM1.fastq,AC3trimmed64offsetM1.fastq,AC4trimmed64offsetM1.fastq,AC5trimmed64offsetM1.fastq,AC6trimmed64offsetM1.fastq -b barcode64offsetM1AC1.fastq,barcode64offsetM1AC2.fastq,barcode64offsetM1AC3.fastq,barcode64offsetM1AC4.fastq,barcode64offsetM1AC5.fastq,barcode64offsetM1AC6.fastq -m mapAC1.txt,mapAC2.txt,mapAC3.txt,mapAC4.txt,mapAC5.txt,mapAC6.txt -o ./demultiplexedM1q29p99 --barcode_type 8 --max_barcode_errors 0.5 -q 29 -p .99 && cd ./demultiplexedM1q29p99 && pick_otus.py -i ./seqs.fna -o ./otus_s97uclust/ 	
split_libraries_fastq.py -i AC1trimmed64offsetM1.fastq,AC2trimmed64offsetM1.fastq,AC3trimmed64offsetM1.fastq,AC4trimmed64offsetM1.fastq,AC5trimmed64offsetM1.fastq -b barcode64offsetM1AC1.fastq,barcode64offsetM1AC2.fastq,barcode64offsetM1AC3.fastq,barcode64offsetM1AC4.fastq,barcode64offsetM1AC5.fastq -m mapAC1.txt,mapAC2.txt,mapAC3.txt,mapAC4.txt,mapAC5.txt -o ./demultiplexedM1q29p99 --barcode_type 8 --max_barcode_errors 0.5 -q 29 -p .99 && cd ./demultiplexedM1q29p99-AC1toAC5 && pick_otus.py -i ./seqs.fna -o ./otus_s97uclust/ 	
split_libraries_fastq.py -i AC1trimmed64offsetM5.fastq,AC2trimmed64offsetM5.fastq,AC3trimmed64offsetM5.fastq,AC4trimmed64offsetM5.fastq,AC5trimmed64offsetM5.fastq -b barcode64offsetM5AC1.fastq,barcode64offsetM5AC2.fastq,barcode64offsetM5AC3.fastq,barcode64offsetM5AC4.fastq,barcode64offsetM5AC5.fastq -m mapAC1.txt,mapAC2.txt,mapAC3.txt,mapAC4.txt,mapAC5.txt -o ./demultiplexedM5q29p99 --barcode_type 8 --max_barcode_errors 0.5 -q 29 -p .99 && cd ./demultiplexedM5q29p99-AC1toAC5 && pick_otus.py -i ./seqs.fna -o ./otus_s97uclust/ 	

To run several commands at once from OTU picking all the way to rarefied otu table. 
pick_otus.py -i ./seqs.fna -s 1 -o ./otus_s100uclust/ && cd ./otus_s100uclust/ && make_otu_table.py -i ./seqs_otus.txt -o ./otu_table.biom && single_rarefaction.py -i ./otu_table.biom -o ./otu_table_rarefied -d 4000 && biom convert -i ./otu_table_rarefied -o ./otu_table_rarefied.txt -b 


Or you can try the following: 
pick_de_novo_otus.py -i ./seqs.fna -o ./otus_s100uclust-v2/ -p ~/Gita/uclus_dereplication_params.txt --derep_fullseq

where the following is in $PWD/dereplication_params.txt:

pick_otus:similarity 1.0

pick_rep_set.py -i ./seqs_otus.txt -f ../seqs.fna -o ./rep_set.fna
align_seqs.py -i ./rep_set.fna -m muscle -o ./muscleAlignment
make_phylogeny.py -i ./muscleAlignment/rep_set_aligned.fasta -o ./phylo.tre

merge_mapping_files.py -m mapAC1.txt,mapAC2.txt,mapAC3.txt,mapAC4.txt,mapAC5.txt -o ./mapAC1to5merged.txt
had to change the human column to subject id 

beta_diversity_through_plots.py -i otu_table.biom -o ./biodiv4000/ -t ./phylo.tre -m ../../mapAC1to5merged.txt -e 4000

beta_diversity_through_plots.py -i otu_table_rarefied_filtered_n10.biom -m allMapsMerged.txt -o ./betaDiversityPlotsDepth4000 -e 4000

beta_diversity.py -i ./filtered_otu_table.biom -m abund_jaccard  -o ./betaDivJaccard -t ../phylo.tre

principal_coordinates.py -i ./abund_jaccard_filtered_otu_table.txt -o ./components

make_2d_plots.py -i ./components -m ../../../../../mapAC1to5merged.txt -o ./pca2dPlots


biom convert -i ../filtered_otu_table.txt -o ./filtered_otu_table.biom --table-type "OTU table"

alpha_rarefaction.py -i ./otu_table_rarefied_filtered_n10.biom -o ./alphaDiversity_n10_v2 -t ./phylo_n10.tre -e 2500 --min_rare_depth 1000 -p ./parameterFile.txt -m ./allMapsMerged_v2.txt -n 5 -a -p para
beta_diversity_through_plots.py -i ./otu_table_rarefied_filtered_n10.biom -o ./betaDiversityManyMetrics_n10 -t ./phylo_n10.tre -e 2500 -p ./parameterFileBetaDiv.txt -a -m ./allMapsMerged_v2.txt

had to change to allMapsMerged_v3 because the country code was including everybody, but instead I only want to analyze samples from the global study (rather than also include siblings and time, etc)
make_distance_boxplots.py -m ../allMapsMerged_v3.txt -d unweighted_unifrac_dm.txt -f CountryCode -o unifracUnweightedBoxPlots_justGlobal/ --save_raw_data


make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d bray_curtis_dm.txt -f StudyTimesSubjID -o ./braycurtisBoxPlots --save_raw_data
make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d spearman_approx_dm.txt -f StudyTimesSubjID -o spearmanBoxPlots --save_raw_data
make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d pearson_dm.txt -f StudyTimesSubjID -o pearsonBoxPlots --save_raw_data
make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d binary_jaccard_dm.txt -f StudyTimesSubjID -o binaryJaccardBoxPlots --save_raw_data
make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d abund_jaccard_dm.txt -f StudyTimesSubjID -o abundanceJaccardBoxPlots --save_raw_data
make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d unweighted_unifrac_dm.txt -f StudyTimesSubjID -o unifracUnweightedBoxPlots/ --save_raw_data
