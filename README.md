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

2) Use custom script seqQualityFilters.py to remove sequences with wrong barcode length, sequences lacking forward or reverse primers at the right positions.
		
3) In QIIME validate the mapping file
	
	validate_mapping_file.py -m  mapping_file_forward_reads.txt -o validate_mapping_file_output
	
4) demultiplex the reads based on split_libraries QIIME script
	
	apparently no error correction for hamming barcodes are available, used :  
	
	split_libraries_fastq.py -i ./trimmed64offsetM1.fastq -o ./demultiplexedM1q29p99/ -b ./barcode64offsetM1.fastq --barcode_type 8 --max_barcode_errors 0.5 -q 29 -p .99 -m ./mapWeek2Marker5 

5) cluster the reads into OTUssplit_libraries_fastq.py

	I used uclust (its default is 97% similarity as the definition of cluster)
	
	pick_otus.py -i ./seqs.fna -o ./otus_s98uclust/ -s 0.98 -m uclust
		
6) Making an OTU table

	make_otu_table.py -i ./seqs_otus.txt -o ./otu_table.biom
	
7) Rarefying the OTU table to contain at least 4000 sequences

	single_rarefaction.py -i ./otu_table.biom -o ./otu_table_rarefied -d 4000

8) Converting OTU table from .biom format to the classic version which is human readable (optional step just for viewing the OTU table)
	
	biom convert -i ./otu_table_rarefied -o ./otu_table_rarefied.txt -b
	
9) Used custom processOtuTable.py
	
	IMPORTANT : when using uclust, remove "denovo" in the rarefied OTU table before running the code
	see step 2 for instruction on how to run your own scripts


alpha_rarefaction.py -i ./otu_table_rarefied_filtered_n10.biom -o ./alphaDiversity_n10_v2 -t ./phylo_n10.tre -e 2500 --min_rare_depth 1000 -p ./parameterFile.txt -m ./allMapsMerged_v2.txt -n 5 -a -p para
beta_diversity_through_plots.py -i ./otu_table_rarefied_filtered_n10.biom -o ./betaDiversityManyMetrics_n10 -t ./phylo_n10.tre -e 2500 -p ./parameterFileBetaDiv.txt -a -m ./allMapsMerged_v2.txt

Examples of commands run to obtain distance metrics in QIIME:

make_distance_boxplots.py -m ../allMapsMerged_v3.txt -d unweighted_unifrac_dm.txt -f CountryCode -o unifracUnweightedBoxPlots_justGlobal/ --save_raw_data

make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d bray_curtis_dm.txt -f StudyTimesSubjID -o ./braycurtisBoxPlots --save_raw_data

make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d spearman_approx_dm.txt -f StudyTimesSubjID -o spearmanBoxPlots --save_raw_data

make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d pearson_dm.txt -f StudyTimesSubjID -o pearsonBoxPlots --save_raw_data

make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d binary_jaccard_dm.txt -f StudyTimesSubjID -o binaryJaccardBoxPlots --save_raw_data

make_distance_boxplots.py -m ../allMapsMerged_v2.txt -d abund_jaccard_dm.txt -f StudyTimesSubjID -o abundanceJaccardBoxPlots --save_raw_data

