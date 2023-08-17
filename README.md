# phageprint
#### Human Phageprints: Individuals can be identified by a small number of globally distributed, temporally persistent resident phage families
<img src="https://github.com/gitamahm/phageprint/blob/master/phageprint.png" width="100%" height="100%">

In the accompanying manuscript, we created an atlas of the human tissue microbiome with cell type resolution, and in this repository, we provide our pipeline called SIMBA for the identification of bacteria, viruses and fungi in human cellular transcriptomes. We used the raw fastqs from the [Tabula Sapiens](https://github.com/czbiohub/tabula-sapiens) single-cell transcriptomic dataset. 

Additionally, we have included several Jupyter notebooks used for further processing SIMBA's outputs and generating figures. Theses notebooks include:

- **`p01_databaseCreationSC_microbeDB_fungalDB.ipynb`**
    - the purpose of this notebook is to create microbeDB and other databases used in the SIMBA pipeline
- **`p02_creating_truth_datasets_for_testing_simba.ipynb`**
	- the purpose of this notebook is to create truth datasets for testing the precision and recall of SIMBA and Kraken2

