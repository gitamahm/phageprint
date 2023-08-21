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

Lastly, we provide examples of the QIIME scripts used in the document called **`examples_of_QIIME_commands`**
