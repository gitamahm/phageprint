#! /usr/bin/env python
import sys
from numpy import *
import matplotlib.pyplot as plt
from pylab import *
import os

# script written on 11/14/2014

# note rows are abundance values for a given OTU across all samples, and columns are abundances of all OTUs across a given sample. 

table = [[0.625,0.375,0,1.5,1.25,1.333333333,0.556666667,1,17,0,0],[2.625,5,0,0,0.625,0,0,0,3.5,0,0],[0.4175,0.25,0,0.166666667,0.875,1.166666667,1.556666667,1,0.333333333,2,2],[1.0825,1.625,3,1,14.5,2.333333333,0.89,14,1.166666667,3,2],[1.1675,0.125,0,0.166666667,0.25,5.333333333,0.223333333,0,0.333333333,1,1],[0.7925,42.25,2,0.333333333,5.75,1.333333333,0.943333333,2.5,1.5,1,0],[0.1675,0.625,0,0.5,0.125,0.5,0.39,0.5,22.66666667,1,0],[0.125,0.25,3,11,0.125,3.333333333,0.5,1,0.666666667,2,0],[0.915,0.25,6,0.166666667,1,8.333333333,2.166666667,0.5,0.166666667,2,1],[1.165,1.125,1,0.333333333,0.75,7.5,0.89,0.5,0.666666667,1,2],[0.1675,0.625,1,0.5,0.625,3,0.166666667,0,0.5,1,0],[49.79,2.875,5,2.666666667,3,2.333333333,3.333333333,4,3.5,5,3],[2.25,4,3,1.666666667,2.375,2.833333333,2.723333333,3,5.833333333,0,0],[32.54,41.5,26,32.16666667,26,39.5,1746.723333,42,41,40,49],[1.75,18.75,0,0.166666667,1.625,2.333333333,1.11,3,2.5,0,2],[5.5425,15.25,3,7.333333333,12.625,313.3333333,7.61,13.5,12.66666667,3,2],[1.375,4.625,596,3,2.625,3.5,4.5,4,1.166666667,10,6],[1.2075,0.625,2,0.333333333,0.25,0.833333333,0.5,1,0.666666667,735,1],[0.54,3,1,0.666666667,4.75,2.166666667,1.556666667,155,0.666666667,15,5],[14.7925,62.125,25,22.16666667,66,24.33333333,12.44333333,45,5.833333333,18,30],[0.9175,0,0,24.16666667,0,0.333333333,0,0,0,0,0],[0.1675,1.625,1,0.666666667,1.375,14,0.276666667,1,0.333333333,0,1],[766.625,6.875,16,10.83333333,7.875,4.666666667,4.223333333,13,4.833333333,3,0],[72.25,1587.75,38,69,131.875,77,61.5,70,91.83333333,94,53],[20,30.25,1,21.83333333,24.125,1171.833333,20.39,15.5,21.66666667,3,6],[2.4175,2.125,3,0.5,6.375,1.5,1.166666667,11,0.833333333,3,0],[34.0825,59.75,33,8.166666667,170.375,19.66666667,13.5,61.5,6.833333333,35,62],[0.29,1.875,1,7.5,1.375,0.5,0.89,9.5,0.833333333,0,1],[0.585,3.375,31,53.5,1.875,2.666666667,1.056666667,45,1.166666667,0,4],[80,2.375,2,0.333333333,2.25,2.5,3.943333333,2,3.833333333,3,5],[1.79,0.875,2,1.166666667,0.25,1.166666667,38.44333333,0,0.333333333,1,0],[0.3325,0.75,4,1,0.5,0,1.056666667,0.5,0.666666667,4,12],[1.29,0.875,16,46.5,1.5,2.333333333,2.276666667,5,0.833333333,2,1],[1.29,6.875,18,8.166666667,2,1,2.276666667,36,1.166666667,2,0],[19.0425,0.25,0,0,0.5,0.166666667,0,0,0.666666667,0,0],[0.415,3.75,0,8.5,1.625,0.833333333,20.61,4,11.33333333,0,0],[0.1675,7,1,2.333333333,0.5,4,1,4.5,0.166666667,0,0],[13.2075,1.5,0,0.833333333,2.25,0.333333333,0.61,1,1.166666667,1,0],[5.5,0.125,1,0,0.125,0.166666667,0,0,0.833333333,0,0],[6.7075,0.125,0,0,0,0,0.166666667,0,0,0,1],[3.96,15.25,28,32.66666667,3.125,9.166666667,8.39,27,5.166666667,3,6],[14.5,1.5,45,0.5,4,1.166666667,0.776666667,0,0.666666667,1,1],[4.415,5.875,19,5.833333333,65,40,5.943333333,50,8.833333333,13,45],[2.5,3.875,1,3.333333333,1.375,2.166666667,158.7233333,4,4.333333333,0,0],[2.335,2.125,1,0.333333333,0.625,1.166666667,0.723333333,0.5,8,10,1],[1.0825,0.375,0,0.333333333,0.5,13.16666667,1,0,0,0,0],[9.5425,15.375,12,7.166666667,9,24,108,15,19.83333333,94,38],[0,0.375,1,0.333333333,0,3.166666667,7.61,0,0.166666667,1,2],[15.2075,24.25,8,3.5,81.75,3.666666667,6.833333333,63.5,9.5,49,6],[1.54,2.25,13,14.16666667,3.125,9.833333333,4.39,2.5,0.666666667,4,53],[24.3325,83.875,122,105.3333333,33.875,16.33333333,49.11,1982,38.33333333,40,37],[1.0825,0.875,1,0.833333333,0.625,115,0.5,0.5,0.666666667,2,1],[3.79,1.625,1,2.5,0.625,1,5.39,1,1.166666667,0,1],[0.6675,0.875,1,0.333333333,0.5,1,0.223333333,0,6.166666667,0,0],[0.25,0.5,0,0.5,0.875,0.5,9.056666667,1,4.833333333,0,0],[1.125,1.875,2,1.666666667,1.625,0.5,118,1,1.5,0,0],[0.2075,0.875,50,72.66666667,0.25,1.333333333,1.833333333,3.5,0.5,0,1],[72.625,1.375,1,1.333333333,1.75,2.5,1.39,1.5,3.333333333,1,1],[9.3325,16.75,1,103.5,10,6,11.83333333,14.5,14.5,2,3],[89.6675,162,58,60.5,1948.25,68.83333333,58.11,98.5,78,73,105],[3.46,0,0,0.166666667,0,1,0.723333333,0,0.333333333,0,0],[8.4175,24.75,92,144.5,18.5,18,15.72333333,34,16.16666667,19,13],[1.9175,31.625,1,0.166666667,3.25,0.166666667,0.276666667,0,1,1,0],[104.2075,116.875,19,48.66666667,82.875,66.66666667,86.83333333,71.5,2070.166667,57,33],[5.585,0.5,0,1,0.25,0.5,0.443333333,0,0.5,0,1],[4.125,36.375,117,113.3333333,53.125,13.33333333,24.22333333,59.5,6.5,11,97],[14.2075,3.125,0,1.833333333,1.5,11.5,1.776666667,2.5,28.16666667,21,4],[3.4575,24.25,0,0.5,4,4.666666667,1.89,0.5,0.666666667,10,1],[20.8325,1.625,0,0.333333333,1.5,1.333333333,1.11,1,3,0,0],[4.21,4.375,7,2.166666667,2,8.166666667,4.39,5.5,13.16666667,16,29],[263.585,10,2,3,6.125,5.833333333,7.39,5,14.5,5,5],[2.2075,0.875,2,0.333333333,0.25,0.666666667,1.61,0.5,0.333333333,26,10],[0.0825,0,0,0,0,0,0.333333333,0,0,14,0],[1.415,1,0,0.5,0.375,0.666666667,1.776666667,1,37.66666667,2,2],[0.25,14.5,1,1,9.375,1.666666667,0.11,2.5,0.666666667,1,2],[0,0.375,1,0,1.125,9.5,1.776666667,0.5,0.166666667,0,1],[3.96,0.5,1,5,0.25,0.333333333,0.11,1,0.666666667,0,0],[31.2075,27.5,6,7,19.5,9.666666667,37.22333333,13.5,13,6,3],[0,0.25,0,0,1.625,0,0.11,0,0.333333333,0,0],[0,0,0,0,2.125,0.166666667,0,0,0,0,0],[8.6675,0.875,0,0.166666667,1.125,0.166666667,0.11,1,3.666666667,0,0],[2.8325,0.75,0,0,0.625,1.5,0,0,1.666666667,0,1],[0.5425,1.125,0,0.333333333,0.875,14.33333333,0.833333333,1,0.5,0,1],[10.96,6.75,3,2.666666667,77.5,3.166666667,2.39,5.5,5.5,5,5],[0.5,0.125,21,1.166666667,0.25,0.5,0,0.5,0.333333333,0,0],[0.4575,5.75,3,0.666666667,0.25,1.5,0.556666667,0,1.833333333,3,3],[10.585,107.625,2,2.666666667,20.875,14.83333333,2,3,4.666666667,2,1],[1.4575,15.375,6,6,14.25,6.333333333,6.556666667,5.5,2.166666667,0,14],[0.7925,1,0,0.333333333,0.75,44.16666667,0.89,0,0.166666667,0,0],[4.2925,67.25,20,18.66666667,38.25,34.66666667,42.94333333,48,15.33333333,31,59],[3.665,3,1,1.333333333,1.375,3.333333333,1.166666667,0,67.66666667,1,5],[1.7925,0.25,0,0,0.375,1.166666667,0.723333333,0.5,14.83333333,0,0],[4.0425,23.75,10,3.666666667,6.625,4,8.443333333,9,6.5,78,15],[7.29,6.375,19,2.833333333,11.125,187.1666667,61,12,30.5,51,87],[3.5425,1.75,2,2.166666667,1.625,0.833333333,1.89,1.5,63.16666667,3,24],[2.21,1,0,0,0,0.166666667,0.11,0.5,0,0,0],[0,1.375,6,12.5,1.875,1.333333333,1.943333333,22,0.333333333,2,1],[1.0825,0.875,1,0.5,1.125,78.5,3.556666667,0.5,1.166666667,4,1],[33.2075,0.625,1,0.333333333,1,2,0.833333333,0.5,1.666666667,0,1],[14.7925,21.875,2,11.33333333,17.375,21.16666667,9.556666667,24.5,90.5,2,15],[0.665,1.875,5,1.833333333,9.375,0.833333333,1.723333333,13,0.5,0,6],[6.75,89.25,13,8,40.5,13,6.5,26.5,3.166666667,7,12],[1270.2925,67,27,28.5,46,73,71.83333333,36.5,90,67,49],[1,0.625,11,0.166666667,0.625,2.166666667,0.833333333,0.5,0.333333333,2,6],[1.04,5.625,6,62,5,5.333333333,4.056666667,22.5,3,4,11],[3.7075,4.125,1,1,2.5,3,2.443333333,3,76,0,0],[0,0.5,7,4.333333333,0.125,2.166666667,0.166666667,0.5,0,0,0],[1.0425,6.375,3,4.666666667,1.25,1.833333333,45.44333333,1,4.5,1,6],[5.4575,0.5,0,0,0.25,1.666666667,0.89,0.5,22,5,3],[9.835,91,31,26.83333333,125.75,39.16666667,32.5,146,17.5,35,22],[0.9175,1,1,0.166666667,1.75,0.5,0.666666667,0,36.16666667,3,1],[1.125,1.75,4,1.833333333,0.875,5,2.666666667,1,57.33333333,5,1],[0.375,0.5,4,4.5,0.375,8.833333333,23.72333333,1.5,1.666666667,4,2],[0.25,0.5,0,0.166666667,0.5,0.5,47.72333333,0,0.666666667,0,0],[4.25,14.75,14,11,24.125,8,5.943333333,22.5,2.666666667,2,1],[46.4575,1,2,1,0.625,1.5,1.443333333,0,0.333333333,0,0],[0.0825,0.125,36,7.666666667,0.125,0.333333333,0,0,0,0,0],[1.9175,6.625,10,25.5,3.625,13,2.333333333,41,0.833333333,2,1],[0,0.125,1,0.166666667,0.25,31.5,0.61,0,0.333333333,1,0],[10.1675,16.125,8,11.33333333,20.5,9.666666667,12.5,11,192.6666667,10,6],[48.25,72.75,125,53.66666667,54.75,57.66666667,87.27666667,55,38.83333333,1578,128],[9.9575,0.75,693,2.166666667,0.875,1,0.776666667,0.5,0.666666667,1,2],[6.335,12.75,13,473.5,9.75,7,13.05666667,15.5,10.33333333,13,14],[3.21,4.625,3,5.166666667,3.75,3.833333333,161.1666667,3.5,4.333333333,12,14],[0.625,0.5,0,4,1.25,0.5,0.833333333,0,1.166666667,0,0],[11.5,319,11,21,15.375,14.5,19.05666667,31,15.33333333,25,9],[2.5425,0.75,4,0.5,0.875,0.5,0.723333333,0.5,0.333333333,1,2],[117.4575,10.125,2,4.5,8.5,7.333333333,7.833333333,6,10,8,9],[136.04,15.25,7,3.666666667,13.875,14.16666667,5.223333333,7,6.333333333,1,1],[36.21,4.75,1,1,5,2.5,1,3,31.16666667,1,2],[6.96,22.625,6,3.166666667,137.75,3.5,4.276666667,28.5,3.833333333,6,8],[33,3.75,10,3.166666667,4.875,19.5,6.833333333,3.5,42.5,6,81],[19.0425,2,1,1.5,1.125,0.666666667,30.61,2,14.5,0,1],[6.3325,0.125,0,0.166666667,0,1,0.166666667,0,0.5,0,0],[4.4575,1.375,0,0.333333333,0.75,8.666666667,1.223333333,1,1,1,1],[0.0825,1,2,16.66666667,1.25,0.666666667,0.333333333,1.5,0.333333333,1,0],[19.6675,27,108,15.5,102.25,79.33333333,41.94333333,95,83.5,205,201],[1.21,1.5,0,1.333333333,0.375,35.16666667,0.333333333,0,0.5,0,0],[10.125,14.5,7,3.833333333,9.5,3.833333333,68.66666667,12,18.83333333,4,7],[0.875,1.875,2,1.5,0.5,15.83333333,1.723333333,1,5.666666667,4,2],[0.5,0.75,10,22.66666667,0.5,4.833333333,6.61,2,1.166666667,0,1],[0.5425,2.125,650,24.5,2.5,2.5,3.5,5.5,1.166666667,5,6],[0.375,0.875,5,1.166666667,7.625,2.5,1.61,0.5,1,6,9],[2.71,61.875,0,1.333333333,29.875,2.333333333,2.223333333,6,3.333333333,3,0],[0.3325,0.875,7,1.166666667,0.25,3.333333333,1.056666667,1,1,18,4],[6.4575,0.375,0,0.166666667,0.375,1,0.166666667,0,1.333333333,1,0],[0.415,1.5,30,45.5,1,1.666666667,1.776666667,3,1.166666667,1,1],[0,0,0,0,0,0,5.666666667,0,0,0,0],[2.085,3.25,1,2.833333333,2.125,1.833333333,187.5566667,0.5,14.66666667,6,0],[0.0825,0,0,0.166666667,0.125,0.166666667,0.276666667,0,3.666666667,0,0],[27.79,42.625,265,37.83333333,43.125,23.5,49.22333333,48,20.33333333,115,1921],[18.7075,26.375,12,13.5,28.25,677.6666667,22.77666667,11.5,26.5,23,18],[3.375,1.875,1,0.333333333,3,0.333333333,3.333333333,0,1.5,1,1],[0.3325,0.625,1,8.166666667,0.625,1.833333333,1,1,1,0,0],[34.3325,0.5,1,0.666666667,1.625,1.666666667,1.666666667,0.5,1.333333333,5,2],[1.9575,45.125,2,2.333333333,3.25,2.833333333,2.056666667,1,1.833333333,2,2],[3.625,6.75,8,16.16666667,4.5,30.66666667,5.943333333,3.5,69.5,3,3],[0.7075,0.375,0,0,0.625,30.16666667,0.166666667,0,0.5,0,0],[4.1675,3.25,0,0.166666667,3.5,3.333333333,0.39,1,0.333333333,0,0],[5.7925,11.375,75,7.333333333,15.875,11.66666667,13.27666667,14.5,5.166666667,25,330],[13.6675,36.5,16,8.166666667,29.25,11.16666667,50.94333333,30,9.166666667,25,18],[0.25,0,3,8.166666667,0.125,0,0.11,0,0.166666667,0,0],[0.165,0,0,0.666666667,0,0,3.666666667,0,0.333333333,0,2],[13.25,38.375,141,1757.666667,24.5,28.16666667,30.05666667,29,15.16666667,30,39],[1.3325,5.75,2,2,4.5,7.166666667,2.443333333,48,2.166666667,2,1],[0.1675,0.5,1,8,0.25,0.333333333,0.666666667,6.5,0,0,0],[0.25,1.25,12,39.66666667,0.5,2.5,1.5,2,0.166666667,0,1],[19.2075,45.875,20,3.333333333,59.25,15.33333333,22.89,19.5,8.833333333,56,44]]
samples= ['H.M5.03','H.M5.06','H.M5.10','H.M5.H12','H.M5.16','H.M5.17','H.M5.22','H.M5.24','H.M5.25','H.M5.37','H.M5.39']
subjects= ['S03','S06','S10','S12','S16','S17','S22','S24','S25','S37','S39'] #this is the list of samples but shortened so every element in the samples list provides subject information. The order is perserved (important). I created this list using a simple grep command on the samples list. Need this info for partitioning the network diagram's edges based on the subject each edge originates from. 
sites = ['NA','NA','NA','NA','NA','NA','NA','NA','NA','NA','NA']  #same as subjects list. Used grep command to extract oral site information from the samples list (order is perserved). 
otus = ['165','133','131','130','137','136','135','134','139','138','166','24','25','26','27','20','21','22','23','160','28','29','0','4','8','120','121','122','123','125','126','127','128','129','69','59','132','54','57','56','51','50','53','52','179','114','88','89','111','110','113','112','82','181','80','119','118','84','85','3','177','7','108','109','102','103','100','101','106','107','104','39','38','33','32','31','30','37','36','35','34','162','86','60','61','62','63','64','65','66','178','68','176','175','174','173','172','171','170','182','183','180','2','186','184','185','6','99','98','168','169','164','90','93','167','95','161','97','96','11','10','13','12','15','17','16','19','18','117','151','150','153','152','155','154','157','156','159','158','83','48','49','47','44','45','42','43','40','87','1','5','9','146','147','145','142','143','140','141','148','77','76','74','73','72','71','79','78']
# just making copies of our otu table so that after we modify it we can still look up the old version if necessary. 

newtable=table[:]

# Previously, this was used as a filter by setting to zero the read count for OTUs with less than 10 sequences per sample, and then getting rid of any OTUs that don't show up at all in samples after this step. 

count=0
for col in range(len(newtable[0][:])): # samples
	for row in range(len(newtable)):  # otus
		if newtable[row][col] < 0.0:   # this basically negates the filter
			newtable[row][col] = 0.0
			count = count+1
print 'this is the count of OTUs each defined by less than 10 sequences within the OTU table: \n', count

#removing the OTUs that now appear in zero samples after the above filter was implemented. 

noise=[]

for row in range(len(newtable)):
	if sum(newtable[row]) == 0.0:
		noise.append(row)  #collecting the indices that correspond to OTUs that need to be filtered out of the new OTU table.		
print 'these are the OTUs that are now filtered out because they do not show up in any sample after the first filter was applied: \n', noise, '\n and the number of them:', len(noise)

bad_otus = set(noise)  #making a set out of the indices that correspond to OTUs that need to be filtered out of the new OTU table.
newtable = [ value for index, value in enumerate(newtable) if index not in bad_otus]
newotus = [value for index, value in enumerate(otus) if index not in bad_otus]
#print 'this is the new otu table: \n', newtable  # the filter is turned off right now
#print 'this is the list of otus that passed the filter: \n', newotus	# the filter is turned off right now
		

# Calculating the sum of each each row (number of sequences assigned to an OTU across all samples). rowsum is a vector holding each OTU's sequence count. 

rowsum=[]
for row in range(len(newtable)):
	rowsum.append(sum(newtable[row]))

# Calculating sum of each column (number of sequences assigned to one sample across all OTUs). allcol is a vector holding each sample's sequence count. 

allcol=[]		
for col in range(len(newtable[0][:])):
	onecol=[]
	for row in range(len(newtable)): 
		onecol.append(newtable[row][col])	
	allcol.append(sum(onecol))
	

#calculating relative abundance (relab) for every element in the OTU matrix. And, counting the number of super rare, rare and abundant OTUs.  

relab=[] 
histcount=[] 
all_sample_otuids=[]
file = open('edge_table.csv','r+') #for edge table to feed into Gephi software
file.write('Source;Target;Type;Id;Label;Weight\n')


for col in range(len(newtable[0][:])):  
	
	onecol=[]
	sample_otuids_rare=[]
	sample_otuids_common=[]
	sample_otuids_abundant=[]
	b_t=0 #below threshold (Otu occupies 0.001 or 0.1% of the total sequence population assigned to a given sample)
	rare =0
	common=0
	abundant=0
	
	for row in range(len(newtable)): 
		
		frac = float(newtable[row][col])/float(allcol[col]) # calculating relative abundance of a given OTU in a given sample
		onecol.append(frac)
		
		if   0.001 < frac < 0.01:
			rare = rare + 1
			sample_otuids_rare.append(newotus[row])
			s = samples[col] + ';' + newotus[row] + ';Directed;' + samples[col] +';'+ subjects[col]+';' + str(frac) +'\n'
			file.write(s)
		
		elif 0.01 <= frac < 0.1:
			common = common + 1
			sample_otuids_common.append(newotus[row])
			
			s = samples[col] + ';' + newotus[row] + ';Directed;' + samples[col] +';'+ subjects[col]+ ';' + str(frac) +'\n'
			file.write(s)
		
		elif 0.1 <= frac:
			abundant = abundant + 1
			sample_otuids_abundant.append(newotus[row])
			
			s = samples[col] + ';' + newotus[row] + ';Directed;' + samples[col] +';'+ subjects[col]+ ';' + str(frac) +'\n'
			file.write(s)
		
		else: 
			b_t = b_t + 1
	
		
	histcount.append([b_t, rare, common, abundant])
	all_sample_otuids.append([sample_otuids_rare,sample_otuids_common,sample_otuids_abundant])
	relab.append(onecol)

#print 'below_threshold, rare, common, and abundant OTU counts in each sample: \n', histcount 
#print 'IDs for the rare, common, and abundant OTUs in each sample: \n', all_sample_otuids

file.close()


file2= open('node_table.csv','r+') #for creating a nodes table
file2.write('Nodes;Id;Label;Attribute;Oral_site\n')

for i,v in enumerate(samples):
	file2.write(v + ';' + v + ';' + v + ';' + subjects[i] + ';' + sites[i] + '\n')
	
for ii, vv in enumerate(newotus):
	file2.write(vv + ';' + vv + ';' + vv + ';OTU;NA\n')
	
file2.close()

