#! /usr/bin/env python
"""Author: Gita Mahmoudabadi
Phillips Lab, Caltech
11/18/2014
"""
import sys
from numpy import *
import matplotlib.pyplot as plt
from pylab import *
import os
import scipy.stats
import re 


######--------------Reading data from an input OTU table----------------------------------


# Reading data from an OTU table. The information of interest are three lists: list of OUT ids (otus), sample ids (samples), and absolute abundance matrix (called table here)

with open('otu_table_rarefied.txt','r') as f:  #counting the number of lines in the file (used to define number of loops in a later for loop)
    numlines = sum(1 for _ in f)

f= open('otu_table_rarefied.txt','r')
table = []
otus =[]



for i in range(2,numlines):  #starting on the third line where the numbers start
	a= str(f.readline())
	b = a.split()
	otus.append((b[0]))  #the first element hold OTU id 

	vec = []
	for j in range(1,len(b)):  #the second element onward holds absolute abundance values for a given otu across all samples. 
		vec.append(float(b[j]))
	table.append(vec)  
f.close()

newtable=table[:] #just creating a copy of our otu table 


######--------------Removing Noise-based OTUs---------------------------------------------

noise=[]
for lsto in range(len(newtable)):     #0:188 for HB1  (otus)
	count=0
	for lsti in range(len(newtable[0])):    #0:16 for HB1 (samples)	
		if newtable[lsto][lsti] <= 4.0:   ############################## Here you can define various thresholds
			count = count+1
	if count == len(newtable[0]):
		noise.append(lsto)
		

# removing the OTUs with lower than 10 sequence counts, across all samples

bad_otus = set(noise)  # making a set out of the indices that correspond to OTUs that need to be filtered out of the new OTU table.
newtable = [ value for index, value in enumerate(newtable) if index not in bad_otus]
newotus = [value for index, value in enumerate(otus) if index not in bad_otus]


samplst = [float(ii) for ii in samples] # converting the sample ids back to numbers so I can sort them

indices = sorted(range(len(samplst)), key=lambda jj:samplst[jj]) # this outputs the original indices of each element in samplist for when it's sorted. 

ord_ids =[]
for d in indices:
	ord_ids.append(str(samplst[d]))

p= len(newotus)
pp = len(indices)


table1= array(newtable[:]) # making an array of the table's copy 
table2= zeros((p,pp))

for ii in range(len(indices)):
	table2[:,ii] = table1[:,indices[ii]] #changing the columns of table2 to table1 columns but according to a particular order defined by indices
	
	
# adding the otu col to the table of numbers so I can just write them to the file more easily
otulst= [int(x) for x in newotus]
otuarray = array(otulst)
otucol = otuarray.reshape(len(newotus),1)
table3 = hstack((otucol, table2))

file0 = open('filtered_otu_table.txt','a')

head ='OTUs\t'
for id in range(len(ord_ids)):
	head = head + ord_ids[id] +'\t'

np.savetxt('filtered_otu_table.txt', table3, header= head, delimiter='\t', fmt = '%1.2f')
file0.close()



with open('filtered_otu_table.txt','r') as g:
    numlines = sum(1 for _ in g)

g= open('filtered_otu_table.txt','r')
table_ra = []
otus_ra =[]

sampleline = g.readline()  #this line holds sample information but starting at the second element after it's split
s_ra = sampleline.split()

samples = s_ra[2:]

for i in range(1,numlines):  #starting on the second line where the numbers start
	a_ra= str(g.readline())
	b_ra = a_ra.split()
	otus_ra.append((b_ra[0]))  #the first element hold OTU id 

	vec_ra = []
	for j in range(1,len(b_ra)):  #the second element onward holds absolute abundance values for a given otu across all samples. 
		vec_ra.append(float(b_ra[j]))
	table_ra.append(vec_ra)  


######--------------Calculating Relative OTU abundance------------------------------------

# Calculating the sum of each each row (number of sequences assigned to an OTU across all samples). rowsum is a vector holding each OTU's sequence count. 

rowsum=[]
for row in range(len(table_ra)):
	rowsum.append(sum(table_ra[row]))

# Calculating sum of each column (number of sequences assigned to one sample across all OTUs). allcol is a vector holding each sample's sequence count. 

allcol=[]		
for col in range(len(table_ra[0][:])):
	onecol=[]
	for row in range(len(table_ra)): 
		onecol.append(table_ra[row][col])	
	allcol.append(sum(onecol))
	
	
# calculating relative abundance (relab) for every element in the OTU matrix. And, counting the number of super rare, rare and abundant OTUs.  

relab=[] 
histcount=[] 
all_sample_otuids=[]

for col in range(len(table_ra[0][:])):  #this is 0:16 for HB1
	
	onecol=[]
	sample_otuids_absent=[]
	sample_otuids_rare=[]
	sample_otuids_common=[]
	sample_otuids_abundant=[]
	absent=0 #otus that have a relative abundance of zero. 
	b_t=0 #below threshold (Otu occupies 0.001 or 0.1% of the total sequence population assigned to a given sample)
	rare =0
	common=0
	abundant=0
	
	
	for row in range(len(table_ra)): #0:77 for HB1 
		frac = float(table_ra[row][col])/allcol[col] # calculating relative abundance of a given OTU in a given sample
		onecol.append(frac)
		
		if   0.001 <= frac < 0.01:
			rare = rare + 1
			sample_otuids_rare.append(newotus[row])
			
		elif 0.01 <= frac < 0.1:
			common = common + 1
			sample_otuids_common.append(newotus[row])
		
		
		elif 0.1 <= frac:
			abundant = abundant + 1
			sample_otuids_abundant.append(newotus[row])	
		
		elif 0 < frac < 0.001: 
			b_t = b_t + 1
			
		elif frac == 0:
			absent = absent + 1
			sample_otuids_absent.append(newotus[row])
	
		
	histcount.append([absent, b_t, rare, common, abundant])
	all_sample_otuids.append([sample_otuids_rare,sample_otuids_common,sample_otuids_abundant])
	relab.append(onecol)



print('absent, below threshold, rare, common, and abundant OTU counts in each sample: \n', histcount) 

#print 'IDs for the rare, common, and abundant OTUs in each sample: \n', all_sample_otuids


######-----------------Writing Relative Abundance matrix tab del formatted file-----------

# Writing a file containing OTU IDs (beginning of every line, starting at the second line) and subject IDs (top line). This file contains the relative abundance values of surviving OTUs in each sample. 

file1 = open('rel_abundance.txt','a')
file1.write('OTUs\t')
for sample in range(len(samples)):
	file1.write(samples[sample] + '\t')
file1.write('\n')

for otu in range(len(newotus)):
	file1.write(newotus[otu]+'\t')
	for jj in range(len(relab[:])):
		file1.write(str(relab[jj][otu]) + '\t')
	file1.write('\n')
file1.close()

######-----------------Calculating pairwise Pearson Correlations--------------------------

# just creating an empty square matrix (n by n, where n is the number of samples for this particular marker). Making empty lists of lists isn't so easy in Python. 

n=len(samples)
lst = [None]*n
r_mat=[]
#r2_mat=[]
for y in range(n):
	r_mat.append(lst[:])  #Note here, I am only appending the COPY of lst, not lst itself because when I did that, in the later section of the code, the indices were linked. r_mat will later be a matrix holding Pearson correlation coefficient (r) values. 
	#r2_mat.append(lst)  # r^2 values


for i in range(n):
	for j in range(i,n):
		r_mat[i][j]=scipy.stats.pearsonr(relab[i], relab[j])[0]
		r_mat[j][i]=r_mat[i][j]
		#r2_mat[i][j]=(scipy.stats.pearsonr(relab[i], relab[j])[0])**2

######-----------------Writing r-matrix file ---------------------------------------------

#Writing a file that will contain Pearson correlation values.

file2 = open('r_matrix.txt','a')
file2.write('\t')
for sample in range(len(samples)):
	file2.write(samples[sample] + '\t')
file2.write('\n')

for sample in range(len(samples)):
	file2.write(samples[sample]+'\t')
	for kk in range(len(r_mat[:])):
		file2.write(str(r_mat[kk][sample]) + '\t')
	file2.write('\n')
file2.close()



######-----------------Plotting phageprints-----------------------------------------------

for d in range(n):  #which is 0:16 
	plt.subplot(n,1,d)
	plt.bar(range(n),relab[d], color="#0f87cf", alpha=.6, width=.5 )
	plt.title(sample[d], fontsize=9)
	plt.axis([0, len(newotus), 0, 0.0005])
	plt.yticks(fontsize = 5)
	plt.xticks([])
	#plt.locator_params(axis = 'x', nbins=2)
	plt.locator_params(axis = 'y', nbins = 1)

	
	#plt.xlabel('OTUs')
	#plt.ylabel('Relative OTU abundance across samples')
	
plt.show()
