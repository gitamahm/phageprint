#! /usr/bin/env python

"""Author: Gita Mahmoudabadi
Email: gita@caltech.edu
Phillips Lab, Caltech
11/18/2014
"""
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import re
#from Bio.Alphabet import IUPAC
import os

""" This code is made up of several parts. 

Part 1 takes in the joined reads, and writes two files: highq_reads.fastq, and
lowq_reads.fastq). The highq reads are those with higher than 29 Phred score 
for all their bases, except for the first two bases in the beginning and
end of the sequence. 

Part 2 takes in the highq.fastq generated in Part 1, and places sequences 
in 6 different bins according to their primer sequences. It also places an additional 
QC filter as it identifies any read that does not have the correct barcode length, 
and the right forward and reverse primer sequences at the right positions within the read. 
These reads are written to bad_primer_barcode_seqs.fastq.
Additionally, the orientation of all surviving reads is coordinated such that 
they all start with a barcode and their forward primer. 

Part 3 takes in the marker3.fastq file generated in step 2, and looks for Checksum 
barcodes that are incorrect. It writes sequences with the correct barcodes to one file
(correct_chksum_marker3.fastq), and all else to another file (incorrect_chksum_marker3.fastq). 
"""



#os.chdir('/Users/octatig88/Desktop/all/phd/rsch/bioinfo/biopython-1.64/ngs/miseq_data/sibling/fastqJoinedP0/')

#PART 1..................................................................................

file_j = 'fastqjoin.join.fastq' #this is the file containing all joined paired end reads
h_j=open(file_j) 
j_records=SeqIO.parse(h_j, 'fastq')                                                                                                                                          

#marker_num=input('Enter the total number of markers to analyze') 

h_highq= open('highq_reads.fastq','a')  #opening a file to keep all high quality reads in
h_lowq = open('lowq_reads.fastq','a')  #opening a file to keep all low quality reads in

c=0
ch = 0
cl = 0
#creating a match list
for j_record in j_records:   #looping through each fastq record
	
	qlist = j_record.letter_annotations["phred_quality"][2:-2] 	 #I am basically creating a list out of the Phred scores, while ignoring the first and last two bases in every sequence (they typically have low Phred scores, but that can be resolved using error detection codes for barcode sequences). 
	value = 0 
	for q in range(len(qlist)): #I am looping through the list of Phred scores for each sequence
		if qlist[q] > 25:   #if the score is greater than 29 (ascii for an error rate of 1 in 1000), value is 1, and we can continue stepping through the loop
			value = 1
		else:       #however, if there's any base with quality score less than or equal to 29, the loop is broken, the sequence is written to the file "lowq_reads.fastq"
			value = 0
			SeqIO.write(j_record, h_lowq, 'fastq')
			cl = cl+1     #keeping count of low quality reads
			break 	 	
	if value == 1:         #at the end of the interior loop, if the value is still 1, that means that the list of phred scores for that record has not contained any low quality bases, and therefore we can write the record to the highq_reads.fastq file. 
		 SeqIO.write(j_record, h_highq, 'fastq')
		 ch = ch + 1     #here we are keeping count of high quality reads

print ('this is the number of high qulaity sequences:' ,ch) 
print ('this is the number of low quality sequences:', cl)


#closing the files we opened. 
h_lowq.close() 
h_j.close()


#PART 2..................................................................................

h_highq = open ('highq_reads.fastq')
hq_records=SeqIO.parse(h_highq, 'fastq') 

#primer sequences and their barcodes 5' to 3'

srch1r='GA.CC.TA.AA.GCCAAAGAGTTCGTAA.'
srch3r='GA.GA..T.CA.GCACACCCAAATAATCA.'
srch5r='CC.ATG.TGGA.GA.AGCACTCCTTAC.'
srch6r='GA.GA.AA.GT.TT.AAGCAACTACGACTCA.'
srch7r='GG.TT.GA.GG.TCGGAAG.'
srch8r='GG.GA.TGGGG.GT.GCT.'

srchlist_r= [srch1r,srch3r,srch5r,srch6r,srch7r,srch8r]

srch1f='........CCGATCTGTC.CA.GG.GA.GA'
srch3f='....GTGCGGCAAC.AA.CA.GA.CA'
srch5f='........CGTGATGGCTG.CT.GA.TT.GA.GA'
srch6f='........CTCAATGT.TCAGG.CTGGT.TT.GA.GA'
srch7f='........TCACAATCGAGAC.CC.AA.GC'
srch8f='........CCTTTG.TTGGC.TGGTT.GA.GA'

srchlist_f = [srch1f,srch3f,srch5f,srch6f,srch7f,srch8f]


#short hand for indexing the srchlists which is made easy because both are of the same length
n= len(srchlist_f)
k= range(n)


#building two lists containing reverse complement sequences of reverse and 
#forward primers...also counting the primer lengths for later application. 

srchlist_rc_r=[]
primerlength_r=[]
srchlist_rc_f=[]
primerlength_f=[]

for l in k:
    srchlist_rc_r.append(str(Seq(srchlist_r[l]).reverse_complement()))
    srchlist_rc_f.append(str(Seq(srchlist_f[l]).reverse_complement()))
    primerlength_r.append(len(srchlist_r[l]))
    primerlength_f.append(len(srchlist_f[l]))


#opening the files that will later hold the good reads for each marker (good reads without the primer or barcode sequence)
h_1=open('trimmed64offsetM1.fastq','a')
h_3=open('trimmed64offsetM3.fastq','a')
h_5=open('trimmed64offsetM5.fastq','a')
h_6=open('trimmed64offsetM6.fastq','a')
h_7=open('trimmed64offsetM7.fastq','a')
h_8=open('trimmed64offsetM8.fastq','a')

h_list=[h_1,h_3,h_5,h_6,h_7,h_8]

#opening the files that will later hold the good reads with the barcode and primer sequences. 
h_1_wbcp=open('withBarcodeAndPrimers64offsetM1.fastq','a')
h_3_wbcp=open('withBarcodeAndPrimers64offsetM3.fastq','a')
h_5_wbcp=open('withBarcodeAndPrimers64offsetM5.fastq','a')
h_6_wbcp=open('withBarcodeAndPrimers64offsetM6.fastq','a')
h_7_wbcp=open('withBarcodeAndPrimers64offsetM7.fastq','a')
h_8_wbcp=open('WithBarcodeAndPrimers64offsetM8.fastq','a')

h_wbcp_list = [h_1_wbcp, h_3_wbcp, h_5_wbcp, h_6_wbcp, h_7_wbcp, h_8_wbcp]

#opening the files that will hold the barcodes (recall, QIIME's demultiplexing
#script wants the sequences and the barcodes separately). 

hbc_1 = open('barcode64offsetM1.fastq', 'a')
hbc_3 = open('barcode64offsetM3.fastq', 'a')
hbc_5 = open('barcode64offsetM5.fastq', 'a')
hbc_6 = open('barcode64offsetM6.fastq', 'a')
hbc_7 = open('barcode64offsetM7.fastq', 'a')
hbc_8 = open('barcode64offsetM8.fastq', 'a')

hbc_list = [hbc_1, hbc_3, hbc_5, hbc_6, hbc_7, hbc_8]

#opening files that will later hold reads with the right barcode and primers, but 
#that are either too short or too long
h_1_wl=open('wronglength64offsetM1.fastq','a')
h_3_wl=open('wronglength64offsetM3.fastq','a')
h_5_wl=open('wronglength64offsetM5.fastq','a')
h_6_wl=open('wronglength64offsetM6.fastq','a')
h_7_wl=open('wronglength64offsetM7.fastq','a')
h_8_wl=open('wronglength64offsetM8.fastq','a')

h_wl_list = [h_1_wl, h_3_wl, h_5_wl, h_6_wl, h_7_wl, h_8_wl]

#opening file to write to bad reads, containing either bad primer sequences or
#the wrong number of barcode bases. 
h_unid = open('bad_primer_barcode_seqs64offset.fastq','a')

#cumulative length of marker, barcode, the reverse and forward primers
m1 = 303
m3 = 237
m5 = 309
m6 = 391
m7 = 468
m8 = 334
length_list= [m1, m3, m5, m6, m7, m8]  
#low_len_list = [item-11 for item in length_list]
#up_len_list = [item+11 for item in length_list] 
bc_length= [8, 4, 8, 8, 8, 8]



for hq_record in hq_records:  		#looping through the high_quality records
	
	counter=0
	for i in k:	 #looping through the 6 markers/primers
	
		a= re.search(srchlist_f[i], str(hq_record.seq[0:primerlength_f[i]]))  # using grep search to see if bases 0 to forward primer length i from hq_record match forward primer i and have the right barcode length before them. Note, we're no specific about the barcode sequence at this step. 
		b= re.search(srchlist_r[i], str(hq_record.seq[-primerlength_r[i]:])) #using grep search to see if the last base to reverse primer length i from hq_record match the reverse primer i sequence.
				
		start = primerlength_f[i]  #we're defining the bases at which we want to extract the sequence when we're getting rid of the forward and reverse primers in the sequences
		end = -primerlength_r[i]
		
		if a and b:   #here we're only accepting sequences with the right forward and reverse primer sequence and position, as well as the right barcode length positioned before the forward primer. 
	
			L = len(hq_record.seq)	#calculating the length of the sequence
	
			if L == length_list[i]:   #here we're applying a length filter on the surviving sequences, using the length of marker i which we defined earlier. 
				
				SeqIO.write(hq_record, h_wbcp_list[i], "fastq-illumina")  # we're going to write this sequence.
				
				hq_bc_record= SeqRecord(hq_record.seq[0:bc_length[i]], id =hq_record.id, letter_annotations= {"phred_quality":hq_record.letter_annotations["phred_quality"][0:bc_length[i]]})	#we're also going to write a new record that holds the same sequence but without the barcode or primers. 
				SeqIO.write(hq_bc_record, hbc_list[i], "fastq-illumina")  
			
				hq_record_new = SeqRecord(hq_record.seq[start:end], id= hq_record.id, letter_annotations = {"phred_quality":hq_record.letter_annotations["phred_quality"][start:end]}) #also writing a separate record for the barcode sequence. This file and the primer-free reads will be fed into the demultiplex script in QIIME that only accepts the barcodes and the reads separately. 
				SeqIO.write(hq_record_new, h_list[i], "fastq-illumina")				
				break	#if the three records are written for this sequence, then the loop is broken, and we move onto the upper loop (i.e. moving on to the next sequence without repeating the same steps for different markers, as we already have a match). 
		
			else: 
				SeqIO.write(hq_record, h_wl_list[i], "fastq-illumina") #otherwise, this is a sequence that has the right barcode and primers, but is not the right length, so it gets written to a different file. 
				break		
	
		c = re.search(srchlist_rc_r[i],str(hq_record.seq[0:-end]))   #however, we can expect half of the sequence to be reverse complements. Here we're checking to see if the reverse complemented version of the reverse primer exists at the beginning of the sequence. 
		d = re.search(srchlist_rc_f[i],str(hq_record.seq[-start:])) #checking to see if the reverse complemented versions of the forward primer and barcode are at the end of the sequence.
	
		if c and d: #if the condition is true, the sequence is a reverse complement of a sequence that we can work with
	
			L = len(hq_record.seq) #again, we want to see if it also has the right length.
	
			if L == length_list[i]:			
				hq_record_rc = hq_record.reverse_complement(id=hq_record.id+"_rc") #for the sequences with the right length, we're going to re-orient them by reverse complementing them. 
				SeqIO.write(hq_record_rc, h_wbcp_list[i], "fastq-illumina")  #then writing this record
				
				hq_bc_record= SeqRecord(hq_record_rc.seq[0:bc_length[i]], id =hq_record_rc.id, letter_annotations= {"phred_quality":hq_record_rc.letter_annotations["phred_quality"][0:bc_length[i]]})	#also, again we're going to remove the primers and barcode and write a new record 
				SeqIO.write(hq_bc_record, hbc_list[i], "fastq-illumina")
								
				hq_record_new = SeqRecord(hq_record_rc.seq[start:end], id= hq_record_rc.id, letter_annotations = {"phred_quality":hq_record_rc.letter_annotations["phred_quality"][start:end]}) #as before, we're going to write their barcodes to a separate file
				SeqIO.write(hq_record_new, h_list[i], "fastq-illumina")
				break
			
			else: 
				hq_record_rc = hq_record.reverse_complement(id=hq_record.id+"_rc") #if the sequence doesn't have the right length, we're going to keep it and write its reverse complement to a file (for later analysis). 
				SeqIO.write(hq_record_rc, h_wl_list[i], "fastq-illumina")
				break		
		else:
			counter = counter +1 	# for marker i, if neither of the if statements is applicable, then the sequence cannot be identified with the marker i	and so we're going to keep a count of how many times this sequence hasn't been identified. If it's been 6 times (i.e. it's not matching any of the markers), then we'll write it to a file containing un-identifiable sequences.	
			if counter == n:
				SeqIO.write(hq_record, h_unid, 'fastq-illumina')


#closing the files I opened
for ii in k:
	h_list[ii].close()
	h_wl_list[ii].close()
	hbc_list[ii].close()
	h_wbcp_list[ii].close()

h_unid.close()
h_highq.close()


#PART 3..................................................................................
'''
file_m3= 'marker3.fastq'
h_m3=open(file_m3)                                                                                                                           

h_error= open('incorrect_chksum_marker3.fastq','a')
h_correct = open('correct_chksum_marker3.fastq','a')

m3_records=SeqIO.parse(h_m3, 'fastq')               

for m3_record in m3_records:
	bc_seq = m3_record[0:4]
	bc_list = []
	for c in bc_seq:
		if c == 'A':
			bc_list.append(0)
		if c == 'C':
			bc_list.append(1)
		if c == 'G':
			bc_list.append(2)
		if c == 'T':
			bc_list.append(3)

	chksum = sum(bc_list[0:3]) % 4 
	if chksum == bc_list[3]:
		SeqIO.write(m3_record, h_correct, 'fastq')
	else:
		SeqIO.write(m3_record, h_error, 'fastq')

h_m3.close()
h_error.close()
h_correct.close()
'''