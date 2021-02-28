from os import listdir
from os.path import exists
from sys import argv, exit
from collections import Counter
from math import log
import numpy as np
import shelve
import pysam
import os

# Function to calculate Gini Coefficient
def gini(array):
    array = np.sort(array) #values must be sorted
    index = np.arange(1,array.shape[0]+1) #index per array element
    n = array.shape[0]#number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) 


def gbd_gini(sortedBam, offset, tree_dict):
    gini_cds = {}
    reads_genome = open("gbd_gini.csv","a")
            
    # Read Lengths Index
    lengthofreads = [25,26,27,28,29,30,31,32,33,34]
    
    # Open File      		
    filenamelist = sortedBam.split("/")
    filename = filenamelist[-1]
    file_id = filename[:-4]
    infile = pysam.Samfile(sortedBam, 'rb')
		    
    # Initialise Variables
    cds_reads = 0
    utr5_reads = 0
    utr3_reads = 0
    pos_list_cds = []
				
    #For every read in the sam/bam file
    for read in infile.fetch():	
	if not read.is_unmapped:
	    #read is mapped
	    readlen = read.query_alignment_length
	    if readlen <= 24 or readlen >=35:
		continue
	    chrom = read.reference_name
	    if not read.is_reverse:
		cpositions = read.positions
		cpositions.sort()
		pos = cpositions[offset[lengthofreads.index(readlen)]]
		for key in tree_dict:
		    if key[-1] == chrom[-1]:
			if tree_dict[key]["CDS+"].overlaps(pos) == True:
			    cds_reads += 1
			    pos_list_cds.append(pos)
			elif tree_dict[key]["UTR5+"].overlaps(pos) == True:
			    utr5_reads += 1
			elif tree_dict[key]["UTR3+"].overlaps(pos) == True:
			    utr3_reads += 1
			else:
			    continue  
		    else:
			pass
	    else:
		cpositions = read.positions
		cpositions.sort()
		pos = cpositions[-1-offset[lengthofreads.index(readlen)]]
		for key in tree_dict:
		    if key[-1] == chrom[-1]:
			if tree_dict[key]["CDS-"].overlaps(pos) == True:
			    cds_reads += 1
			    pos_list_cds.append(pos)
			elif tree_dict[key]["UTR5-"].overlaps(pos) == True:
			    utr5_reads += 1
			elif tree_dict[key]["UTR3-"].overlaps(pos) == True:
			    utr3_reads += 1
			else:
			    continue  
		    else:
			pass
		    
    infile.close()
    
    # Get CDS count list
    CDS_count_dict = Counter(pos_list_cds)
    count_list_cds = list(CDS_count_dict.values())
    
    # Get GBD proportions
    ###
    prop_5UTR = round(float(utr5_reads) / sum([float(utr5_reads),float(cds_reads),float(utr3_reads)]),3)
    prop_CDS = round(float(cds_reads) / sum([float(utr5_reads),float(cds_reads),float(utr3_reads)]),3)
    prop_3UTR = round(float(utr3_reads) / sum([float(utr5_reads),float(cds_reads),float(utr3_reads)]),3)
    GBD = round(log((prop_5UTR + prop_CDS)/prop_3UTR,2),2)
    if GBD < 0:
	STD_GBD = 0.00
    else:
	STD_GBD = round((GBD/10),2)
    
    # Gini Coefficient: CDS Reads    
    count_list_cds[0] = float(count_list_cds[0])
    count_array_cds = np.array(count_list_cds)
    cds_gini = gini(count_array_cds)
    STD_GINI = round(cds_gini,2)
			    
    #Write to dictionary
    gini_cds[file_id] = count_list_cds
	    
    # Output file
    reads_genome.write("{},{},{},{},{},{}\n".format(file_id,prop_5UTR,prop_CDS,prop_3UTR,STD_GBD,STD_GINI))
    reads_genome.close()
    gini_shelve = shelve.open("gini.shelf")
    gini_shelve.update(gini_cds)
    gini_shelve.close()	
    return
 
 
         
            
