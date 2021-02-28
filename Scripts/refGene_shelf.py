from intervaltree import Interval, IntervalTree
import shelve
import sys

def refGeneLocations(pathToRefgene):
    infile = open(pathToRefgene,"r")

    all_cds_plus = {}
    all_cds_minus = {}
    UTR5_plus = {}
    UTR5_minus = {}
    UTR3_plus = {}
    UTR3_minus = {}

    # parse cds and urt regions from exons
    for line in infile:
        splitline = line.split("\t")
        tran = splitline[1]
        if tran.startswith("NM_"):
	    strand = splitline[3]
	    chrom = splitline[2].lower()
	    exon_number = int(splitline[8]) 
	    cds_start = int(splitline[6])
	    cds_end = int(splitline[7]) 
	    exon_starts = splitline[9].split(",")
	    exon_ends = splitline[10].split(",")
            if cds_start == cds_end:
                continue
	    if strand == '+':
	        for i in range(exon_number):
	            exon_start = int(exon_starts[i])
	            exon_end = int(exon_ends[i])
	            if exon_end > cds_start:
	               exon_end = cds_start-1
	            if exon_start>= exon_end:
	                continue                
		    if chrom in UTR5_plus:
		        UTR5_plus[chrom].append((exon_start, exon_end))
		    else:
		        UTR5_plus[chrom] = [(exon_start, exon_end)]
	        for i in range(exon_number):
	            exon_start = int(exon_starts[i])
	            exon_end = int(exon_ends[i])
	            if exon_start < cds_start:
	                exon_start = cds_start
	            if exon_end > cds_end:
	               exon_end = cds_end
	            if exon_start>= exon_end:
	                continue                
		    if chrom in all_cds_plus:
		        all_cds_plus[chrom].append((exon_start, exon_end))
		    else:
		        all_cds_plus[chrom] = [(exon_start, exon_end)]
	        for i in range(exon_number):
	            exon_start = int(exon_starts[i])
	            exon_end = int(exon_ends[i])
	            if exon_start < cds_end:
	               exon_start = cds_end+1
	            if exon_start>= exon_end:
	                continue                
		    if chrom in UTR3_plus:
		        UTR3_plus[chrom].append((exon_start, exon_end))
		    else:
		        UTR3_plus[chrom] = [(exon_start, exon_end)]
	    else:
	        for i in range(exon_number):
	            exon_start = int(exon_starts[i])
	            exon_end = int(exon_ends[i])
	            if exon_end > cds_start:
	               exon_end = cds_start-1
	            if exon_start>= exon_end:
	                continue                
		    if chrom in UTR3_minus:
		        UTR3_minus[chrom].append((exon_start, exon_end))
		    else:
		        UTR3_minus[chrom] = [(exon_start, exon_end)]
	        for i in range(exon_number):
	            exon_start = int(exon_starts[i])
	            exon_end = int(exon_ends[i])
	            if exon_start < cds_start:
	                exon_start = cds_start
	            if exon_end > cds_end:
	               exon_end = cds_end
	            if exon_start>= exon_end:
	                continue                
		    if chrom in all_cds_minus:
		        all_cds_minus[chrom].append((exon_start, exon_end))
		    else:
		        all_cds_minus[chrom] = [(exon_start, exon_end)]
	        for i in range(exon_number):
	            exon_start = int(exon_starts[i])
	            exon_end = int(exon_ends[i])
	            if exon_start < cds_end:
	               exon_start = cds_end+1
	            if exon_start>= exon_end:
	                continue                
		    if chrom in UTR5_minus:
		        UTR5_minus[chrom].append((exon_start, exon_end))
		    else:
		        UTR5_minus[chrom] = [(exon_start, exon_end)]

        
    print "refGene parsed"
    infile.close()
    tree_dict = {}

    key_list = []
    for key in all_cds_plus:
       if key not in key_list:
           key_list.append(key)
    for key in all_cds_minus:
       if key not in key_list:
           key_list.append(key)
        
    #Create the interval trees, keys are chromosomes
    for key in key_list:
        # initialise with empty trees in case utr keys are empty
        tree_dict[key] = {"CDS+":IntervalTree([Interval(-1, 0)]),"CDS-":IntervalTree([Interval(-1, 0)]),"UTR5+":IntervalTree([Interval(-1, 0)]),"UTR5-":IntervalTree([Interval(-1, 0)]),"UTR3+":IntervalTree([Interval(-1, 0)]),"UTR3-":IntervalTree([Interval(-1, 0)])}
    for key in all_cds_plus:        
        tree = IntervalTree.from_tuples(all_cds_plus[key])
        tree_dict[key]["CDS+"] = tree

    for key in all_cds_minus:
        tree = IntervalTree.from_tuples(all_cds_minus[key])
        tree_dict[key]["CDS-"] = tree
        
    for key in UTR5_plus:
        tree = IntervalTree.from_tuples(UTR5_plus[key])
        tree_dict[key]["UTR5+"] = tree
        
    for key in UTR5_minus:
        tree = IntervalTree.from_tuples(UTR5_minus[key])
        tree_dict[key]["UTR5-"] = tree

    for key in UTR3_plus:
        tree = IntervalTree.from_tuples(UTR3_plus[key])
        tree_dict[key]["UTR3+"] = tree
        
    for key in UTR3_minus:
        tree = IntervalTree.from_tuples(UTR3_minus[key])
        tree_dict[key]["UTR3-"] = tree

    # Shelf: Tree
    tree_shelve = shelve.open("refGene_dict.shelf")
    tree_shelve.update(tree_dict)
    tree_shelve.close()

    print "tree made"
    return
