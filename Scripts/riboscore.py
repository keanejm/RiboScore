import shelve
from sys import argv, exit
from os import mkdir, listdir, remove
from os.path import exists
from shutil import rmtree
from refGene_shelf import refGeneLocations
from offset_shelf import offset_set
from parallel import rustforker, gbdforker
from rustfasta import parsefasta
from riboscore_plot import RiboScore_PDF, RiboScore_Logo
import pandas as pd
from math import sqrt

def wpearson(vec_1, vec_2, weights, r = 2):
	
	list_length = len(vec_1)
        
	try:
		weights = list(map(float,weights))
	except:
		print('Invalid weights.')
		sys.exit(1)
	
	try:
		vec_1 = list(map(float,vec_1))
		vec_2 = list(map(float,vec_2))
		if len(vec_2) != list_length or len(weights) != list_length:
		    print('Vector/Weight sizes not equal.')
		    sys.exit(1)
	except:
		print('Invalid vectors.')
		sys.exit(1)
		        
	w_sum = sum(weights)
	vec1_sum = 0.0
	vec2_sum = 0.0
	for x in range(len(vec_1)):
		vec1_sum += (weights[x] * vec_1[x])
		vec2_sum += (weights[x] * vec_2[x])
	vec1_avg = (vec1_sum / w_sum)
	vec2_avg = (vec2_sum / w_sum)

	sum_top = 0.0
	sum_bottom1 = 0.0
	sum_bottom2 = 0.0

	for x in range(len(vec_1)):
		dif_1 = (vec_1[x] - vec1_avg)
		dif_2 = (vec_2[x] - vec2_avg)
		sum_top += (weights[x] * dif_1 * dif_2)
		sum_bottom1 += (dif_1**2)*(weights[x])
		sum_bottom2 += (dif_2**2)*(weights[x])

	cor = sum_top / (sqrt(sum_bottom1 * sum_bottom2))

	return round(cor,r)


def runrust(pathToFastaFile, pathToRefgene, pathToStudy, pathToOffsets, pathToPDF, StudyName):
    # Check if paths are valid
    if not exists(pathToFastaFile):
        exit("Error: %s does not exist" % (pathToFastaFile))
    if not exists(pathToRefgene):
        exit("Error: %s does not exist" % (pathToRefgene))
    if not exists(pathToStudy):
        exit("Error: %s does not exist" % (pathToStudy))
    if not exists(pathToPDF):
        mkdir(pathToPDF)

    # Ensure one end slash
    if not pathToStudy[-1] == "/":
        pathToStudy += "/"

    # Make temporary directory to store intermediate files
    outputDir = "temporary_rust_folder_" + pathToStudy.split("/")[-2]
    if exists(outputDir):
        exit("Error: %s already exists - please rename or remove" % (outputDir))
    mkdir(outputDir)

    # Parse the fasta/annotation/genome file
    seqDict = parsefasta(pathToFastaFile)
                                    
    # Shelf Offsets
    offset_set(pathToStudy,pathToOffsets)

    # Open Study Offsets Dictionary
    study_offset_shelve = shelve.open("offset.shelf")
    offsets = dict(study_offset_shelve)
    study_offset_shelve.close()
    
    # Collate utr and cds locations to dictionary
    refGeneLocations(pathToRefgene)
    
    # Open tree_dict from shelf
    tree_shelve = shelve.open("refGene_dict.shelf")
    tree_dict = dict(tree_shelve)
    tree_shelve.close()
    
    # Process GBD/Gini  in parallel
    gbdforker(pathToStudy, offsets, tree_dict, 20)
        
    # Run rust in parallel
    rustforker(seqDict, pathToRefgene, pathToStudy, offsets, "25:34", outputDir, 20)
        
    ### Calculate Scores
    feature_file = open("rust_tp.csv", "w")
    temp_dir = listdir(outputDir)
    for aFile in temp_dir:
	if "RUST_codon_file_" in aFile:
	    filenamelist1 = aFile.split("_")
	    file_id = filenamelist1[6]
	    if file_id.endswith(".bam"):
		file_id = file_id[:-4]
	    filename1 = " ".join(filenamelist1[3:-4])
	    filename1 = filename1[:-4]
	    for bFile in temp_dir:
	        if "periodicity_" in bFile:
		    filenamelist2 = bFile.split("_")
		    filename2 = " ".join(filenamelist2[1:-1])
		    filename2 = filename2[:-4]
	            if filename1 == filename2:
		        codon_file = open(outputDir + "/"+ aFile,"r")
			for line in codon_file:
			    linesplit = line.split(",")
			    if linesplit[0] == "Kullback Leibler divergence":
				linesplit[1] = file_id
				kld_values = map(float, linesplit[2:62])
			# Calculate RUST Score
			vec_peak = 0.5
			weight_peak = 4
			ideal_vec = [0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
			weight_vec = [1,1,1,1,1,1,1,1,1,1,1] 
			test_vec = [float(i) for i in linesplit[37:48]]
			Asite_peak = [float(linesplit[41]),float(linesplit[42]),float(linesplit[43])]
			peak_5end = [float(linesplit[37]),float(linesplit[38]),float(linesplit[39])]
			peak_3end = [float(linesplit[47]),float(linesplit[46]),float(linesplit[45])]
			peak_ends = peak_5end + peak_3end
			if Asite_peak[0] == max(Asite_peak) and Asite_peak[1] == max(Asite_peak):
			    ideal_vec[Asite_peak.index(max(Asite_peak))+5] = vec_peak
			    weight_vec[Asite_peak.index(max(Asite_peak))+5] = weight_peak
			elif Asite_peak[0] == max(Asite_peak) and Asite_peak[2] == max(Asite_peak):
			    ideal_vec[Asite_peak.index(max(Asite_peak))+5] = vec_peak
			    weight_vec[Asite_peak.index(max(Asite_peak))+5] = weight_peak
			else:
			    ideal_vec[Asite_peak.index(max(Asite_peak))+4] = vec_peak
			    weight_vec[Asite_peak.index(max(Asite_peak))+4] = weight_peak
			weight_vec[peak_5end.index(max(peak_5end))] = weight_peak
			if peak_3end[0] == max(peak_3end):
			    weight_vec[peak_3end.index(max(peak_3end))+10] = weight_peak
			elif peak_3end[1] == max(peak_3end):
			    weight_vec[peak_3end.index(max(peak_3end))+8] = weight_peak
			else:
			    weight_vec[peak_3end.index(max(peak_3end))+6] = weight_peak
			WP_coef_4 = wpearson(ideal_vec,test_vec,weight_vec)
			if WP_coef_4 < -0.50:
			    WP_coef_4 = -0.50
			rust_sc = round((WP_coef_4 - -0.5) / 1.5, 2)
			linesplit.append(" {}".format(rust_sc))
			newline1 = ",".join(linesplit[1:])
			codon_file.close()
	                
                        # Calculate Periodicity Score
			triplet_file = open(outputDir + "/"+ bFile,"r")
			list_0, list_1, list_2 = [], [], []
			for line in triplet_file:
			    linesplit = line.split("\t")
			    if linesplit[0] != "mapped_reads" and int(linesplit[0]) in range(25,35):
				list_0.append(float(linesplit[1]))
				list_1.append(float(linesplit[2]))
				list_2.append(float(linesplit[3]))
			pos_0, pos_1, pos_2 = sum(list_0), sum(list_1), sum(list_2)
			tri_total = sum([pos_0,pos_1,pos_2])
			prop_pos_0, prop_pos_1, prop_pos_2 = pos_0 / tri_total, pos_1 / tri_total, pos_2 / tri_total
			prop_list = [prop_pos_0,prop_pos_1,prop_pos_2]
			max_prop = max(prop_list)
			min_prop = min(prop_list)
			for i in prop_list:
			    if i != max_prop and i != min_prop:
				mid_prop = i
			gini_purity = (max_prop**2 + mid_prop**2 + min_prop**2)
			OldRange = (0.667)
			NewRange = (1.000)  
			gini_purity_std = (((gini_purity - 0.333) * NewRange) / OldRange) 
			triplet_score = ((max_prop - min_prop) + (mid_prop - min_prop) + (min_prop - min_prop))
			tp_sc = round(sum([gini_purity_std,triplet_score])/2,2)
			parameterlist = [str(max_prop),str(mid_prop),str(min_prop),str(tp_sc)]
			newline2 = ",".join(parameterlist)
			
			triplet_file.close() 
		    
			feature_file.write(newline1 + "," + newline2 + "\n")
    feature_file.close()
			###
    
    # Set up dataframe
    df1 = pd.read_csv("rust_tp.csv",header=None)
    df1.columns = ["C{}".format(x) for x in range(66)]
    df1.set_index('C0',inplace=True)
    df2 = pd.read_csv("gbd_gini.csv",header=None)
    df2.columns = ["C0","C66","C67","C68","C69","C70"]
    df2.set_index('C0',inplace=True)
    df = df1.join(df2)

    # Open dictionary gini from shelf
    gini_shelve = shelve.open("gini.shelf")
    gini_dict = dict(gini_shelve)
    gini_shelve.close()
    
    #Plot Functions
    RiboScore_PDF(df,gini_dict,pathToPDF,StudyName)
    RiboScore_Logo(df,pathToPDF)
    
    # Remove temporary folder
    rmtree(outputDir)
    remove("rust_tp.csv")
    remove("gbd_gini.csv")
    remove("gini.shelf")
    remove("refGene_dict.shelf")
    remove("offset.shelf")
    remove("RiboScore.csv")

def main():
    if len(argv) != 7: 
        exit("Usage: 6 arguments required\n1: Path to fasta file of genome\n2: Path to annotation file of coding regions in refGene format\n3: Path to folder containing folders containing bed and sorted bam files\n4: Path to Output Folder\n5: Name of Study")

    # Get paths
    pathToFastaFile = argv[1]
    pathToRefgene = argv[2]
    pathToStudy = argv[3]
    pathToOffsets = argv[4]
    pathToPDF = argv[5]
    StudyName = argv[6]

    runrust(pathToFastaFile, pathToRefgene, pathToStudy, pathToOffsets,pathToPDF,StudyName)

if __name__ == "__main__":
    main()
