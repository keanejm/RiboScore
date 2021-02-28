from os import listdir
import shelve
import pysam


def offset_set(pathToStudy,pathToOffsets):
	# Initialize Dictionary
	offset_dict = {}
	# Ensure one end slash
	if not pathToStudy[-1] == "/":
		pathToStudy += "/"
	if not pathToOffsets[-1] == "/":
		pathToOffsets += "/"

	FileCount = 0

	bam_files = listdir(pathToStudy)
	off_files = listdir(pathToOffsets)

	for treatmentDirContent in bam_files:
	    if treatmentDirContent[-4:] == ".bam":
		FileCount += 1
                infile = pysam.Samfile(pathToStudy + treatmentDirContent, 'rb')
		total_reads = 0
		mapped_reads = 0
		lengthofreads = [25,26,27,28,29,30,31,32,33,34]
		lengthofreadscount = [0,0,0,0,0,0,0,0,0,0]

		#For every read in the sam/bam file
		for read in infile.fetch(until_eof=True):
		    total_reads += 1
		    if not read.is_unmapped:
			#read is mapped
			mapped_reads += 1
			readlen = read.query_alignment_length
		        if readlen in lengthofreads:
			    lengthofreadscount[lengthofreads.index(readlen)] += 1
			else:
			    pass
		    else:
			pass

		max_count = max(lengthofreadscount)	  
		optimum_readlength = lengthofreads[lengthofreadscount.index(max_count)]
		optimum_range = str(optimum_readlength) + ":" + str(optimum_readlength)
			  
					  
		for offsetsDirContent in off_files:
	   	    if treatmentDirContent[:-4] == offsetsDirContent[:-14] and offsetsDirContent[-3:] == "txt":
			    offset_file = open(pathToOffsets + offsetsDirContent,"r")
			    readlength_offsets = []
			    for line in offset_file:
				if line[0] != "#" and line[0] != "l":
				    splitline = line.split("\t")
				    line_offset = int(splitline[1]) + 3 
				    readlength_offsets.append(line_offset)
				    if splitline[0] == str(optimum_readlength):
					optimum_offset = int(splitline[1]) + 3
			    readlength_offsets.append(optimum_range)
			    readlength_offsets.append(optimum_offset)
			    offset_file.close()
			    break
			
		if treatmentDirContent[:-4] not in offset_dict:
			offset_dict[treatmentDirContent[:-4]] = readlength_offsets
				        
		print "FileCount: {}".format(FileCount)
			
			
	# Shelf: Study offset dictionary
	offset_shelve = shelve.open("offset.shelf")
	offset_shelve.update(offset_dict)
	offset_shelve.close()
	return


















