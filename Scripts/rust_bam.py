from os.path import exists
from os import mkdir
from pysam import Samfile
from sys import exit
from math import log  
import time

class Gene:
    def __init__(self, name, seq_chrom):
        self.name = name
        self.chrom_seq = seq_chrom

    def chrom(self):
        return self.name[2]

    def strand( self ):
        if self.name[3] in ["-","+"]:
            return self.name[3]
        else :
            exit("Error: strand name neither '-' nor '+'")

    def coordinates(self) :
        exon_number = int(self.name[8])
        exon_starts = self.name[9].split(",")
        exon_ends = self.name[10].split(",")
        list_coorr = []
        for number in range(exon_number):
            exon_start = int(exon_starts[number])
            exon_end = int(exon_ends[number])
            for number in range(exon_start, exon_end):
                list_coorr.append(number)

        list_coorr.sort()
        return set(list_coorr)

    def coordinates_50utr(self) :
        exon_number = int(self.name[8])
        exon_starts = self.name[9].split(",")
        exon_ends = self.name[10].split(",")
        list_coorr = []
        for number in range(exon_number):
            exon_start = int(exon_starts[number])
            exon_end = int(exon_ends[number])
            if number == 0 :
                exon_start = int(exon_starts[number])-50

            if number == exon_number-1:
                exon_end = int(exon_ends[number])+50

            for number in range(exon_start, exon_end):
                list_coorr.append(number)

        list_coorr.sort()
        if self.name[3] == "-":
            list_coorr.reverse()

        return list_coorr

    def coordinates_coding_50utr(self) :
        cds_start50,cds_end50 = int(self.name[6]),int(self.name[7])
        if cds_start50 == cds_end50: return []
        exon_number = int(self.name[8])
        exon_starts = self.name[9].split(",")
        exon_ends = self.name[10].split(",")
        list_coorr = []
        for number in range(exon_number):
            exon_start = int(exon_starts[number])
            exon_end = int(exon_ends[number])
            if exon_start <= cds_start50 :
                exon_start = cds_start50-50

            if exon_start> exon_end :
                continue

            if exon_end >= cds_end50:
                for number in range(exon_start, cds_end50+50):
                    list_coorr.append(number)

                break    

            for number in range(exon_start, exon_end):
                list_coorr.append(number)

        list_coorr.sort()
        if self.name[3] == "-":
            list_coorr.reverse()

        return list_coorr

    def sequence_cds(self) :
        cds_start50,cds_end50 = int(self.name[6]),int(self.name[7])
        exon_number = int(self.name[8])
        exon_starts = self.name[9].split(",")
        exon_ends = self.name[10].split(",")
        list_coorr = []
        seq = ""
        for number in range(exon_number):
            exon_start = int(exon_starts[number])
            exon_end = int(exon_ends[number])
            if exon_start < cds_start50 :
                exon_start = cds_start50

            if exon_end > cds_end50:
               exon_end = cds_end50

            if exon_start> exon_end :
                continue

            seq += self.chrom_seq[self.name[2]][exon_start:exon_end]

        if self.name[3] == "+": 
            return seq

        else : 
            return seq.reverse_complement()


    def sequence(self) :
        exon_number = int(self.name[8])
        exon_starts = self.name[9].split(",")
        exon_ends = self.name[10].split(",")
        list_coorr = []
        seq = ""
        for number in range(exon_number):
            exon_start = int(exon_starts[number])
            exon_end = int(exon_ends[number])
            seq += self.chrom_seq[self.name[2]][exon_start:exon_end]
        if self.name[3] == "+":
            return seq

        else :
            return seq.reverse_complement()

    def tx_start_end(self) :
        return int(self.name[4]),int(self.name[5])

    def cds_start_end(self) :
        return int(self.name[6]),int(self.name[7])

    def alternative_name(self):
        return self.name[11]

def stop_err(msg):
    exit("%s\n" % msg)                  # Prints msg to stderr and exits with code 1

def mean_value(input_list):
    return(sum(map(float, input_list))/len(input_list))  

universal_code = {
"Ala" : ["GCT", "GCC", "GCG", "GCA"],
"Gly" : ["GGT", "GGC", "GGG", "GGA"],
"Pro" : ["CCT", "CCC", "CCG", "CCA"],
"Thr" : ["ACT", "ACC", "ACG", "ACA"],
"Val" : ["GTT", "GTC", "GTG", "GTA"],
"Ser" : ["TCT", "TCC", "TCG", "TCA", "AGT", "AGC"],
"Arg" : ["CGT", "CGC", "CGG", "CGA", "AGG", "AGA"],
"Leu" : ["CTT", "CTC", "CTG", "CTA", "TTG", "TTA"],
"Phe" : ["TTT", "TTC"],
"Asn" : ["AAT","AAC"],
"Lys" : ["AAG", "AAA"],
"Asp" : ["GAT", "GAC"],
"Glu" : ["GAG", "GAA"],
"His" : ["CAT","CAC"],
"Gln" : ["CAG", "CAA"],
"Ile" : ["ATT", "ATC", "ATA"],
"Met" : ["ATG"],
"Tyr" : ["TAT","TAC"],
"Cys" : ["TGT", "TGC"],
"Trp" : ["TGG"],
"Stop" : ["TGA","TAG","TAA"]}

def writerustfiles(seqDict, refGene, sortedBam, offset, readlen_range, outDir):     
    chromosome_list = seqDict.keys()

    readlen_rangesplit = readlen_range.split(":")
    if len(readlen_rangesplit) == 1: 
        accepted_read_lengths = [int(readlen_rangesplit[0])]
        length_values ="%s"%int(readlen_rangesplit[0])

    elif len(readlen_rangesplit) == 2 :
        accepted_read_lengths = [readlen for readlen in range(int(readlen_rangesplit[0]),int(readlen_rangesplit[1])+1)]
        length_values ="%s_%s"%(int(readlen_rangesplit[0]),int(readlen_rangesplit[1]))

    else : 
        stop_err("Lengths of footprints parameter not in correct format, it should be either colon seperated with the second value greater or equal to the first, (28:32) or a single interger (31)")
    if len(accepted_read_lengths) == 0  :
        stop_err( "Lengths of footprints parameter not in correct format, it should be either colon seperated with the second value greater or equal to the first, (28:32) or a single interger (31)")

    nts = ["A","G","C","T"]
    alignments_A1 = Samfile(sortedBam, 'rb')
    codon_enrichment_dict = {}
    codon_enrichment_expected_dict= {}
    for nt in nts :
        for nt2 in nts :
            for nt3 in nts :
                codon = "%s%s%s"%(nt, nt2, nt3)
                codon_enrichment_dict[codon] = {}
                codon_enrichment_expected_dict[codon] = []
                for number in range(0,60,1) :
                    codon_enrichment_dict[codon][number] = [0.0,0.0]

    infileopen = open(refGene)
    infileopen.seek(0)
    infileopen.readline()

    subcodonp1 = {}
    for n in range(20,51) :
      subcodonp1[n] = {0:0,1:0,2:0}

    longest_transcript = {}
    used_transcripts = []
    skip_list = []
    for line in infileopen:
        linesplit = line[:-1].split("\t")
        gene_name = linesplit[1]
        gene_name2 = linesplit[12]
        G = Gene(linesplit, seqDict)
        cdsstart, cdsend =  G.cds_start_end()
        coding_codordinates2 = G.coordinates_coding_50utr()
        
        if gene_name in used_transcripts: 
            skip_list.append(gene_name)
            continue
        used_transcripts.append(gene_name)
        
        if G.chrom() not in chromosome_list : 
            continue
            #skip_list.append(gene_name)
        if len(coding_codordinates2)  != len(set(coding_codordinates2)):
            continue
            #skip_list.append(gene_name)
        #if len(G.chrom()) > 6 : 
            #continue
            

        coding_codordinates2set = set(coding_codordinates2)
        if len(coding_codordinates2set) < 330 or len(coding_codordinates2set)%3 != 1 or cdsstart < 50 :
            continue                                # takes only those whose ORF is divisible by three

        longest_transcript.setdefault(gene_name2,[]).append((len(coding_codordinates2set),gene_name))

    list_longest = []
    for key, value in longest_transcript.items() :
        value.sort()
        list_longest.append(value[-1][1])

    list_longest= set(list_longest)     

    infileopen.seek(0)
    infileopen.readline()
    list_chrom = []
    for line in infileopen:
        linesplit = line[:-1].split("\t")
        gene_name = linesplit[1]
        if gene_name not in list_longest:
            continue

        if gene_name in skip_list :
            continue

        G = Gene(linesplit, seqDict)
        #if len(G.chrom()) > 6 :
        #print G.chrom()
            #continue

        cdsstart, cdsend = G.cds_start_end()
        sequence = str(G.sequence_cds().seq)
        list_chrom.append(G.chrom())

        coding_codordinates2 = G.coordinates_coding_50utr()
        coding_codordinates2set = set(coding_codordinates2)

        selected_region = sequence[120:-60]

        all_reads = alignments_A1.fetch(G.chrom(), cdsstart-50,cdsend+50)
        profile_list = [0.0 for n in range(len(coding_codordinates2))]

        period_dict = {}
        for n in range(20,51) :
            period_dict[n] = profile_list[:]

        for read in all_reads :
            readlen = read.qlen
            if readlen < 20 :
                continue 
            if readlen > 50 :
                continue 
            if G.strand() == "+" :
                if read.is_reverse :
                    continue

                cpositions = read.positions
                cpositions.sort()
                                
                # Specific Inferred Offset for Metafootprint Profile
                shelf_readlengths = [25,26,27,28,29,30,31,32,33,34]
                if readlen in shelf_readlengths:
                    coordinate = cpositions[offset[shelf_readlengths.index(readlen)]]
                else:
                    coordinate = cpositions[offset[10]]
                                                
                if coordinate in coding_codordinates2set:
                    if readlen in accepted_read_lengths : 
                        profile_list[coding_codordinates2.index(coordinate)] += 1
                
                    period_dict[read.qlen][coding_codordinates2.index(coordinate)] += 1

            else:
                if not read.is_reverse :
                    continue
                cpositions = read.positions
                cpositions.sort()
                                
                # Specific Inferred Offset for Metafootprint Profile
                shelf_readlengths = [25,26,27,28,29,30,31,32,33,34]
                if readlen in shelf_readlengths:
                    coordinate = cpositions[-1-offset[shelf_readlengths.index(readlen)]]
                else:
                    coordinate = cpositions[-1-offset[10]]
                
                if coordinate in coding_codordinates2set:
                    if readlen in accepted_read_lengths : 
                        profile_list[coding_codordinates2.index(coordinate)] += 1
                        
                    period_dict[read.qlen][coding_codordinates2.index(coordinate)] += 1

        for n in range(20,51) :
            cds1 = period_dict[n][50:-50]
            for subcodon, value in enumerate(cds1 ):
                subcodonp1[n][(subcodon) % 3 ] += value

        profile_list = profile_list[170:-110]

        average_gene_density = float(sum(profile_list))/len(profile_list)
        if average_gene_density != 0 :
            num_codon = len([1 for number88 in range(0,len(profile_list),3) if ((profile_list[number88]+profile_list[number88+1]+profile_list[number88+2])/3)>average_gene_density])
            expected_codon_density = float(num_codon)/(len(profile_list)/3) 

            codon_sequence = sequence
            codon_start = 0
            for sliding_w_n in range(0,len( selected_region),3) :
                codon_window = str(codon_sequence[codon_start:codon_start+180])
                if len(set(codon_window) - set(["A","T","G","C"])) != 0: 
                    codon_start     += 3
                    continue

                if (profile_list[sliding_w_n]+profile_list[sliding_w_n+1]+profile_list[sliding_w_n+2])/3 > average_gene_density:
                    for number in range(0,60) :
                        codon = codon_window[number*3:(number+1)*3]
                        codon_enrichment_dict[codon][number][0] += 1
                        codon_enrichment_dict[codon][number][1] += 1

                else :
                    for number in range(0,60) :
                        codon = codon_window[number*3:(number+1)*3]
                        codon_enrichment_dict[codon][number][0] += 1

                codon = codon_window[40*3:(40+1)*3]
                codon_enrichment_expected_dict[codon].append(expected_codon_density)
                codon_start += 3

    # Temporary fix for division by zero in mean_value(codon_enrichment_expected_dict[codon])
    #if len(codon_enrichment_expected_dict[codon]) == 0:
        #exit("Error for %s: Length of codon_enrichment_expected_dict[codon] is zero" % (sortedBam))

    if not exists(outDir):
        mkdir(outDir)

    alignment_filename = "%s_%s_%s" % (sortedBam.split("/")[-3],sortedBam.split("/")[-2],sortedBam.split("/")[-1])

    periodicity_file = open("%s/periodicity_%s"%(outDir, alignment_filename),"w")
    for n in range(20,51):
      periodicity_file.write("%s\t%s\t%s\t%s\n"%(n, subcodonp1[n][0],subcodonp1[n][1],subcodonp1[n][2]))

    list_chrom = list(set(list_chrom))
    read_count = 0
    for chrom in list_chrom:
        all_reads = alignments_A1.fetch(chrom)
        for read in all_reads: read_count += 1

    periodicity_file.write("mapped_reads\t%s\n"%(read_count))
    periodicity_file.close()

    outfile = open("%s/RUST_codon_file_%s_%s_%s" % (outDir, alignment_filename, offset, length_values), "w")
    outfile.write("codon, expected value")
    for number106 in range(-40, 20): 
        outfile.write(", %s" % number106)

    outfile.write("\n")

    list_codons = []
    codons = codon_enrichment_dict.keys()
    codons.sort()
    rust_expected = []
    rust_observed_metafootprint = []
    for codon in codons:
        if codon in list_codons:
            continue

        if codon in ["TAA","TGA","TAG"]:
            continue

        list_codons.append(codon)
        outfile.write("%s" % codon)
        if codon_enrichment_expected_dict[codon] != []:
            outfile.write( ", %s" % mean_value(codon_enrichment_expected_dict[codon]))
            rust_expected.append(mean_value(codon_enrichment_expected_dict[codon]))
        else:
            outfile.write( ", ND")
            rust_expected.append("ND")

        list_data = []
        for number in range(0, 60) :
            if codon_enrichment_dict[codon][number][0] != 0:
                outfile.write(", %s"%(codon_enrichment_dict[codon][number][1] / codon_enrichment_dict[codon][number][0]))
                list_data.append(codon_enrichment_dict[codon][number][1] / codon_enrichment_dict[codon][number][0])

            else :
                outfile.write(", 0")
                list_data.append(0)

        outfile.write("\n")
        rust_observed_metafootprint.append(list_data)

    shannon_values = []
    if "ND" in rust_expected:
        shannon_values = ["ND" for loc_i in range(60)]
    else:
        rust_expected_sum = sum(rust_expected)
        q_values = [n/rust_expected_sum for n in rust_expected]

        for loc_i in range(60) :
            rust_observed = [n[loc_i] for n in rust_observed_metafootprint ]
            rust_observed_sum = sum(rust_observed)
            rust_observed_min = min(rust_observed)
            if rust_observed_min == 0:
                shannon_values.append("ND")

            else:
                p_values = [n/rust_observed_sum for n in rust_observed]
                shannon = []
                list_normalised = []
                for p_value,q_value in zip(p_values, q_values):
                    shannon.append(abs(p_value * log((p_value/q_value), 2)))
                    list_normalised.append(p_value / q_value)

                shannon_values.append(sum(shannon))

    outfile.write("\nKullback Leibler divergence,")       
    for value in shannon_values:
        outfile.write(", %s"%value)

    outfile.close()

    return 0

