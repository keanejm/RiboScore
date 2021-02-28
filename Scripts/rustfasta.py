from Bio.SeqIO import to_dict, parse

def parsefasta(mRNA_sequences):
    in_seq_handle = open(mRNA_sequences)
    seqDict = to_dict(parse(in_seq_handle, "fasta"))
    in_seq_handle.close()
    return seqDict
