# RiboScore
RiboScore is a tool for quantitatively accessing the quality of ribosome profiling data. The software takes Ribo-seq data that is mapped to the genome and assigns scores to data features that reflect the quality of the dataset, including the proportion of reads mapping to coding regions compared to untranslated regions, the heterogeneity of read distribution along coding regions, triplet periodicity and sequencing bias. The software is written in python (version 2.7.6). 

## Riboscore Inputs
1. Refseq genome fasta file: The complete human genome file (hg38) can be downloaded [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/). A compressed version of the file is also available in the inputs/sequence folder. The complete fasta file can be merged using the following command:
```sh 
cat hg38.gz* | zcat > hg38.fa 
```
 2. Annotation file in RefGene format: 

