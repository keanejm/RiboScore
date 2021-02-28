# RiboScore
RiboScore is a tool for quantitatively accessing the quality of ribosome profiling data. The software takes Ribo-seq data that is mapped to the genome and assigns scores to data features that reflect the quality of the dataset, including the proportion of reads mapping to coding regions compared to untranslated regions, the heterogeneity of read distribution along coding regions, triplet periodicity and sequencing bias. The software is written in python (version 2.7.6). 

## RiboScore Inputs
1. Refseq genome fasta file: The complete human genome file (hg38) can be downloaded [here](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/). A compressed version of the file is also available in the inputs/sequence folder. The complete fasta file can be merged using the following command:
```sh 
cat hg38.gz* | zcat > hg38.fa 
```
 2. Annotation file in RefGene format: Available in the Inputs/Annotations folder. Files in this format can be downloaded [here](https://genome.ucsc.edu/cgi-bin/hgTables) by selecting 'all fields from selected table' as the output format.
 3. Path to sorted Bam files and their indexes (must be mapped to the genome): As an example, a compressed version of the data from [Calviello et al. 2016](https://www.nature.com/articles/nmeth.3688), mapped to the hg38 genome and sorted using samtools is available in the Inputs/Bam_Files folder. The raw data from this study can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73136). The complete bam file can be merged using the following command:
```sh 
cat SRR2433794.bam.gz* | zcat > SRR2433794.bam 
```
4. Path to offset files produced by plastid package: The file names must match those of the bam files apart from the suffix, which should be "_p_offsets.txt"
5. Path to output directory
6. Assign name to output files

## RiboScore Execution
Run RiboScore using the following command:
```python
python Scripts/riboscore.py /Inputs/Sequence/hg38.fa ./Inputs/Annotations/refGene.txt ./Inputs/Bam_Files/ ./Inputs/Offsets/ ./Outputs/Calviello16/ Calviello16
```


