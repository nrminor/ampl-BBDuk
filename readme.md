# BBDuk primer matching
- This repository is for those conducting sequencing using matched paired end read primers.  
  - It filters reads that do not contain matched primers on both of the paired reads
- This repository is a wrapper for the bbmap/bbtools bbduk to improvement of identifying valid reads based on the primers
  - reference: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
  - reference: BBMap – Bushnell B. – sourceforge.net/projects/bbmap/
- As ran per the documentation, the bbduk will do a good job of identifying (and optionally trimming) primers based on an inputed primer list
- However, for paired reads the existing workflow does not identify the reads based on: 
  1. primer being evident in both of the paired reads per the matching criteria
  2. The paired reads both having the correct matched-paired primer 
- This workflow loops through the primer pair(s) and runs bbduk/repair.sh in order to meet the above criteria.

## Purpose
- This mitigates keeping reads that have 1 or 2 matched primers, but primers come from different primer pair/groupss.
- Most often these reads are a PCR chimera artifact, and should be excluded
- Since many are PCR chimeras, they are repeated artifacts
  - These invalid reads can be as high as 20% in some systems.
  - This means calling artifacts as variants at the >10% allele fraction.
- Filtering by length does not always work (setting min and max lengths) 
  - amplicons can be different lengths 
  - difficult to set a per-amplicon filtering criteria
  - pcr chimera can have a similar length.

## Use Case
- Using a matched-pair primer system to conduct genomic sequencing.
- Using matched-pair reads
- With the proper settings this workflow can be used for a variety sequencing techniques, and is not bound to a single vendor.
- 
# Requirements
- Linux or MacOsx (Windows OS will need modifications to run)
- python 3.6-3.8 (later versions are untested)
- bbtools (bbduk, bbmerge, bbmap, and repair.sh) installed to an excecutable location
  - https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
  - 
## How to Run:
- Download the files from this repository
- bbduk and repair.sh settings are hardcoded in the matched_pairs_bbduk_wrapper.py. Please modify the script as needed.
- Shown below is example files downloaded to /Volumes/drive1/bbduk_primer_matching/example

```shell
# bbduk and reformat.sh executables in /bin/bbmap, 8 GB ram allocated
python3 /Volumes/drive1/bbduk_primer_matching/matched_pairs_bbduk_wrapper.py \
--primer_fasta /Volumes/drive1/bbduk_primer_matching/example/KX601166.2_primers.fasta \
--in_path /Volumes/drive1/bbduk_primer_matching/example/zikv_dak_R1.fastq.gz \
--in2_path /Volumes/drive1/bbduk_primer_matching/example/zikv_dak_R2.fastq.gz \
--outm /Volumes/drive1/bbduk_primer_matching/example/zikv_dak_R1_matched.fastq.gz \
--outm2 /Volumes/drive1/bbduk_primer_matching/example/zikv_dak_R2_matched.fastq.gz \
--ktrim f
# bbmerge executables in /bin/bbmap, 8 GB ram allocated
java -ea -Xmx8000m -Xms8000m -Djava.library.path=/bin/bbmap/jni/ \
-cp /bin/bbmap/current/ jgi.BBMerge \
in=/Volumes/drive1/bbduk_primer_matching/example/zikv_dak_R1_matched.fastq.gz \
in2=/Volumes/drive1/bbduk_primer_matching/example/zikv_dak_R2_matched.fastq.gz \
out=/Volumes/drive1/bbduk_primer_matching/example/zikv_dak_matched_merge.fastq.gz
# bbmap executables in /bin/bbmap, 8 GB ram allocated
java -ea -Xmx8000m -Xms8000m -cp /bin/bbmap/current/ align2.BBMap build=1 \
in=/Volumes/drive1/bbduk_primer_matching/example/zikv_dak_matched_merge.fastq.gz \
ref=/Volumes/drive1/bbduk_primer_matching/example/KX601166.2.fasta \
outm=/Volumes/drive1/bbduk_primer_matching/example/zikv_dak_matched_merge.bam \

```
## required python packages:
- biopython
- os
- argparse
- subprocess
- gzip
## Format of the the input primer fasta
- The primer pairs must be indicated with the following nomenclature in the headers of each read:
  - <primer_name>_<primer_group>_<primer_direction>_<optional_additional_primer>
  - OR
  - <primer_name>_<primer_group>_<primer_direction>
  - Underscores are not allowed in the primer_name, primer_group
  - the primer_direction must be 'LEFT' or 'RIGHT'

### Example of input primer fasta:
```text
>ZIKVDAK_40_RIGHT
GCGTCAATATGCTGTTTTGCGT
>ZIKVDAK_40_LEFT
TGGTGCGTAGGATCATAGGTGA
>ZIKVDAK_39_RIGHT
CTGACTAGCAGGCCTGACAACA
>ZIKVDAK_39_LEFT
GCTGTGCCAATTGACTGGGTAC
```
## Example of input primer fasta using an alternate primer:
  - Some systems have a possible mutation in one or more primer regions
    - The use of alternative primers is required to properly amplify the amplicon.
```text
>ZIKVDAK_40_RIGHT
GCGTCAATATGCTGTTTTGCGT
>ZIKVDAK_40_LEFT_1
TGGTGCGTAGGATCATAGGTGG
>ZIKVDAK_40_LEFT_2
TGGTGCGTAGGATCATAGGTGA
>ZIKVDAK_39_RIGHT
CTGACTAGCAGGCCTGACAACA
>ZIKVDAK_39_LEFT
GCTGTGCCAATTGACTGGGTAC
```

## Notes on Efficiency
- This workflow is a looping wrapper around BBduk and repair.sh.
- It is VERY inefficient, introducing the matched pair criteria directly into bbduk would greatly increase the speed substantially
  - The speed increase would be about (1.5 * primer_pair_count) x, therefore a 40 primer pair system can expect to in <2% of the processing time
- maybe a match_primers flag
# Example of matched and mismatched-primer paired reads 
  - using Zikv KX601166.2 strain laboratory sample data
    - reference: https://www.ncbi.nlm.nih.gov/nuccore/KX601166.2
      - Shabman,R., Puri,V., Dilley,K., Fedorova,N., Shrivastava,S.,
            Amedeo,P., Williams,M., Hu,L. and Suthar,M.S
      - Submitted (25-JUL-2016) J. Craig Venter Institute, 9704 Medical
            Center Drive, Rockville, MD 20850, USA
  - Sequenceing conducted on a Illumina MiSeq instrument
    - "TruSeq" setting with read length = 251
    - Input Fastq Files have been randomly downsampled by 95% prior to uploading, due to github file size limits.
  - Publication: 
    - Block, L.N., Schmidt, J.K., Keuler, N.S. et al. Zika virus impacts extracellular vesicle composition and cellular gene expression in macaque early gestation trophoblasts. Sci Rep 12, 7348 (2022). https://doi.org/10.1038/s41598-022-11275-9
## Matched-primers paired-read example
    
Primers used:
```text
>ZIKVDAK_3_RIGHT
CCGGGGCAATCAACAATATCATGA
>ZIKVDAK_3_LEFT
TACAGATCATGGACCTCGGGCA
```
Forward Matched-primers read:
```text
@M01472:427:000000000-JWD6B:1:1111:9803:23202 1:N:0:1
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTCTCACTCCACGAGGAAGCTGCAAACGCGATCGCAGACTTGGCTAGAGTCAAGAGAATACACA
+
DDDDDFFFFFFFGGGGGGGGGGGHHHFHHHHHGGGGGHHHHHHHHHHHHHHGHHHGGHHHHHHHHGEGGGGGGHGHHHGHHHHFHHGGGGGGHGHHHHHHHHHHHGGGGGGHGHGGGHHHGGGGGHHHHGGGGGHHHHHHHHHHHHHHGGGGHHHHHHHGGGGGGHGHHHHHHHHHGFGGGGGFGGGGGGGGGGGGGGGGGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

```
Reverse complement Matched-primers read:
```text
@M01472:427:000000000-JWD6B:1:1111:9803:23202 2:N:0:1
CCGGGGCAATCAACAATATCATGACCAAGTATATAACTTTTTGGCTCGTCGAGCTTCCCAGGAGCCAGGCTATAGCTACAGCCACTAGCGCAAACCCGGGGTTCCTGAATATCCAATTTTCAACCTTGATCAGGTGCTTTGTGTATTCTCTTGACTCTAGCCAAGTCTGCGATCGCGTTTGCAGCTTCCTCGTGGAGTGAGAAGGAAGCGTCACGGCTCTTCTGGATCGTCGTGCTTCACCTTTTTTATGA
+
BCBBBCCCCFFFGGGGGGGGGGHHHHHGGHHHHHGHHHHHHHGGHHHGFGGGGGGHHHHHHGHGHGHHHHHHGHFHHHHHGHFHHHHFHFEGFGHGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHGHHHHHFHHHHHHHHHHHHHHHHHHHHHHHEFHHGHGGHGCFGFGGCHHHHHHGGFGGGGG?FGGBFGGGGGGGGGGDDFFBDDFFFFFFFFEFBFFAFDAEFFFFFFFFFFFEFF0

```
Merged Matched-primers read (length = 390):
```
@M01472:427:000000000-JWD6B:1:1111:9803:23202 1:N:0:1
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTCTCACTCCACGAGGAAGCTGCAAACGCGATCGCAGACTTGGCTAGAGTCAAGAGAATACACAAAGCACCTGATCAAGGTTGAAAATTGGATATTCAGGAACCCCGGGTTTGCGCTAGTGGCTGTAGCTATAGCCTGGCTCCTGGGAAGCTCGACGAGCCAAAAAGTTATATACTTGGTCATGATATTGTTGATTGCCCCGG
+
DDDDDFFFFFFFGGGGGGGGGGGHHHFHHHHHGGGGGHHHHHHHHHHHHHHGHHHGGHHHHHHHHGEGGGGGGHGHHHGHHHHFHHGGGGGGHGHHHHHHHHHHHGGGGGGHGHGGGHHHGGGGGHHHHGGGGGHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHHHHGHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHGGGGGGGGHGFGEFHFHHHHFHGHHHHHFHGHHHHHHGHGHGHHHHHHGGGGGGFGHHHGGHHHHHHHGHHHHHGGHHHHHGGGGGGGGGGFFFCCCCBBBCB
```
## Mismatched primer paired-read example
- This has a primer from two different primer groups on the ends
- It is likely a pcr chimera.
Primers used:
```text
>ZIKVDAK_11_RIGHT
CCTTTCCCTTGACAGCTGTTCC
>ZIKVDAK_3_LEFT
TACAGATCATGGACCTCGGGCA
```
Forward Mismatched primer read:
```
@M01472:427:000000000-JWD6B:1:1102:12301:25336 1:N:0:1
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTGTGGAGGATCACGGGTTTGGGATCTTCCACACCAGTGTTTGGCTGAAGGTCAGAGAGGACTA
+
CCCCBFFFFFFFGGGGGGGGGGGHHHHHHHHHGGGGGHHHHHHHHHHHGHHHHHHGGHHHHHHHHGGGGGEGGHGHHHGHHHHHHHGGGGGGHGHHHGGGHHHHHGGGGGGHGHFGGHHHDFGGGHHHHGGGGGHHHHHHHHHHHHHHCDGFGHHFGHHGGGGGGHFCHHHHHHGHGGAEGGGGGGGGGGGGGGGGGGGGFFFF@AFFEDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
Reverse complement Mismatched primer read
```
@M01472:427:000000000-JWD6B:1:1102:12301:25336 2:N:0:1
CCTTTCCCTTGACAGCTGTTCCTATGACGGCTGGGTCACACTCTAGTGAGTAGTCCTCTCTGACCTTCAGCCAAACACTGGTGTGGAAGATCCCAAACCCGTGATCCTCCACAAGGAAGCGTCACGGCTCTTCTGGATCGTCGTGCTTCACCTTTTTTATGATGACAGGTTCCGTACACAACCCAAGTCGATGTCGTGTTGCACCAGCAATCGACGTCATCTGGCTCCACTCCCTCGTCTAGCATGGGGCC
+
ABBCCFFFFFFFGGFGGGGGGGGHHHGHGGGGGGGHHHHHHHHHGFHAFEGHHHGHHHHHHHEHHHHHHGHGHGHGHGFHHHHHGGHHGHHFHHHHHHGGGGFGHHHHHHHGHEHEHHHHGGGGGGGGGGHHHHHHGHHGHHGGGGGHGHHHHHEHHGHHHHHGHHHHGGHHHGHHGGHFFGGGGHGHHGGGHHGHGCEFEFF0CFGGGGGFEFBCDGGGGGFFFFFFFF/FBE/B>FFFEFFFFFFFFF-

```
Merged Mismatched primer read=301:
```text
@M01472:427:000000000-JWD6B:1:1102:12301:25336 1:N:0:1
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTGTGGAGGATCACGGGTTTGGGATCTTCCACACCAGTGTTTGGCTGAAGGTCAGAGAGGACTACTCACTAGAGTGTGACCCAGCCGTCATAGGAACAGCTGTCAAGGGAAAGG
+
CCCCBFFFFFFFGGGGGGGGGGGHHHHHHHHHGGGGGHHHHHHHHHHHGH<JJJJJJJJJJJJJJJJJJJHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJEFAHFGHHHHHHHHHGGGGGGGHGHHHGGGGGGGGFGGFFFFFFFCCBBA
```
### Evidence of PCR Chimera Artifact in mismatched read:
```text
Merged Read
2915 - 3035
# Forward
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTGTGGAGGATCACGGGTTTGGGATCTTCCACACCAGTGTTTGGCTGAAGGTCAGAGAGGACTA
# Reverse complement
                                                  GGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTGTGGAGGATCACGGGTTTGGGATCTTCCACACCAGTGTTTGGCTGAAGGTCAGAGAGGACTACTCACTAGAGTGTGACCCAGCCGTCATAGGAACAGCTGTCAAGGGAAAGG
# Merged read
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTTGTGGAGGATCACGGGTTTGGGATCTTCCACACCAGTGTTTGGCTGAAGGTCAGAGAGGACTACTCACTAGAGTGTGACCCAGCCGTCATAGGAACAGCTGTCAAGGGAAAGG                                                                                          
# KX601166.2 551-740
TACAGATCATGGACCTCGGGCACATGTGTGACGCCACCATGAGTTATGAGTGCCCCATGCTAGACGAGGGAGTGGAGCCAGATGACGTCGATTGCTGGTGCAACACGACATCGACTTGGGTTGTGTACGGAACCTGTCATCATAAAAAAGGTGAAGCACGACGATCCAGAAGAGCCGTGACGCTTCCTT
# KX601166.2 2915-3035
                                                                                                                                                                                     GCTTCCTTGTGGAGGATCACGGGTTTGGGATCTTCCACACCAGTGTTTGGCTGAAGGTCAGAGAGGACTACTCACTAGAGTGTGACCCAGCCGTCATAGGAACAGCTGTCAAGGGAAAGG

```
