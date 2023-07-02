# Python BBDuk Wrapper for Filtering Paired Reads to Complete, 2-primer Amplicons

More documentation coming soon, but for now, here is an example invocation of the script:

```
ampl-BBDuk.py \
--primer_fasta primer_sequences.fasta \
--in_path sample_L001_R1_001.fastq.gz \
--in2_path sample_L001_R2_001.fastq.gz \
--outm sample_amplicons_r1.fastq.gz \
--outm2 sample_amplicons_r2.fastq.gz \
--ktrim t
```
