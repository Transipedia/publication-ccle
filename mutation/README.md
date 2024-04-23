## Mutation analysis
This Snakemake pipeline performs mutation analysis in manuscript "Exploring a large cancer cell line RNA-sequencing dataset with k-mers".

_NB: Reindeer queries are performed using Reindeer service (https://github.com/Bio2M/rdeer-service) running Reindeer 1.02. 
This allows multiple queries running in real time on a memory resident index. To run this pipeline on a local Reindeer index, merge all queries into a single query file, and run reindeer query once, to avoid incurring multiple index loads).  

### Main Steps:
- Download Depmap files (rule: dnld_depmap)
- Extract useful information and links to SSR IDs (rule: depmap_csv_to_tsv)
- Split results as one file per gene (rule: vcf_by_gènes)
- Generate 61nt probes (rule: vcf2seq_61)
- Launch kmerator for masking probes and retain all unmasked k-mers (rule: kmerator)
- Apply low complexity filter (rule: complexity)
- Run Reindeer query over CCLE index using filtered k-mers (rule: reindeer)
- Return Reindeer output if at least 3-kmers are matched, compute mean count per probe (rule: merge_kmer)
- Discard results with mean counts <=5 (as in DepMap)
- Compute VAFs, create result table and graphics. 

### Files

- ``Snakefile``: the pipeline definition file
- ``config.yaml``: pipeline launch parameters and list of genes studied
- ``bin/``: directory containing the scritps required for the pipeline
- ``config.yml``: file used to set the pipeline parameters and list of genes to be studied.
 
### Config.yml Parameters

```
### General parameters
bin_dir:     # location of scripts
base_dir : # result base directory
thread :    # the multithreaded commands will use this parameters
output_dir : # relative to base_dir, change it for multiple analysis
version :  # some output file use the version in their name

### rules parameters

depmap:
  filter_af  :  # to set a python test in the column 30 (named RNASeq_AC)
  isCOSMIChotspot : #  if True, keep only probes found in the  'isCOSMIChotspot' column

complexity:
  enabled:   # if True, complexity rule is applied

filter:
  enabled:	# if True, filter rule is applied
  recur-max:  # maximum recurrency to considere probes as positive.
  abund-min: # minimum count of probe to considere it as positive

merge_kmer:
  min_pos:  # Minimum number of positive counts when merging kmer as probe

kmerator:
  release : # transcriptome release used for kmerator app
  jellyfish:  # jellyfish genome index

genome:     # genome fasta file

annot:
  version:  # genome annotation version (ex: GRCh37)
  release:  # genome annotation release (ex: 42)

reindeer:
  index:  # Reindeer index directory name
  server: # Reindeer server name
  options: # additionnal options

genes:		# list of interrest genes
  - APC
  - ARID1A
  - ATM
  - ATRX
```
### create the DAG

Create the pipeline's DAG to visualize process steps

```
snakemake --rulegraph | dot -Tsvg > dag.svg
```

### try a blank test

Simulating the pipeline launch allows you to check which rules will be executed and how often (and sometimes raise syntax errors)
 
 ```
 snakemake -n 
 ```

### run the pipeline

Finally, the pipeline will be started by

```
snakemake -j 12
```
