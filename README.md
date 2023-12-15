# publication-ccle
scripts and other files to be able to reproduce  the CCLE Transipedia publication

## Prerequisite



### The CCLE Reindeer index

The core of the study is the CCLE index built by Reindeer, comprising 1019 samples. It is to huge to be stored here (236 GB). The better way is to ask for us to transfer ours. But if you want build your own index, you have to 

- **download** all the 1019 fastq of the study  (11 TB), the [sra-CCLE_metadata.tsv](./sra-CCLE_metadata.tsv)  file contain the list and links to dowloads the selection of fastq files.
- **trim** the fastq with cutadatp, options: **<TODO\>**
- **build the index** (prefer --disk-mode) using trimmed fastq. 
- in the index directory, create a tsv file named **fos.txt**. It must contains the list of the samples in the first column and the normalized value of the sample (computed using the kmer number of the fastq files of the sample) (needed by rdeer-service)
- If you have used the  ``--disk-query`` Reindeer index option in the directory index, create an empty file named **disk-query** (needed by rdeer-service).

### softwares and python environment

The applications below are needed:

- kmerator [https://github.com/Transipedia/kmerator](https://github.com/Transipedia/kmerator)
- rdeer-service [https://github.com/Bio2M/rdeer-service](https://github.com/Bio2M/rdeer-service)
- snakemake [https://github.com/snakemake/snakemake](https://github.com/snakemake/snakemake)

In addition, some python scripts needs external libraries:

- PyYAML
- pyfaidx
- pandas
- plotly
- seaborn 

All these are listed in the ``requirements.txt`` file. 
So, you can install python softwares and libraries with pip:

```
pip install -r requirements.txt
```

### R environment

R script needs the library ``Biostring``, show the documentation at [https://bioconductor.org/packages/release/bioc/html/Biostrings.html](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

## Load the CCLE index

Thanks to rdeer-service, you can load the index in memory, and run multiple queries against the index very quicly, rdeer-service embed a special version of  Reindeer capable of running in socket mode.
 
 We recommand to take a look on the README.md of rdeer-service to be able to load the CCLE index: [https://github.com/Bio2M/rdeer-service](https://github.com/Bio2M/rdeer-service).
 

## Using snakemake

As the mutation and chimera pipelines are launched with snakemake, it's needed to know 
some option of it. 

### create the DAG

create the pipeline DAG to visualize the sequence of process steps

```
snakemake --rulegraph | dot -Tsvg > dag.svg
```

### try blank test

 Simulating the pipeline launch allows you to check which rules will be executed and how often (and somtimes raise syntax errors)
 
 ```
 snakemake -n 
 ```

### run the pipeline

Finally, the pipeline will be started by

```
snakemake -j 12
```

## Mutation part
 
The mutations part is based on the Snakemake pipeline manager. 

- ``Snakefile``: the pipeline definition file
- ``config.yaml``: pipeline launch parameters and list of genes studied
- ``bin/``: directory containing the scritps required for the pipeline
 
### the config.yaml file
 
 The file ``config.yml`` is designed to modify the pipeline parameters, and also includes the list of genes to be studied.
 
#### Parameters

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


 
 



