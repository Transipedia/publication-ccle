# Exploring a large cancer cell line RNA-sequencing dataset with k-mers
scripts and other files to be able to reproduce the publication

## 1. Index preparation


### The CCLE Reindeer index

The core of the study is the CCLE index built by Reindeer, comprising 1019 samples. It is too large to be stored here (236 GB). A better way is to ask us to transfer it to you. But if you want build your own CCLE index, here is the procedure.

_NB: this procedure is suited to Reindeer V 1.02 or earlier._

- **download** all 1019 fastq files (11 TB), the [sra-CCLE_metadata.tsv](./sra-CCLE_metadata.tsv) file contain the list and links to download the selection of fastq files.
- **trim** the fastq with cutadatp, options: **<TODO\>**
- **build the index** (prefer --disk-mode) using trimmed fastq. 
- in the index directory, create a tsv file named **fos.txt**. It must contain the list of the samples in the first column and the normalized value of the sample (computed using the kmer number of the fastq files of the sample) (needed by rdeer-service). Follow the link below to rdeer-service for more details.
- If you have used the  ``--disk-query`` Reindeer index option in the directory index, create an empty file named **disk-query** (needed by rdeer-service).

### software and python environment

Applications below are needed:

- kmerator [https://github.com/Transipedia/kmerator](https://github.com/Transipedia/kmerator)
- rdeer-service [https://github.com/Bio2M/rdeer-service](https://github.com/Bio2M/rdeer-service)
- snakemake [https://github.com/snakemake/snakemake](https://github.com/snakemake/snakemake)

In addition, certain python scripts need external libraries:

- PyYAML
- pyfaidx
- pandas
- plotly
- seaborn 

All these are listed in the ``requirements.txt`` file, so, you can install all python software and libraries with pip:

```
pip install -r requirements.txt
```

### R environment

The R script needs the library ``Biostring``, see documentation at [https://bioconductor.org/packages/release/bioc/html/Biostrings.html](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

## Load the CCLE index

Rdeer-service loads the index in memory and enables running multiple queries against the index in real time. Rdeer-service embeds a special version of Reindeer capable of running in socket mode.
 
We recommand taking a look at the README.md of rdeer-service to learn how to load the CCLE index: [https://github.com/Bio2M/rdeer-service](https://github.com/Bio2M/rdeer-service).
 

## Using snakemake

The mutation pipeline is launched with snakemake. Thus some basic Snakemake usage is needed.  

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

## 2. Quantification analysis

The quantification folder holds multiple scripts evaluating Reindeer's quantification of gene expression or transposable elements. Please refer to the README.md in `quantification/` folder for more detailed information.

## 3. Mutation analysis
 
Mutations analysis is provided as a Snakemake pipeline. 

See README for detail. 

## 4. Fusion transcripts / chimera analysis

The fusion/chimera analysis is computed with a bash script as follows. 

- ``depmap_fusion_to_bed.sh``: main script to generate chimera probes, specific 31-mers from the input CCLE fusion annotation file [https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv](https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv) and count matrix by fusion events in all the samples from a reindeer index 
- ``config.sh``: script to declare annotation files
	- ``genome``= reference genome in fasta format
	- ``exons``= exon annotations from gencode in bed format
	- ``prime3 ``= bedfile of the 3' end of exons generated from the exons bed file
	-  ``prime5 ``= bedfile of the 5' end of exons generated from the exons bed file
	-  ``server ``= put the name of the server where the index is loaded
	-  ``index ``= name of the reindeer index
	
Once variables are configured, use the ``depmap_fusion_to_bed.sh`` as follow :

``depmap_fusion_to_bed.sh -a annotfile -c cclefile [-t threshold]``

With :

- ``annotfile``= corresponding depmap / SRR cell line names. A file is included into the chimera directory (`PRJNA523380_CCLE_metadata_RNAseq-name-to-srr-and-depmap.csv`)
- ``cclefile``= depmap fusion annotation file [https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv](https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv)
- ``threshold``= integer to filter events from depmap having less than x reads supporting the junction (0 by default)


 
 



