# Exploring a large cancer cell line RNA-sequencing dataset with k-mers
Scripts and data files for reproducing the publication results

## 1. Index preparation


### The CCLE Reindeer index

The core of the study is the CCLE index built by Reindeer, comprising 1019 samples. It is too large to be stored here (236 GB). A better way is to ask us to transfer it to you. But if you want build your own CCLE index, here is the procedure.

_NB: this procedure is suited to Reindeer V 1.02 or earlier._

- **download** all 1019 fastq files (11 TB), the [sra-CCLE_metadata.tsv](./sra-CCLE_metadata.tsv) file contain the list and links to download the selection of fastq files.
- **trim** the fastq with ``cutadatp -q 10,10 -m 31 -p out2.fastq in1.fastq in2.fastq``
- **build the index** (prefer --disk-mode) using trimmed fastq: ``reindeer --index  --disk-query -f fof.unitigs.fa``
- in the index directory, create a tsv (tab-separated) file named **fos.txt**. It must contain the list of the samples in the first column and the normalized value of the sample (computed using the kmer number of the fastq files of the sample) (needed by rdeer-service). Follow the link below to rdeer-service for more details.
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

The R script complexity.R needs the library ``Biostring``, see documentation at [https://bioconductor.org/packages/release/bioc/html/Biostrings.html](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

### Loading the CCLE index in memory

Rdeer-service loads the index in memory and enables running multiple queries against the index in real time. Rdeer-service embeds a special version of Reindeer capable of running in socket mode.
 
We recommand taking a look at the README.md of rdeer-service to learn how to load the CCLE index: [https://github.com/Bio2M/rdeer-service](https://github.com/Bio2M/rdeer-service).
 

## 2. Quantification analysis

The quantification folder holds scripts for evaluating Reindeer's quantification of gene expression. See `quantification/README` for more detail.

## 3. Mutation analysis
 
Mutations analysis is provided as a Snakemake pipeline. 

See mutations/README for detail. 

## 4. Fusion / chimera analysis

The fusion/chimera analysis is computed with a bash script. See chimera/README for detail. 

## 5. Transposable element analysis

See transposable_element/README for detail

 



