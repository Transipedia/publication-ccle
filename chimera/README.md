## Fusion / chimera analysis

This bash script performs the fusion analysis from manuscript "Exploring a large cancer cell line RNA-sequencing dataset with k-mers".

_NB: Reindeer queries are performed using Reindeer service (https://github.com/Bio2M/rdeer-service) running Reindeer 1.02. This allows multiple queries running in real time on a memory resident index. To run this pipeline on a local Reindeer index, merge all queries into a single query file, and run reindeer query once, to avoid incurring multiple index loads).

Main Steps of the script

* You first need to download the mandatory files (``annotfile`` and ``cclefile``) given as parameters
* Generate a file with corresponding sample IDs (SRR) + fusion events information
* Generate bed files for the left and right parts of the chimera + generate the corresponding 51nt long probes
* Filter to exon borders events
* Launch kmerator to mask non-specific probes / retain all unmasked k-mers
* Apply low complexity filter
* Run Reindeer query over CCLE index with filtered k-mers
* Return Reindeer output if at least 3 kmers per probe have a positive count + compute mean count per probe

Files

- ``depmap_fusion_to_bed.sh``: main script to generate chimera probes, specific 31-mers from the input CCLE fusion annotation file [https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv](https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv) and count matrix by fusion events in all the samples from a reindeer index 
- ``config.sh``: script to declare annotation files
- ``genome``= reference genome in fasta format
- ``exons``= exon annotations from gencode in bed format
- ``prime3 ``= bedfile of the 3' end of exons generated from the exons bed file
-  ``prime5 ``= bedfile of the 5' end of exons generated from the exons bed file
-  ``server ``= name of the reindeer server (where the index is loaded)
-  ``index ``= name of the reindeer index directory
	
Once variables are configured, use the ``depmap_fusion_to_bed.sh`` as follow :

``depmap_fusion_to_bed.sh -a annotfile -c cclefile [-t threshold]``

With :

- ``annotfile``= corresponding depmap / SRR cell line names. A file is included into the chimera directory (`PRJNA523380_CCLE_metadata_RNAseq-name-to-srr-and-depmap.csv`)
- ``cclefile``= depmap fusion annotation file [https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv](https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q2&filename=CCLE_fusions.csv)
- ``threshold``= integer to filter events from depmap having less than x reads supporting the junction (0 by default)


