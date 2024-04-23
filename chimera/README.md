## Fusion / chimera analysis


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


