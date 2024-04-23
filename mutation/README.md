## This Snakemake pipeline performs mutation analysis in manuscript "Exploring a large cancer cell line RNA-sequencing dataset with k-mers".
# Main Steps:
- Download Depmap files (rule: dnld_depmap)
- Extract useful information and links to SSR IDs (rule: depmap_csv_to_tsv)
- Split results as one file per gene (rule: vcf_by_gènes)
- Generate 61nt probes (rule: vcf2seq_61)
- Launch kmerator for masking probes and retain all unmasked k-mers (rule: kmerator)
- Apply low complexity filter (rule: complexity)
- Run Reindeer query over CCLE index using filtered k-mers (rule: kmerator)
- Return Reindeer output if at least 3-kmers are matched, compute mean count per probe (rule: merge_kmer)
- Discard results with mean counts <=5 (as in DepMap)
- Compute VAFs, create result table and graphics. 
