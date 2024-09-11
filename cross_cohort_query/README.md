# Cross-cohort Query

This folder holds scripts to query probes across cohorts, to validate `REINDEER`'s ability in transferring knowledge from one data set to another.

The basic idea is to query a list of probes acquired from our publication in two independent studies:
- Single-base substitutions in 77 lung adenocarcinoma samples analysed in [Seo, et al.](https://doi.org/10.1101/gr.145144.112);
- Gene fusions in 148 acute myeloid leukaemia samples analysed in [XXX]().

For single-base substitution query (folder `query_Seo_et_al/`), we:
1. Extracted the shared probes between our list from Tab S1 and the Tab S3 of publication [Seo, et al.](https://doi.org/10.1101/gr.145144.112), using the script `extract_shared_probe.R`;
2. Generated a fasta file with the shared probe list, with `make_shared_fa.R`;
3. Queried the probes on [Transipedia.org](Transipedia.org), on the index titled `LUAD SEO (154 experiments)` and with `Raw counts` as counting method;
4. Compared the query results with the ground truth with a post-query filter of $min\_{hit} \geq 3$, as detailed in `compare_shared_query.R`.

For gene fusion query (folder `query_leucegene/`), we:
...

Compared the query results with the ground truth with a post-query filter of $min\_hit \geq 3$, as detailed in `compare_query.R`.
