# Cross-cohort Query

This folder holds scripts to query probes across cohorts, to validate `REINDEER`'s ability in transferring knowledge from one data set to another.

The basic idea is to query a list of probes acquired from our publication in two independent studies:
- Single-base substitutions in 77 lung adenocarcinoma samples analysed in [Seo, et al.](https://doi.org/10.1101/gr.145144.112);
- Gene fusions in 148 acute myeloid leukaemia samples analysed in [MacRae, et al.](https://doi.org/10.1371/journal.pone.0072884), [Pabst, et al.](https://doi.org/10.1182/blood-2015-11-683649) and [Lavallée, et al.](https://doi.org/10.1182/blood-2016-03-703868).

For single-base substitution query (folder `query_Seo_et_al/`), we:
1. Extracted the shared probes between our list from Tab S1 and the Tab S3 of publication [Seo, et al.](https://doi.org/10.1101/gr.145144.112), using the script `extract_shared_probe.R`, and generated a fasta file with the shared probe list, with `make_shared_fa.R`;
2. Queried the probes on [Transipedia.org](Transipedia.org), on the index titled `LUAD SEO (154 experiments)` and with `Raw counts` as counting method;
3. Compared the query results with the ground truth with a post-query filter of min_hit $\geq 3$, as detailed in `compare_shared_query.R`.

For gene fusion query (folder `query_leucegene/`), we:
1. Extracted the shared probes between our list from Tab S6 and the publications [MacRae, et al.](https://doi.org/10.1371/journal.pone.0072884), [Pabst, et al.](https://doi.org/10.1182/blood-2015-11-683649) and [Lavallée, et al.](https://doi.org/10.1182/blood-2016-03-703868);
2. Queried the probes on [Transipedia.org](Transipedia.org), on the index titled `Leucegene with metadata (148 experiments)` and with `Raw counts` as counting method;
3. Compared the query results with the ground truth with a post-query filter of $min\_hit \geq 3$, as detailed in `compare_query.R`.
