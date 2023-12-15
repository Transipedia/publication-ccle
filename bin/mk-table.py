#!/usr/bin/env python3


"""
From the fasta file by gene from Mytailorisrich

1. `Total Kmer-Alt`: total count of different kmers in the fasta file by gene
        nombre de kmers de la sortie de mytailorisrich (nombre de ligne d'entête)
2. `non-overlap Kmer-Alt`: count of kmer-Alt in exons
        nombre de kmers après avoir retirer les kmers qui sont chevauchant exon/intron
3. `overlap Kmer-Alt`: count of kmer-Alt in exon/intron boundaries
        nombre de kmers chevauchant exon/intron
4. `SRR hit max`:
        nombre de SRR maximum représentant un seul kmer
5. `Total SRR`: number of different SRR in the fasta file by gene
        nombre de SRR impliqués par gène (SRR distincts des fichiers 'by-gene', sachant qu'un même
        SRR peut-être présent plusieurs fois pour un même fichier)

From the count table of CCLE INDEx (raw counts)

6. `Recurrent kmer (>100)`: kmer-Alt with positive count in >100 cell lines /1019 (to exclude)
        Nombre de variants dont AU MOINS 100 LIGNEES > 0
7. `KmerAlt > 0`: kmer-Alt with positive count
        Nombre de variants dont AU MOINS UNE LIGNEE > 0
8. `KmerAlt =0`: kmer-Alt with negative count
        Nombre de variants dont TOUTES LES LIGNEES == 0
9. `Kmer>0 (1hit)`: kmer-Alt with positive count in one cell line
        Nombre de variants dont UNE SEULE LIGNEE > 0
10. `Kmer Alt>0 (hit =1 OK)`: kmer-Alt with positive count in the expected cell line (depMAP info)
        Nombre de variants dont UNE SEULE LIGNEE > 0 A ÉTÉ TROUVéE dans DEPMAP (pour le même variant)
11. `Kmer alt>0 (hit>1)`: kmer-Alt with positive count in more than one cell line
        Nombre de lignées
12. `Kmer alt>0 (hit>1 OK)`: kmer-Alt with positive count in more than one cell line and expected cell lines (DepMAp info)
        si aucun SRR>0 ne matche avec DepMap pour un kmer>0,
        retrourne le nombre de SRR>0 pour incrémenter un compteur pour le gene.
        A la place, on peut retourner 1 pour donner le nombre de kmer>0 mais dont aucun SRR ne matche
13. `Kmer alt>0 (hit>1 not OK)`: kmer-Alt with positive count in more than one cell line and not in the expected cell line (DepMAp info)
        Lorsque un seul SRR>O est trouvé pour un kmer>0,
        retourne 1 pour incrémenter un compteur pour le gène
14. `CELL_LINE >0`: number of cell lines with a positive counts (to compare with the number of total SRR from fasta file)
15. `Kmer-Ref =0`: number of corresponding kmer-alt + 0 (to compare with the overlap Kmer-ALT number)
"""


import os
import sys
import argparse
import yaml
import subprocess
import pandas as pd


__version__     = "0.2.0"
__author__      = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    args = usage()
    ### open config file and configure some path
    with open("config.yaml") as fh:
        config = yaml.full_load(fh)
    subdir = os.path.join(args.dir, f"a{args.abund_min}-r{args.recur_max}-minpos{args.min_pos}")
    if config["filter"]["enabled"]:
        rdeer_alt_dir = os.path.join(subdir, "filter_kmer/alt") 
    else:
        rdeer_alt_dir = os.path.join(args.dir, f"merge_kmer-minpos{config['merge_kmer']['min_pos']}/alt")
    ### 
    var = Vars(args, subdir)
    a_min = args.abund_min
    ### handle recurrency value when '%' is set
    if args.recur_max.endswith("%"):
        r_max = int(args.recur_max[:-1])
        root,dirs,files = list(os.walk(rdeer_alt_dir))[0]
        file = f"{os.path.splitext(files[0])[0]}.tsv"
        with open(os.path.join(root, file)) as fh:
            header = fh.readline()
            ncol = len(header.split('\t')) -1
        r_max = round(r_max/100*ncol)
    else:
        r_max = int(args.recur_max)


    table = {
        ### SEQTAILOR/DEPMAP
        "Gene": var.genes,                                 # seqtailor-results-exonic (but any results directory is fine)
        "Total kmer": [],                                  # seqtailor-results
        "SRR hit max": [],                                 # seqtailor-results-exonic
        "Distinct SRR": [],                                # seqtailor-results-exonic
        "Total SRR": [],                                   # seqtailor-results-exonic
        ### REINDEER
        "Recurrent kmer (>100)": [],                       # reindeer-alt-results
        "Kmer-alt>0": [],                                  # reindeer-results-alt
        "Kmer-alt=0": [],                                  # reindeer-results-alt
        "Kmer-alt>0 (hit<100)": [],                        # reindeer-results-alt
        "Kmer-alt>0 (hit<100 ok)": [],                     # reindeer-results-alt + seqtailor-results-exonic
        "Kmer-alt>6 (hit<100)": [],                        # reindeer-results-alt
        "Kmer-alt>6 pairs (hit<100)": [],                  # reindeer-results-alt
        f"Kmer-alt>{a_min} (hit<{args.recur_max} ok)": [], # merge_kmer / alt
        "Kmer-alt>0 (ok)": [],                             # merge_kmer / alt
        "Kmer-alt>0 (hit>1 not ok)": [],                   # reindeer-results-alt + seqtailor-results-exonic
        "Kmer-alt>10 (ok)": [],                            # reindeer-results-alt + seqtailor-results-exonic
        "CELL_LINE >0": [],                                # reindeer-results-alt
        "CELL_LINE >10": [],                               # reindeer-results-alt
        "Kmer-ref=0": [],                                  # reindeer-results-ref
        "Kmer-alt>0 (hit=1)": [],                          # reindeer-results-alt
        "Kmer-alt>0 (hit=1 ok)": [],                       # reindeer-results-alt + seqtailor-results-exonic
        "Kmer-alt>0 (hit=1 not ok)": [],                   # reindeer-results-alt + seqtailor-results-exonic
        }


    for i,file in enumerate(var.base_files):

        ### build dataframes for depmap and reindeer
        depmap = depmap_to_dataframe(var, file)
        reindeer_alt = reindeer_to_dataframe(file, dir=rdeer_alt_dir)
        reindeer_ref = reindeer_to_dataframe(file, dir=var.rdeer_ref_dir)
        # ~ reindeer_merged = reindeer_to_dataframe(file, dir=var.rdeer_merged_dir)
        # ~ reindeer_filtered = reindeer_to_dataframe(file, dir=var.rdeer_filtered_dir)

        ### get infos

        ############################### DEPMAP ###############################

        ### depmap total: get total kmer count from output of vcf2seq
        table["Total kmer"].append(get_fasta_nbseq(file, dir=var.vcf2seq_total_dir))

        ### depmap exonic: get largest number of SRRs with the same kmer
        srr_max = 0
        for meta in depmap["meta"]:
            srr_max = max(srr_max, len(meta))
        table["SRR hit max"].append(srr_max)

        ### depmap exonic: get number of distinct SRR in the fasta file
        srr_all = set()
        srr_distinct = 0
        for meta in depmap["meta"]:
            srr_all.update([item.split("|")[0] for item in meta])
        table["Distinct SRR"].append(len(srr_all))

        ### depmap exonic: get total SRR count in the fasta file
        table["Total SRR"].append(int(depmap["nb_meta"].sum()))

        ############################### REINDEER ###############################

        ### Reindeer: recurrent kmer (>100): kmer-Alt with positive count in >100 cell lines /1019 (to exclude)
        table["Recurrent kmer (>100)"].append(reindeer_alt.gt(0).sum(axis=1).gt(100).sum())

        ### Reindeer: KmerAlt > 0: kmer-Alt with at least one positive count
        table["Kmer-alt>0"].append(reindeer_alt.gt(0).sum(axis=1).gt(0).sum())

        ### Reindeer: KmerAlt = 0: kmer-Alt with all negative count
        table["Kmer-alt=0"].append(reindeer_alt.eq(0).sum(axis=1).eq(len(reindeer_alt.columns)).sum())

        ### kmer-Alt with count > n in one cell line
        table["Kmer-alt>0 (hit=1)"].append(reindeer_alt.gt(0).sum(axis=1).eq(1).sum())

        table["Kmer-alt>0 (hit<100)"].append(reindeer_alt.gt(0).sum(axis=1).between(1,100).sum())
        table["Kmer-alt>6 (hit<100)"].append(reindeer_alt.gt(6).sum(axis=1).between(1,100).sum())
        table["Kmer-alt>6 pairs (hit<100)"].append(sum([x for x in reindeer_alt.gt(6).sum(axis=1) if x < 100]))

        ### Number of cell lines with a positive counts (to compare with the number of total SRR from fasta file)
        table["CELL_LINE >0"].append(reindeer_alt.gt(0).sum(axis=1).sum() or 0)
        table["CELL_LINE >10"].append(reindeer_alt.gt(10).sum(axis=1).sum() or 0)

        ### `Kmer-Ref =0`: number of corresponding kmer-alt + 0 (to compare with the overlap Kmer-ALT number)
        table["Kmer-ref=0"].append(reindeer_ref.eq(0).sum(axis=1).eq(len(reindeer_alt.columns)).sum())

        ########################### DEPMAP + REINDEER ##########################

        ### `Kmer Alt>n (OK)`: kmer-Alt with count > n in the expected cell line (depMAP info)
        ### if uniq = True the 'variant > n' is present once for the variant in reindeer file
        table["Kmer-alt>0 (ok)"].append(intersect(depmap, reindeer_alt, gene=var.genes[i]))
        table["Kmer-alt>10 (ok)"].append(intersect(depmap, reindeer_alt, nb_kmer=10))
        table["Kmer-alt>0 (hit=1 ok)"].append(intersect(depmap, reindeer_alt, op_hit="==", nb_hit=1))
        table["Kmer-alt>0 (hit<100 ok)"].append(intersect(depmap, reindeer_alt, op_hit="<",nb_hit=100))
        table[f"Kmer-alt>{a_min} (hit<{args.recur_max} ok)"].append(intersect(depmap, reindeer_alt,
                                            op_kmer="gt", nb_kmer=int(a_min), op_hit="<",nb_hit=int(r_max)))
        table["Kmer-alt>0 (hit>1 not ok)"].append(intersect(depmap, reindeer_alt, op_hit=">",
                                            nb_hit=1, not_ok=True, gene=var.genes[i]))
        table["Kmer-alt>0 (hit=1 not ok)"].append(intersect(depmap, reindeer_alt, op_hit="==",
                                            nb_hit=1, not_ok=True))


    output(args, table)


def depmap_to_dataframe(var, file):
    """
    arguments: output of mytailorisrich
    return: dataframe with columns:
        mut: chr_pos_ref_alt
        meta: [SRR|AC|VAF,...]
        nb: meta element count (SRR with this mutation)
        seq: sequence
    """
    variant = []
    meta = []
    nb_meta = []
    seq = []
    with open(os.path.join(var.vcf2seq_alt_dir, f"{file}.fa")) as fh:
        for line in fh:
            if line.startswith(">"):
                mut, id, smeta, nb = line.rstrip("\n").split(" ")
                mut = "_".join(mut[1:].split("_")[:4])
                lmeta = smeta.split(",")
                variant.append(mut)
                meta.append(lmeta)
                nb_meta.append(int(nb))
            else:
                seq.append(line.rstrip())

    data = pd.DataFrame({
                        "meta": meta,
                        "nb_meta": nb_meta,
                        "sequence": seq,
                        },
                        index=variant)
    return data


def reindeer_to_dataframe(file, dir=None):
    """ Function doc """
    file_path = os.path.join(dir, f"{file}.tsv")
    try:
        data = pd.read_csv(file_path, sep="\t", index_col=0)
        data.index = [ '_'.join(x.split("_")[:4]) for x in data.index ]
    except pd.errors.EmptyDataError:
        data = pd.DataFrame()
    return data


def get_fasta_nbseq(file, dir=None):
    infile = os.path.join(dir, f"{file}.fa")
    cmd = ["wc", "-l", infile]
    nb_line = subprocess.check_output(cmd)
    nb_line = int(nb_line.decode().split(" ")[0]) // 2
    return nb_line


def intersect(depmap, reindeer, op_kmer="gt", nb_kmer=0, op_hit='>', nb_hit=0, not_ok=False, gene=None):
    """
        op_kmer: gt, lt, ge, eq, ...
        op_hit: ==, >, <=, ...
    """

    count = 0
    nocount = 0

    ### reindeer dict from dataframe keys: variants, values: [srr,]]
    # ~ print("reindeer length:", len(reindeer))
    ### filtered by  'kmer > 0' by default
    if reindeer.empty:
        return 0
    rdeer_dic = dict(getattr(reindeer, op_kmer)(nb_kmer).dot(reindeer.columns + '@'))  # idem ~ rdeer_dic = dict(reindeer.gt(nb).dot(reindeer.columns + '@'))
    ### values as list (SRR)
    rdeer_dic = {k:v.rstrip("@").split("@") for k,v in rdeer_dic.items() if v }
    ### filtered by 'hit > 0' by default
    rdeer_dic = {k:v for k,v in rdeer_dic.items() if eval(f"{len(v)}{op_hit}{nb_hit}")}

    ### depmap dict from dataframe keys: variants, values: [srr,]]
    depm_dic = dict(depmap["meta"])
    for key in depm_dic:
        depm_dic[key] = [ a.split("|")[0] for a in depm_dic[key] ]

    ### find matching variants with same srr
    for variant, lsrr in rdeer_dic.items():
        srr_not_in_depm_variant = 0
        try:
            intersect = [value for value in depm_dic[variant] if value in lsrr]
        except KeyError:
            continue
        if intersect:
            count += len(intersect)         # 1 => matching variants counts, len(intersect) => matching srr counts
        else:
            nocount += 1 # len(lsrr)

    if not_ok:
        return nocount
    return count


def output(args, table):
    table = pd.DataFrame(table)

    ### add totals row
    table.loc['total'] = table.sum(numeric_only=True)
    table.loc['total', 'Gene'] = 'Total'

    if args.outfile:
        table.to_csv(args.outfile, sep="\t", index=False)
    else:
        print(table.to_string(index=False))



class Vars:
    """ Class doc """

    def __init__(self, args, subdir):
        """ Class initialiser """
        base_dir = args.dir
        self.vcf2seq_total_dir   = os.path.join(base_dir, "vcf2seq/alt")
        self.vcf2seq_alt_dir     = os.path.join(base_dir, "vcf2seq/alt")
        # ~ self.rdeer_alt_dir    = os.path.join(subdir, "filter_kmer/alt")
        self.rdeer_ref_dir    = os.path.join(base_dir, f"merge_kmer-minpos{args.min_pos}/ref")
        self.base_files       = [os.path.splitext(file)[0] for file in sorted(os.listdir(self.vcf2seq_total_dir)) if not file.startswith('.')]
        self.genes            = [file.split("_")[-1] for file in self.base_files]


def usage():
    doc_sep = '=' * min(72, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("dir",
                        help="directory",
                        metavar=('directory'),
                       )
    parser.add_argument("-r", "--recur-max",
                        default=100,
                       )
    parser.add_argument("-a", "--abund-min",
                        default=6,
                       )
    parser.add_argument("-m", "--min-pos",
                        default=3,
                       )
    parser.add_argument("-o", "--outfile",
                        help="output file",
                        metavar=('outfile.tsv'),
                        nargs="?",
                       )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    return parser.parse_args()


if __name__ == "__main__":
    main()
