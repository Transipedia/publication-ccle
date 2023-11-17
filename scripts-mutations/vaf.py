#!/usr/bin/env python3

"""
Compute Reindeer Allelic Frequency and compare with DepMap
inputs:
    - output of reindeer on REF
    - output of reindeer on ALT after rescuing null counts
    - output of DepMap
output:
    - tsv file
"""

import sys
import os
import argparse
import collections


__appname__   = "vaf"
__shortdesc__ = "Compute Reindeer Allelic Frequency and compare with DepMap."
__licence__   = "GPL3"
__version__   = "0.1.1"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"



def main():
    args = usage()
    results = compute(args)
    output(args, results)


def compute(args):
    samples = []
    results = []
    results.append("\t".join(["gene", "mutation", "sample", "rdeer_AF", "rdeer_VAF", "depmap_AF", "depmap_VAF", "rdeer_alt_AF/depmap_alt_AF", "rdeer_VAF/depmap_VAF"]))

    ### Analyze depmap file and load as dict
    depmap = collections.defaultdict(dict)
    for row in args.depmap_file:
        items = row.rstrip('\n').split('\t')
        mut_srr = f"{items[0]}_{items[1]}_{items[3]}_{items[4]}_{items[8]}"
        gene = items[7]
        af = items[21]
        if not af:
            af = vaf = "unknown"
        else:
            al, re = [int(a) for a in af.split(':')]
            vaf = round(al / (al+re) * 100)
        depmap[mut_srr] = {"gene": gene, "af": af, "vaf": vaf}


    ### Load ALT and REF file as dict
    samples = args.reindeer_ref.readline().rstrip("\n").split("\t")[1:]
    args.reindeer_alt.readline().rstrip("\n").split("\t")[1:]
    ref_dict = {}
    alt_dict = {}

    for row in args.reindeer_ref:
        id, *counts = row.rstrip().split('\t')
        id = '_'.join(id.split('_')[:4])
        ref_dict[id] = counts

    for row in args.reindeer_alt:
        id, *counts = row.rstrip().split('\t')
        id = '_'.join(id.split('_')[:4])
        alt_dict[id] = counts

    # ~ print(next(iter(ref_dict)))
    # ~ print(next(iter(alt_dict)))
    # ~ print("nb ref:" ,len(ref_dict))
    # ~ print("nb alt:" ,len(alt_dict))

    ### For each mutation of ALT dict
    for mut, counts_alt in alt_dict.items():
        ### get counts in REF dict
        try:
            counts_ref = ref_dict[mut]
        except KeyError as err:
            # ~ sys.exit("\nError: mutation ALT not found in REF:\n"
                     # ~ f"  gene: {gene}\n  mutation: {mut}\n")
            continue

        ### Compute VAF for each sample for each mutation
        counts_alt = [float(i) for i in counts_alt]
        counts_ref = [float(i) for i in counts_ref]
        for i, sample in enumerate(samples):
            try:
                vaf = round(counts_alt[i] / (counts_ref[i] + counts_alt[i]) * 100)
            except ZeroDivisionError:
                vaf = 100 if counts_alt[i] else 0

            ### create new row with firsts values
            item = f"{gene}\t{mut}\t{sample}\t{round(counts_alt[i])}:{round(counts_ref[i])}\t{vaf}"

            ### Add depmap vaf to item
            dep_key = f"{mut}_{sample}"

            if dep_key in depmap:
                # ~ print("depkey is in depmap:", dep_key)
                item += f"\t{depmap[dep_key]['af']}\t{depmap[dep_key]['vaf']}"
                ### Add some ratio between Reindeer and DepMap results
                if depmap[dep_key]['vaf'] != "unknown":
                    rei_alt_af = float(counts_alt[i])
                    dep_alt_af = int(depmap[dep_key]['af'].split(':')[0])
                    alt2alt = round(rei_alt_af / dep_alt_af, 2)
                    vaf2vaf = round(vaf / depmap[dep_key]['vaf'], 2)
                    item += f"\t{alt2alt}\t{vaf2vaf}"
                else:
                    item += f"\t_\t_"
            else:
                # ~ print("depkey is NOT in depmap:", dep_key)
                item += f"\t_\t_\t_\t_"

            ### filter  depmap_AF == '_' AND (rdeer_AF == '0:0' OR rdeer_AF.split(':')[1] == 0)
            # ~ gene, mut, sample, rd_af, rd_vaf, dm_af, *rest = item.split('\t')
            # ~ if dm_af == '_' and (rd_af == '0:0' or rd_af.split(':')[0] == '0'):
                # ~ continue

            results.append(item)

    return results


def output(args, results):
    if not args.outfile:
        print(*results, sep='\n')
    else:
        for row in results:
            args.outfile.write(f"{row}\n")



def usage():
    """
    Help function with argument parser.
    https://docs.python.org/3/howto/argparse.html?highlight=argparse
    """
    doc_sep = '=' * min(72, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("reindeer_ref",
                        help="Reindeer on REF output",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("reindeer_alt",
                        help="Reindeer on ALT output",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("depmap_file",
                        help="depmap file",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("--order-by",
                        type=str,
                        choices=["gene", "mutation", "sample"],
                        default="gene"
                       )
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType('w'),
                       )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1 and sys.stdin.isatty():
        parser.print_help()
        sys.exit()
    return parser.parse_args()


if __name__ == "__main__":
    main()
