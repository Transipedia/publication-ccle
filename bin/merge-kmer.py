#!/usr/bin/env python3

"""
Merge 61 splited kmers by mutation
input : output of kmerator
output: single tsv file with average of positive counts
"""

import sys
import os
import argparse
from statistics import mean


__appname__   = "application_name"
__shortdesc__ = "Short description."
__licence__   = "none"
__version__   = "0.1.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    """ Function doc """
    args = usage()
    # ~ header, reindeer_dict = mk_reindeer_dict(args)
    # ~ print(reindeer_dict.keys())
    # ~ print(len(reindeer_dict))
    header, pos_counts = mean_pos_counts(args)
    output(args, header, pos_counts)


def mean_pos_counts(args):
    ### reporting
    # ~ dir = f"{os.path.dirname(args.output) or '.'}/report"
    # ~ os.makedirs(dir, exist_ok=True)
    # ~ report = []

    ### step 1/3: load rdeer tsv file as dict
    reindeer_dict = {}
    header = args.rdeer_file.readline().rstrip('\n').split('\t')[1:]
    for raw in args.rdeer_file:
        name, *counts = raw.rstrip('\n').split('\t')
        name = name.split('.kmer')[0]
        counts = [int(i) for i in counts]
        ### keep only raw with at list one positive count
        # ~ if sum(counts) > 0:
        reindeer_dict.setdefault(name, []).append(counts)
    ### 2/3: compute average of positives counts
    pos_counts = []
    for k,v in reindeer_dict.items():
        pos_counts_by_srr = {}
        for counts in v:
            for i,count in enumerate(counts):
                if count > 0:
                    pos_counts_by_srr.setdefault(header[i], []).append(count)
        pos_counts.append({k: pos_counts_by_srr})
    ### Rebuilds counts, then insert positives counts and filter raw with positives counts
    rebuilded_rdeer_format = []
    for name_dict in pos_counts:
        for name, SRRs in name_dict.items():
            counts = [ '0' for _ in header]
            for SRR, pos_counts in SRRs.items():
                if len(pos_counts) >= args.min_pos:
                    idx = header.index(SRR)
                    counts[idx] = str(mean(pos_counts))
                # ~ else: print("LOST:", name, SRR, pos_counts, len(pos_counts))
            counts.insert(0,name)
            rebuilded_rdeer_format.append(counts)

    return header, rebuilded_rdeer_format


def output(args, header, pos_counts):
    ###     print rdeer counts (from average of positive counts)
    args.outfile.write("seq_name\t{}\n".format('\t'.join(header)))
    for kmer_counts in pos_counts:
        raw = "{}\n".format("\t".join(kmer_counts))
        args.outfile.write(raw)



def usage():
    """
    Help function with argument parser.
    https://docs.python.org/3/howto/argparse.html?highlight=argparse
    """
    doc_sep = '=' * min(72, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    ### OPTION
    parser.add_argument("rdeer_file",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-m", "--min-pos",
                        type=int,
                        help="minimum positive counts",
                        default=1,
                       )
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType('w'),
                        required=True,
                       )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    return parser.parse_args()

if __name__ == "__main__":
    main()
