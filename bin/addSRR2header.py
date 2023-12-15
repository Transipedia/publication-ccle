#!/usr/bin/env python3

"""
STEP 4 of Snakefile-v4
rule vcf2seq: it add to fasta file information of SRR and VAF
"""

import sys
import os
import argparse

__appname__   = "addSRR2header"
__shortdesc__ = "part of Snakefile-v4: Add SSR and VAF information on fasta headers."
__licence__   = "none"
__version__   = "0.1.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


COLUMNS = ("SRR", "RNAseq_AC", "AF_%")


def main():
    """ Function doc """
    args = usage()

    ### Determine columns to check
    depmap_header = args.depmap.readline().rstrip().split('\t')
    column_ids = []
    for col in COLUMNS:
        column_ids.append(depmap_header.index(col))

    ### build dict with COLUMS of interest
    depmap_dict = {}
    for row in args.depmap:
        row = row.rstrip('\n').split('\t')
        ch, pos, ref, alt = row[0], row[1], row[3], row[4]
        key = f">{ch}_{pos}_{ref}_{alt}"
        add = [row[i] for i in column_ids]
        add = f"{add[0]}|{add[1] or 'unknown'}|{add[2] or 'unknown'}"
        depmap_dict.setdefault(key, []).append(add)


    ### file as data
    data_new = []
    with open(args.fasta) as fh:
        data = fh.readlines()
    for i,row in enumerate(data):
        if row.startswith('>'):
            _row = '_'.join(row.split('_')[:4])
            if _row in depmap_dict:
                data_new.append("{} {} {}\n".format(
                                                row.rstrip('\n'),
                                                ','.join(depmap_dict[_row]),
                                                len(depmap_dict[_row]))
                )
        else:
            data_new.append(row)

    ### Erase existing file by new
    with open(args.fasta, "w") as fh:
        for a in data_new:
            fh.write(a)



def usage():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta",
                       )
    parser.add_argument("-d", "--depmap",
                        type=argparse.FileType('r'),
                        required=True,
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
