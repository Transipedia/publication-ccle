#!/usr/bin/env python3

"""
Like seqtailor, give a VCF file, it return genomic sequence (default length: 31)

Nota:
- When a insertion is larger than '--size' option, only first '--size' nucleotides are outputed.
- header ID are formated like "<chr>_<position>_<ref>_<alt>", and are UNIQUE in the result file.

Specifications: https://github.com/samtools/hts-specs/blob/master/VCFv4.4.pdf
"""

import sys
import os
import argparse
import ascii
try:
    import pyfaidx
except ModuleNotFoundError as err:
    sys.exit(err)


__appname__   = "vcf2seq"
__shortdesc__ = "like seqTailor, give a VCF file, it return genomic sequence."
__licence__   = "GPL3"
__version__   = "0.2.9"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    """ Function doc """
    args = usage()
    try:
        chr_dict = pyfaidx.Fasta(args.genome) # if fai file doesn't exists, it will be automatically created
    except pyfaidx.FastaNotFoundError as err:
        sys.exit(f"FastaNotFoundError: {err}")
    except OSError as err:
        sys.exit(f"\n{COL.RED}WriteError: directory {os.path.dirname(args.genome)!r} may not be "
                  "writable.\nIf you can't change the rights, you can create a symlink and target "
                  f"it. For example:\n  ln -s {args.genome} $HOME\n{COL.END}")
    vcf_ctrl(args, chr_dict)
    res_ref, res_alt = compute(args, chr_dict)
    write(args, res_ref, res_alt)


def vcf_ctrl(args, chr_dict):
    with open(args.vcf.name) as fh:
        for row in fh:
            if row.startswith('#'):
                continue
            chr, pos, id, ref, alt, *rest = row.rstrip('\n').split('\t')

            ### Check some commonly issues
            if not pos.isdigit():
                sys.exit(f"{COL.RED}ErrorVcfFormat: second column is the position. It must be a "
                         f"digit (found: {pos!r}).\n"
                          "A commonly issue is that headers are not commented by a '#' ")
            if chr not in chr_dict:
                sys.exit(f"{COL.RED}ErrorChr: Chromosomes are not named in the same way in the "
                          "query and the genome file. Below the first chromosome found: \n"
                         f" your query: {chr}\n"
                         f" genome: {next(iter(chr_dict.keys()))}\n"
                         f"Please, correct your request (or modify the file '{args.genome}.fai').")
            break


def compute(args, chr_dict):
    res_ref = {}
    res_alt = {}
    num_row = 0
    cols_id = ascii.get_index(args.add_columns)    # columns chars are converted as index, ex: AA -> 27
    for variant in args.vcf:
        num_row += 1
        if not variant.rstrip('\n') or variant.startswith('#'):
            continue
        fields = variant.rstrip('\n').split('\t')
        chr, position, id, ref, alts = fields[:5]

        ### check if --add-columns is compatible with number of columns
        if args.add_columns and max(cols_id) > len(fields):
            sys.exit(f"\n{COL.RED}Error: vcf file has {len(fields)} columns, but you asked for "
                  f"{max(args.add_columns)}.{COL.END}\n")

        alts = alts.split(',')
        for alt in alts:

            header = f">{chr}_{position}_{ref}_{alt}"
            pos = int(position)-1

            ### event bigger than kmer size
            if max(len(ref), len(alt)) > args.size :
                print(f"{COL.PURPLE}Warning: large alteration ({chr}_{pos}_{ref}_{alt}), truncated in output.")


            ### insertion case
            if ref == args.blank:
                # ~ print("-------\n", chr, pos, ref, alt)
                corr = 0 if len(alt) & 1 else 1  # correction when ref is pair

                ref_offset = args.size // 2
                ref_start = pos -ref_offset - 1
                ref_end = pos + ref_offset
                ref_seq = chr_dict[chr][ref_start:ref_end]

                alt_offset = (args.size - len(alt)) // 2
                alt_start = pos - alt_offset - corr
                alt_end = pos + alt_offset
                part1 = chr_dict[chr][alt_start:pos]
                part3 = chr_dict[chr][pos:alt_end]
                alt_seq = f"{part1}{alt}{part3}"

                alt_seq = alt_seq[:args.size]   # reduce alt seq if insertion > --size

                # ~ print("corr : ", corr)
                # ~ print("offset_ref :",ref_offset)
                # ~ print("offset_alt :",alt_offset)
                # ~ print(ref_seq)
                # ~ print(part1, alt, part3)

            ### deletion case
            elif alt == args.blank:
                # ~ print("-------\n", chr, pos, ref, alt)
                corr = 0 if len(ref) & 1 else 1  # correction when ref is unpair

                ref_offset = (args.size - len(ref)) // 2
                ref_start = pos -ref_offset - corr
                ref_end = pos + ref_offset + len(ref)
                ref_seq = chr_dict[chr][ref_start:ref_end]

                alt_offset = args.size // 2
                alt_start = pos - alt_offset - corr
                alt_end = pos + len(ref) + alt_offset - corr + 1
                part1 = chr_dict[chr][alt_start:pos]
                part3 = chr_dict[chr][pos+len(ref):alt_end]
                alt_seq = f"{part1}{part3}"


            ### indel without blank character (first pos is common to ref and alt)
            elif len(alt) != len(ref):

                ### deletion or insertion
                event = "ins" if len(alt) > len(ref) else "del"

                ### It's à little bit tricky for insertion, deletion and paired or not difference
                if not (len(alt)-len(ref)) & 1:
                    corr_diff -= 1
                if event == "ins":
                    corr_diff = -1
                    if not (len(alt)-len(ref)) & 1:
                        corr_diff -= 1
                else:
                    corr_diff = len(ref) - len(alt)

                ### find REF sequence
                corr_ref = 0 if len(ref) & 1 else 1  # correction when ref is pair
                ref_offset = (args.size - len(ref)) // 2
                ref_start = pos -ref_offset - corr_ref
                ref_end = pos + ref_offset + len(ref)
                ref_seq = chr_dict[chr][ref_start:ref_end]

                ### find ALT sequence
                corr_alt = 0 if len(alt) & 1 else 1  # correction when ref is pair
                alt_offset = (args.size - len(alt)) // 2
                alt_start = pos - alt_offset - corr_alt
                alt_end = pos + len(alt) + alt_offset + corr_diff
                part1 = chr_dict[chr][alt_start:pos]
                part3 = chr_dict[chr][pos+len(alt)+corr_diff:alt_end]
                alt_seq = f"{part1}{alt}{part3}"


            ### mutation case
            else:
                NUC = ["A", "T", "C", "G", args.blank]
                for nuc in (ref[0], alt[0]):
                    if nuc not in NUC:
                        sys.exit(f"{COL.RED}Error: The character {nuc!r} is not handle "
                                f"(line {num_row}).\n       You might add the '--blank {nuc}' "
                                "option or check your VCF file.")

                # ~ print("-------\n", chr, pos, ref, alt)
                corr = 0 if len(ref) & 1 else 1  # correction when ref is pair

                ref_offset = (args.size - len(ref)) // 2
                ref_start = pos -ref_offset - corr
                ref_end = pos + ref_offset + len(ref)
                ref_seq = chr_dict[chr][ref_start:ref_end]

                alt_start = pos - ref_offset - corr
                alt_end = pos + len(ref) + ref_offset
                part1 = chr_dict[chr][alt_start:pos]
                part3 = chr_dict[chr][pos+len(ref):alt_end]
                alt_seq = f"{part1}{alt}{part3}"

                ### REF nucleotide must the same than calculated position
                gen_nuc = chr_dict[chr][pos:pos+len(ref)]
                if not ref == gen_nuc:
                    print(f"{COL.PURPLE}Warning: mismatch between REF and genome "
                             f"at line {num_row} ({chr}:{pos+1}).\n"
                             f"          - REF on the vcf file: {ref!r}\n"
                             f"          - found on the genome: '{gen_nuc}'\n"
                              "Please check if the given genome is appropriate.\n")

            col_txt =  ' '.join([fields[num-1] for num in cols_id])
            if len(ref_seq) == args.size == len(alt_seq):
                res_ref[header] = f"{header}_ref {col_txt}\n{ref_seq}"
                res_alt[header] = f"{header}_alt {col_txt}\n{alt_seq}"
            else:
                print(f"{COL.PURPLE}Warning: sequence size not correct at line {num_row}, ignored.{COL.END}")


    return res_ref, res_alt


def write(args, res_ref, res_alt):
    """ Function doc """
    ### define output file
    if not args.output:
        name, ext = os.path.splitext(os.path.basename(args.vcf.name))
        args.output = f"{name}-vcf2seq.fa"
    ### write results in file
    with open(args.output, 'w') as fh:
        if not res_ref:
            return
        if args.type == 'alt':
            fh.write('\n'.join([a for a in res_alt.values()]) + "\n")
        elif args.type == 'ref':
            fh.write('\n'.join([a for a in res_ref.values()]) + "\n")
        elif args.type == 'both':
            res_ref = sorted(list(res_ref.values()))
            res_alt = sorted(list(res_alt.values()))
            for i, _ in enumerate(res_alt):
                fh.write(res_ref[i] + '\n')
                fh.write(res_alt[i] + '\n')



class COL:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def usage():
    doc_sep = '=' * min(80, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("vcf",
                        help="vcf file",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-g", "--genome",
                        help="genome as fasta file",
                        metavar="genome",
                        required=True,
                       )
    parser.add_argument('-s', '--size',
                        type=int,
                        help="size of the output sequence",
                        default=31,
                       )
    parser.add_argument("-t", "--type",
                        type=str,
                        choices=['alt', 'ref', 'both'],
                        default='alt',
                        )
    parser.add_argument("-b", "--blank",
                        type=str,
                        help="Missing nucleotide character",
                        default='-',
                        )
    parser.add_argument("-a", "--add-columns",
                        # ~ type=int,
                        help="Add one or more columns to header (ex: -a 3 AA' will add columns "
                             "3 and 27).",
                        nargs= '+',
                        )
    parser.add_argument("-o", "--output",
                        type=str,
                        help=f"Output file (default: <input_file>-{__appname__}.tsv)",
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
