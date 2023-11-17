#!/usr/bin/env python3

"""
Aims : filter results
- recurrency
- noise

input : tsv file with raw results
output : tsv file with filtered results

"""


import sys
import os
import argparse


__appname__   = "filter"
__shortdesc__ = "Short description."
__licence__   = "none"
__version__   = "0.1.0"
__author__    = "Benoit Guibert <benoit.guibert@free.fr>"


def main():
    """ Function doc """
    args = usage()
    # ~ print(args)
    results = Filter(args)


class Filter:
    """ Class doc """

    def __init__(self, args):
        """ Class initialiser """
        self.args = args
        self.counts = []
        self.report = ""
        self.header = args.infile.readline()
        self.columns = self.header.split('\t')[1:]

        ### build list of variants/counts
        for row in args.infile:
            variant, *counts = row.rstrip('\n').split()
            self.counts.append({"variant":variant, "counts":[float(c) for c in counts]})
        ### recurrency filter
        self.recurrency()

        ### noise filter
        self.noise()

        ### output
        self.write_output()

        ### report
        self.write_report()


    def noise(self):
        reports = []
        ### filtering
        for item in self.counts:
            for i,count in enumerate(item['counts']):
                if 0 < count < self.args.abund_min:
                    item['counts'][i] = 0
                    reports.append(f"{item['variant'].rstrip('.alt')}/{self.columns[i]}: {count:g}")

        ### report
        self.report += f"\nMINIMUM ABUNDANCE  - {len(reports)} found (threshold: {self.args.abund_min})\n\n"
        if reports:
            for report in reports:
                self.report += f"    {report}\n"
        else:
            self.report += "    no noise found\n"


    def recurrency(self):
        recurrents = []
        ### handle recurrency value when '%' is set
        if self.args.recur_max.endswith("%"):
            recur = int(self.args.recur_max[:-1])
            ncol = len(self.header.split()) -1
            recur = round(recur/100*ncol)
        else:
            recur = int(self.args.recur_max)
        ### filtering
        for item in self.counts:
            pos_counts = [c for c in item['counts'] if c > 0]
            if len(pos_counts) > recur:
                item['counts'] = [0 for _ in item['counts']]
                recurrents.append(f"{item['variant'].rstrip('.alt')}: {len(pos_counts)}") # reporting
        ### report
        self.report += f"\nRECURRENCY - {len(recurrents)} found (threshold: {recur})\n\n"
        self.report +=  "  Recurent mutations:\n\n"
        if recurrents:
            for recurrent in recurrents:
                self.report += f"    {recurrent}\n"
        else:
            self.report += "    no recurrent mutation found\n"







    def write_output(self):
        if not self.args.outfile:
            print(self.header, end='')
            for item in self.counts:
                print("{}\t{}".format(item['variant'], '\t'.join([str(f"{i:g}") for i in item['counts']])))
        else:
            self.args.outfile.write(self.header)
            for item in self.counts:
                self.args.outfile.write("{}\t{}\n".format(item['variant'], '\t'.join([str(f"{i:g}") for i in item['counts']])))


    def write_report(self):
        with open(self.args.report, "w") as fic:
            fic.write(f"Gene: {self.args.infile.name.split('_')[-1].rstrip('.tsv')}\n")
            fic.write(self.report)



def usage():
    """
    Help function with argument parser.
    https://docs.python.org/3/howto/argparse.html?highlight=argparse
    """
    doc_sep = '=' * min(72, os.get_terminal_size(2)[0])
    parser = argparse.ArgumentParser(description= f'{doc_sep}{__doc__}{doc_sep}',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,)
    parser.add_argument("infile",
                        help="input tsv file",
                        type=argparse.FileType('r'),
                       )
    parser.add_argument("-o", "--outfile",
                        help="output tsv file",
                        type=argparse.FileType('w'),
                       )
    parser.add_argument('-r', '--recur-max',
                        help="maximum recurrency, as absolute or percentage with %% (default: 50%)",
                        default="50%",
                       )
    parser.add_argument('-a', '--abund-min',
                        type=int,
                        help="minimum abundance (default: 5)",
                        default=5,
                       )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f"{parser.prog} v{__version__}",
                       )
    ### Go to "usage()" without arguments or stdin
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    ### Set default output file name
    if not args.outfile:
        args.outfile = open(f"{os.path.splitext(os.path.basename(args.infile.name))[0]}.filtered.tsv", "w")
    ### Set default report file name
    args.report = f"{os.path.splitext(args.outfile.name)[0]}.log"
    return args


if __name__ == "__main__":
    main()
