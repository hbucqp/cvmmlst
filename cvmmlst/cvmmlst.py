#!/usr/bin/python3

# -*- coding:utf-8 -*-

import os
import sys
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from .mlst import mlst


def args_parse():
    "Parse the input argument, use '-h' for help."
    parser = argparse.ArgumentParser(
        usage='cvmmlst -i <genome assemble directory> -o <output_directory> \n\nAuthor: Qingpo Cui(SZQ Lab, China Agricultural University)\n')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-i", help="<input_path>: the PATH to the directory of assembled genome files. Could not use with -f")
    group.add_argument(
        "-f", help="<input_file>: the PATH of assembled genome file. Could not use with -i")
    parser.add_argument("-o", help="<output_directory>: output PATH")
    parser.add_argument('-minid', default=90,
                        help="<minimum threshold of identity>, default=90")
    parser.add_argument('-mincov', default=60,
                        help="<minimum threshold of coverage>, default=60")
    parser.add_argument('-init', action='store_true',
                        help='<initialize the reference database>')
    parser.add_argument(
        '-t', default=8, help='<number of threads>: default=8')
    parser.add_argument('-v', '--version', action='version',
                        version='Version: ' + get_version("__init__.py"), help='Display version')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def read(rel_path: str) -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path)) as fp:
        return fp.read()


def get_version(rel_path: str) -> str:
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


def initialize_db():
    print("Creating mlst blast database...")
    subprocess.run("bash mlstdb_setup.sh", shell=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, encoding='utf-8')
    print("Done")


def main():
    df_all = pd.DataFrame()
    args = args_parse()
    if args.init:
        initialize_db()
    # elif args.init:
    #     initialize_db()
    else:
        # threads
        threads = args.t
        # print(threads)

        minid = args.minid
        mincov = args.mincov

        # check if the output directory exists
        if not os.path.exists(args.o):
            os.mkdir(args.o)
        output_path = os.path.abspath(args.o)
        # print(output_path)

        # get the database path
        database_path = os.path.join(
            os.path.dirname(__file__), os.path.join(os.path.abspath(os.path.dirname(__file__)), 'db/blast/mlst.fa'))

        files = []

        if args.i is not None:
            # get the input path
            files = os.listdir(os.path.abspath(args.i))
            input_path = os.path.abspath(args.i)
        else:
            files.append(os.path.abspath(args.f))
            input_path = os.path.dirname(os.path.abspath(args.f))

        for file in files:
            file_base = str(os.path.basename(os.path.splitext(file)[0]))
            output_filename = file_base + '_tab.txt'
            # print(output_path)
            # print('xxx')
            # print(file_base)
            outfile = os.path.join(output_path, output_filename)
            # print(outfile)
            file_path = os.path.join(input_path, file)
            if os.path.isfile(file_path):
                # print("TRUE")
                if mlst.is_fasta(file_path):
                    print(f'Processing {file}')
                    result = mlst(file_path, database_path, output_path,
                                  threads, minid, mincov).biopython_blast()
                    # print(result) # for debug
                    if len(result) != 0:
                        # sch = mlst.best_scheme(result)
                        df = mlst.get_st(result)
                        if len(df) != 0:
                            df['FILE'] = file_base
                            order = list(reversed(df.columns.to_list()))
                            df = df[order]
                            # print(df)
                            df.to_csv(outfile, sep='\t', index=False)
                            print(
                                f"Finishing process {file}: writing results to " + str(outfile))
                        # else:
                        #     df = pd.DataFrame()
                        #     df['Note'] = 'Could not matching any loci in all schemes, next...'
                        #     df['ST'] = '-'
                        #     df['Scheme'] = '-'
                        #     df['FILE'] = file_base
                        #     print(
                        #         f'Could not found similar scheme of {file_base}, writing result to ' + str(outfile))
                    else:
                        df = pd.DataFrame()
                        df['Note'] = 'Could not matching any loci in all schemes, next...'
                        df['ST'] = '-'
                        df['Scheme'] = '-'
                        df['FILE'] = file_base
                        print(
                            f'Could not matching any loci in all schemes, next...')
                        order = list(reversed(df.columns.to_list()))
                        df = df[order]
                        df.to_csv(outfile, sep='\t', index=False)

                    df_all = pd.concat([df_all, df])

        # output final pivot dataframe to outpu_path
        summary_file = os.path.join(output_path, 'mlst_summary.csv')
        df_all.to_csv(summary_file, index=False)


if __name__ == '__main__':
    main()
