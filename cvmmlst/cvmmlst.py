#!/usr/bin/python3

# -*- coding:utf-8 -*-

import os
import sys
import argparse
import subprocess
import pandas as pd
import time
import shutil
from tabulate import tabulate
from Bio import SeqIO
from cvmcore.cvmcore import cfunc
from .mlst import mlst  # remember add dot


def args_parse():
    "Parse the input argument, use '-h' for help."
    parser = argparse.ArgumentParser(
        usage='cvmmlst -i <genome assemble directory> -o <output_directory> \n\nAuthor: Qingpo Cui(SZQ Lab, China Agricultural University)\n')

    # Add subcommand
    subparsers = parser.add_subparsers(
        dest="subcommand", title="cvmmlst subcommand")
    init_db_parser = subparsers.add_parser(
        'init', help='<initialize the reference database>')
    # init_db_parser.set_defaults(func=initialize_db)
    # init_db_parser.add_argument('-init', help= '<initialize the reference database>')

    show_schemes_parser = subparsers.add_parser(
        'show_schemes', help="<show the list of all available schemes>")
    # show_schemes_parser.set_defaults(func=show_db_list)
    # show_schemes_parser.add_argument('-show_schemes', help="<show the list of schemes>")

    add_scheme_parser = subparsers.add_parser(
        'add_scheme', help='<add custome scheme, use cvmmlst add_scheme -h for help>')
    add_scheme_parser.add_argument('-name', help="<the custome scheme name>")
    add_scheme_parser.add_argument('-path', help='<the path of scheme>')

    # Add mutually exclusively parameters
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "-i", help="<input_path>: the PATH to the directory of assembled genome files. Could not use with -f")
    group.add_argument(
        "-f", help="<input_file>: the PATH of assembled genome file. Could not use with -i")

    # Add options
    parser.add_argument("-o", help="<output_directory>: output PATH")
    parser.add_argument(
        "-scheme", help="<mlst scheme want to use>, cvmmlst show_schemes command could output all available schems")
    parser.add_argument('-minid', default=90,
                        help="<minimum threshold of identity>, default=90")
    parser.add_argument('-mincov', default=60,
                        help="<minimum threshold of coverage>, default=60")
    parser.add_argument(
        '-t', default=8, help='<number of threads>: default=8')
    parser.add_argument('-v', '--version', action='version',
                        version='Version: ' + get_version("__init__.py"), help='Display version')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args()


def get_rel_path():
    """
    Get the relative path
    """
    here = os.path.abspath(os.path.dirname(__file__))
    return here


def read(rel_path: str) -> str:
    # here = os.path.abspth(os.path.dirname(__file__))
    here = get_rel_path()
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
    """
    Initialize mlst database
    """
    DIR = get_rel_path()
    MLSTDIR = os.path.join(DIR, 'db', 'pubmlst')
    BLASTDIR = os.path.join(DIR, 'db', 'blast')
    BLASTFILE = os.path.join(BLASTDIR, 'mlst.fa')
    # print(BLASTFILE)

    os.makedirs(BLASTDIR, exist_ok=True)

    if os.path.exists(BLASTFILE):
        os.remove(BLASTFILE)

    for root, dirs, files in os.walk(MLSTDIR):
        for scheme in dirs:
            scheme_path = os.path.join(MLSTDIR, scheme)
            scheme_blastfile = os.path.join(scheme_path, 'mlst.fa')
            print('-' * 30)
            print(f"Adding: {scheme}")
            if os.path.exists(scheme_blastfile):
                os.remove(scheme_blastfile)
            for file in os.listdir(scheme_path):
                if file.endswith('.tfa'):
                    with open(os.path.join(scheme_path, file), 'r') as f:
                        for line in f:
                            if 'not a locus' not in line:
                                if line.startswith('>'):
                                    line = f">{scheme}.{line[1:]}"
                                with open(BLASTFILE, 'a') as blast_file:
                                    blast_file.write(line)
                                with open(scheme_blastfile, 'a') as scheme_blast_file:
                                    scheme_blast_file.write(line)
            # print(scheme_blastfile)
            print(f"Created scheme specific BLAST database for {scheme}")
            subprocess.run(['makeblastdb', '-hash_index', '-in', scheme_blastfile, '-dbtype',
                            'nucl', '-title', 'PubMLST', '-parse_seqids'], shell=True, text=True)
            # print('------------------------------')

    subprocess.run(['makeblastdb', '-hash_index', '-in', BLASTFILE, '-dbtype',
                    'nucl', '-title', 'PubMLST', '-parse_seqids'], shell=True, text=True)
    print('-' * 30)
    print(
        f"Created merged BLAST database for all availbale schemes using {BLASTFILE}")


def add_scheme(name: str, path: str):
    """
    Add user defined scheme to database
    The name is the scheme name you want to sue when running with cvmmlst -scheme {name}
    The path is the fasta files and profile download from pubmlst or user defined
    """
    # name = 'custome'
    # path = './custome'

    path = os.path.abspath(path)
    ref_path = get_rel_path()
    dest_path = os.path.join(ref_path, 'db/pubmlst', name)

    if not os.path.exists(dest_path):
        os.mkdir(dest_path)

    extensions = ('.fasta', '.fsa', 'fa', 'fna', 'tfa')

    for file in os.listdir(path):
        file = os.path.join(path, file)
        if os.path.isfile(file):
            file_basename = os.path.splitext(os.path.basename(file))[0]
            if file.endswith(extensions):
                if mlst.is_fasta(file):
                    dest_file = os.path.join(dest_path, f'{file_basename}.tfa')
                    # print(dest_file)
                    shutil.copy(file, dest_file)
                else:
                    print(
                        f'Wrong fasta format of {file} \n, please check your input files in {path}')
                    sys.exit(1)
            elif file.endswith('.txt'):
                dest_file = os.path.join(dest_path, f'{name}.txt')
                # print(dest_file)
                shutil.copy(file, dest_file)


def show_db_list():
    """
    Convert the mlst database scheme to tidy dataframe
    Paramters
    ----------
    # rel_path:
    #     database path
    Returns
    ----------
    A tidy dataframe contains the scheme_name, No. of STs, No. of Locus and the last modified date
    """
    # number_locus = []
    # number_STs = []
    # update_date = []
    # print(args)
    here = get_rel_path()
    # print(here)
    rel_path = os.path.join(here, 'db/pubmlst')
    schemes_list = []
    schemes = os.listdir(rel_path)
    # print(schemes)
    for file in schemes:
        # print(file)
        files_path = os.path.join(rel_path, file)
        if os.path.isdir(files_path):
            # print(files_path)
            scheme_dict = {}
            scheme_path = os.path.join(rel_path, files_path)
            loci_num = 0
            scheme_lmd_time = []
            for scheme_file in os.listdir(scheme_path):
                # print(scheme_file)
                scheme_file = os.path.join(scheme_path, scheme_file)
                if scheme_file.endswith('.tfa'):
                    loci_num += 1
                    # Get the last modified time
                    # First tfa files modified time
                    scheme_update_date = cfunc.get_mod_time(scheme_file)
                if scheme_file.endswith('.txt'):
                    df = pd.read_csv(scheme_file, sep='\t')
                    STs_number = df.shape[0]
            # number_locus.append(loci_num)
            # update_date.append(scheme_lmd_time[0])
            scheme_dict['Schemes'] = file
            scheme_dict['No. of STs'] = STs_number
            scheme_dict['No. of Locus'] = loci_num
            scheme_dict['Update_date'] = scheme_update_date
            schemes_list.append(scheme_dict)
        else:
            next
    db_df = pd.DataFrame(schemes_list)
    # print(db_df)
    tidy_schemes_df = tabulate(db_df, headers='keys')
    return print(tidy_schemes_df)


def main():
    args = args_parse()
    if args.subcommand is None:
        # Parser options
        df_all = pd.DataFrame()  # output summary dataframe
        threads = args.t
        minid = args.minid
        mincov = args.mincov

        # Check if the output directory exists
        if args.o is not None:
            if not os.path.exists(args.o):
                os.mkdir(args.o)
            output_path = os.path.abspath(args.o)
        # print(output_path)

        # Get the corresponding database path
        if args.scheme is not None:
            database_path = os.path.join(os.path.abspath(
                get_rel_path()), 'db/pubmlst', args.scheme, 'mlst.fa')
            if not os.path.exists(database_path):
                print(
                    f"Your input {args.scheme} scheme doesn't exist, please check available schemes using <cvmmlst show_schemes> command.")
                # print(database_path)
                sys.exit(1)
        else:
            # database_path = os.path.join(os.path.dirname(__file__), os.path.join(os.path.abspath(get_rel_path()), 'db/blast/mlst.fa'))
            database_path = os.path.join(os.path.join(
                os.path.abspath(get_rel_path()), 'db/blast/mlst.fa'))
            # print(database_path)

        files = []

        if args.i is not None:
            # get the input path
            files = os.listdir(os.path.abspath(args.i))
            input_path = os.path.abspath(args.i)
        elif args.f is not None:
            files.append(os.path.abspath(args.f))
            input_path = os.path.dirname(os.path.abspath(args.f))

        # Run mlst
        for file in files:
            file_base = str(os.path.basename(os.path.splitext(file)[0]))
            output_filename = file_base + '_tab.txt'
            # print(output_path)
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
    elif args.subcommand == 'show_schemes':
        # print('subcommand test')
        show_db_list()
        # sys.exit(1)
    elif args.subcommand == 'init':
        initialize_db()
        # schemes_dir = os.path.join(get_rel_path(), 'db/pubmlst')
        # initialize_schemes(schemes_dir)
        # sys.exit(1)
    elif args.subcommand == 'add_scheme':
        # print(args.name)
        # print(args.path)
        add_scheme(args.name, args.path)
        print(f'Adding {args.name} to reference database...')
        print(f'Initializing reference data...')
        initialize_db()
    else:
        print(
            f'{args.subcommand} do not exists, please using "cvmmlst -h" to show help massage.')


if __name__ == '__main__':
    main()
