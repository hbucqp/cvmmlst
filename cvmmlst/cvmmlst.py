import os
import sys
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from mlst import mlst


def args_parse():
    "Parse the input argument, use '-h' for help."
    parser = argparse.ArgumentParser(
        usage='cvmmlst -i <genome assemble directory> -o <output_directory> \n\nAuthor: Qingpo Cui(SZQ Lab, China Agricultural University)\n')
    parser.add_argument(
        "-i", help="<input_path>: the PATH to the directory of assembled genome files")
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
    subprocess.run("bash mlstdb_setup.sh", shell=True,
                   capture_output=True, encoding='utf-8')


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

        # get the input path
        input_path = os.path.abspath(args.i)

        # check if the output directory exists
        if not os.path.exists(args.o):
            os.mkdir(args.o)

        output_path = os.path.abspath(args.o)

        # get the database path
        database_path = os.path.join(
            os.path.dirname(__file__), f'db/blast/mlst.fa')

        # print(database_path)

        for file in os.listdir(input_path):
            file_base = str(os.path.splitext(file)[0])
            output_filename = file_base + '_tab.txt'
            outfile = os.path.join(output_path, output_filename)
            # print(file_base)
            file_path = os.path.join(input_path, file)
            if os.path.isfile(file_path):
                # print("TRUE")
                if mlst.is_fasta(file_path):
                    print(f'Processing {file}')
                    result = mlst(file_path, database_path, output_path,
                                  threads, minid, mincov).biopython_blast()
                    print(
                        f"Finishing process {file}: writing results to " + str(outfile))
                    sch = mlst.best_scheme(result)

                    df = mlst.get_st(result, sch)

                    df['FILE'] = file_base
                    df.to_csv(outfile, sep='\t', index=False)
                    # change all tab results to pivot table fomat
                    df_all = pd.concat([df_all, df])

                # if args.store_arg_seq:
                #     Blaster.get_arg_seq(file_base, result_dict, output_path)

        # output final pivot dataframe to outpu_path
        summary_file = os.path.join(output_path, 'mlst_summary.csv')
        df_all.to_csv(summary_file, index=False)


if __name__ == '__main__':
    main()
