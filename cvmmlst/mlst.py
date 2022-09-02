import os
import re
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

# from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline


class mlst():
    def __init__(self, inputfile, database, output, threads, minid=90, mincov=60):
        self.inputfile = os.path.abspath(inputfile)
        self.database = database
        self.minid = int(minid)
        self.mincov = int(mincov)
        self.temp_output = os.path.join(os.path.abspath(output), 'temp.txt')
        self.threads = threads

    def biopython_blast(self):
        cline = NcbiblastnCommandline(query=self.inputfile, db=self.database, dust='no', ungapped=True,
                                      evalue=1E-20, out=self.temp_output, culling_limit=1,
                                      outfmt="6 sseqid slen length nident",
                                      perc_identity=self.minid, max_target_seqs=10000,
                                      num_threads=self.threads)
        stdout, stderr = cline()
        df = pd.read_csv(self.temp_output, sep='\t', names=[
            'sseqid', 'slen', 'length', 'nident'])

        result = {}
        for i, row in df.iterrows():
            sch, gene, num = re.match(
                '^(\w+)\.(\w+)[_-](\d+)', row['sseqid']).group(1, 2, 3)
            hlen = row['slen']
            alen = row['length']
            nident = row['nident']
            if nident * 100 / hlen >= self.mincov:
                if sch not in result.keys():  # check if sch is the key of result
                    result[sch] = {}
                if hlen == alen & nident == hlen:
                    if gene not in result[sch].keys():
                        result[sch][gene] = num
                    else:
                        if num <= result[sch][gene]:
                            result[sch][gene] = num
                        else:
                            next
                elif alen == hlen:
                    result[sch][gene] = f'~{num}'
                    # result[sch] = mlst
                else:
                    result[sch][gene] = f'{num}?'
        os.remove(self.temp_output)
        return result

    @staticmethod
    def build_genotype(scheme):
        """
        get the corresponding allels and genotype frofiles in the following format:
        col = ['lociA', 'lociB','lociC','lociD','lociE','lociF','lociG']

        {scheme:{'nloci':7,
                 'profiles':{profile:ST}
                 }
        }
        """
        count = 0
        col = []
        genotype = {}
        db_path = os.path.join('db/pubmlst', scheme)
        for file in os.listdir(db_path):
            # print(file)
            if file.endswith('.tfa'):
                base = os.path.splitext(file)[0]
                col.append(base)
                # print(base)
                count += 1
                # if file.endswith()
        profile_path = os.path.join('db/pubmlst', scheme, scheme + '.txt')
        df_profile = pd.read_csv(profile_path, sep='\t')
        df_profile['profile'] = df_profile[col].apply(
            lambda x: '-'.join(x .astype(str)), axis=1)
        sig = dict(zip(df_profile['profile'], df_profile['ST']))
        genotype[scheme] = {'nloci': count, 'profiles': sig}
        return col, genotype

    @staticmethod
    def best_scheme(result):
        """
        get the best scheme base on the number of found loci
        """
        schemes = []
        length = []
        for item in result.keys():
            schemes.append(item)
            length.append(len(result[item]))
        scheme = schemes[length.index(max(length))]
        return scheme

    # process result
    # {'listeria_2': {'abcZ': '2', 'cat': '11', 'lhkA': '7', 'dat': '3', 'dapE': '3', 'ldh': '1', 'bglA': '1'}}

    @staticmethod
    def get_st(result, scheme):
        """
        get sequence type
        """
        col, genotype = mlst.build_genotype(scheme)
        loci = len(result[scheme])
        sig = ''
        if loci < genotype[scheme]['nloci']:
            st = 'NA'
        else:
            alleles = []
            for i in col:
                alleles.append(result[scheme][i])
            alleles_str = '-'.join(alleles)
            if re.search('[\?~]', alleles_str):
                st = 'NA'
            else:
                st = genotype[scheme]['profiles'][alleles_str]

            df = pd.DataFrame.from_dict(
                dict(zip(col, alleles)), orient='index').T
            df['ST'] = st
            df['Scheme'] = scheme
        return df

    @staticmethod
    def is_fasta(file):
        """
        chcek if the input file is fasta format
        """
        try:
            with open(file, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                # False when `fasta` is empty, i.e. wasn't a FASTA file
                return any(fasta)
        except:
            return False
