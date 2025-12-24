#!/usr/bin/python3
# -*- coding:utf-8 -*-

import os
import re
import sys
import subprocess
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np

from cvmcore.cvmcore import cfunc


class MLST:  # Renamed to follow Python naming conventions
    def __init__(self, inputfile: str, database: str, output: str, threads: int,
                 minid: float = 90, mincov: float = 60):
        """Initialize MLST analysis.

        Args:
            inputfile: Path to input FASTA file
            database: Path to MLST database
            output: Output directory path
            threads: Number of threads to use
            minid: Minimum identity percentage (default: 90)
            mincov: Minimum coverage percentage (default: 60)
        """
        if not os.path.exists(inputfile):
            raise FileNotFoundError(f"Input file not found: {inputfile}")
        if not cfunc.is_fasta(inputfile):
            raise ValueError(f"Input file is not in FASTA format: {inputfile}")

        self.inputfile = os.path.abspath(inputfile)
        self.database = database
        self.minid = float(minid)
        self.mincov = float(mincov)
        self.temp_output = os.path.join(os.path.abspath(output), 'temp.txt')
        self.threads = threads
        self.blast_type = 'blastn'  # Add blast type attribute

    def process_blast_results(self, df: pd.DataFrame) -> Dict:
        """Process BLAST results and identify alleles.

        Args:
            df: DataFrame containing BLAST results

        Returns:
            Dictionary containing processed MLST results
        """
        result = {}
        length_filter = {}

        for _, row in df.iterrows():
            try:
                sch, gene, num = re.match(
                    r'^(\w+)\.(\w+)[_-](\d+)', row['sseqid']).group(1, 2, 3)
            except AttributeError:
                continue

            hlen = row['slen']
            alen = row['length']
            nident = row['nident']

            if nident * 100 / hlen < self.mincov:
                continue

            if sch not in result:
                result[sch] = {}
                length_filter[sch] = {}

            exact_match = hlen == alen and nident == hlen
            if exact_match:
                self._handle_exact_match(
                    sch, gene, num, hlen, result, length_filter)
            elif alen == hlen and nident != hlen:
                self._handle_new_allele(sch, gene, num, result)
            elif alen != hlen and nident == hlen:
                self._handle_partial_match(sch, gene, num, result)

        return result

    def _handle_exact_match(self, sch: str, gene: str, num: str, hlen: int,
                            result: Dict, length_filter: Dict) -> None:
        """Handle exact matches in BLAST results."""
        if gene not in length_filter[sch]:
            length_filter[sch][gene] = hlen

        if gene in result[sch]:
            if not re.search(r'[~\?]', result[sch][gene]):
                if hlen > length_filter[sch][gene]:
                    result[sch][gene] = num
                    length_filter[sch][gene] = hlen
                elif hlen == length_filter[sch][gene]:
                    result[sch][gene] = f"{result[sch][gene]}, {num}"
            else:
                result[sch][gene] = num
                length_filter[sch][gene] = hlen
        else:
            result[sch][gene] = num
            length_filter[sch][gene] = hlen

    def _handle_new_allele(self, sch: str, gene: str, num: str, result: Dict) -> None:
        """Handle new alleles in BLAST results."""
        if gene not in result[sch]:
            result[sch][gene] = f'~{num}'

    def _handle_partial_match(self, sch: str, gene: str, num: str, result: Dict) -> None:
        """Handle partial matches in BLAST results."""
        if gene not in result[sch]:
            result[sch][gene] = f'{num}?'

    def biopython_blast(self):
        # cline = NcbiblastnCommandline(query=self.inputfile, db=self.database, dust='no', ungapped=True,
        #                               evalue=1E-20, out=self.temp_output,  # delete culling_limit parameters
        #                               outfmt="6 sseqid slen length nident",
        #                               perc_identity=self.minid, max_target_seqs=10000,
        #                               num_threads=self.threads)
        # stdout, stderr = cline()

        cline = [self.blast_type, '-query', self.inputfile, '-db', self.database, '-dust', 'no', '-ungapped',
                 '-evalue', '1E-20', '-out', self.temp_output,
                 '-outfmt', '6 sseqid slen length nident', '-perc_identity', str(
                     self.minid), '-max_target_seqs', '10000',
                 '-num_threads', str(self.threads)]

        # print(cline)

        # Run the command using subprocess
        cline_result = subprocess.run(
            cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Capture the output and error
        stdout = cline_result.stdout
        stderr = cline_result.stderr

        # Print or handle the output and error as needed
        # print(stdout)
        if stderr:
            print(f"Error: {stderr}")

        df = pd.read_csv(self.temp_output, sep='\t', names=[
            'sseqid', 'slen', 'length', 'nident'])
        # df.to_csv('test.csv')
        # print(df)

        result = self.process_blast_results(df)
        # remove temp blastn output file
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
        scheme_path = os.path.abspath(os.path.dirname(__file__))
        count = 0
        col = []
        genotype = {}
        db_path = os.path.join(os.path.join(scheme_path, 'db/pubmlst'), scheme)
        # print(db_path)
        for file in os.listdir(db_path):
            # print(file)
            if file.endswith('.tfa'):
                base = os.path.splitext(file)[0]
                col.append(base)
                # print(base)
                count += 1
                # if file.endswith()
        profile_path = os.path.join(db_path, scheme + '.txt')
        df_profile = pd.read_csv(profile_path, sep='\t')
        df_profile['profile'] = df_profile[col].apply(
            lambda x: '-'.join(x.astype(str)), axis=1)
        sig = dict(zip(df_profile['profile'], df_profile['ST']))
        genotype[scheme] = {'nloci': count, 'profiles': sig}
        return col, genotype

    @staticmethod
    def get_best_scheme(result: Dict) -> List[str]:
        """Get the best scheme based on number and quality of matches.

        Args:
            result: Dictionary containing MLST results

        Returns:
            List of best matching scheme names
        """
        schemes = []
        scores = []

        for scheme, gene_locus_dict in result.items():
            score = len(gene_locus_dict)
            for allele_num in gene_locus_dict.values():
                if '~' in allele_num:
                    score -= 0.5
                elif '?' in allele_num:
                    score -= 1

            schemes.append(scheme)
            scores.append(score)

        max_score = max(scores)
        return [s for s, score in zip(schemes, scores) if score == max_score]

    # process result
    # {'listeria_2': {'abcZ': '2', 'cat': '11', 'lhkA': '7', 'dat': '3', 'dapE': '3', 'ldh': '1', 'bglA': '1'}}

    @staticmethod
    def get_st(result):
        """
        get sequence type
        """

        # # Get best match scheme
        # schemes = []
        # scores = []
        # for item in result.keys():
        #     schemes.append(item)
        #     gene_locus_dict = result[item]
        #     # solve could not found best scheme bug when scheme (have novel or approximate loci)
        #     # have same loci compared to best scheme
        #     nloci = len(gene_locus_dict)
        #     score = nloci
        #     for gene in gene_locus_dict.keys():
        #         allele_num = gene_locus_dict[gene]
        #         if re.search('~',allele_num):
        #             score -= 0.5
        #         elif re.search(r'\?', allele_num):
        #             score -= 1
        #         else:
        #             score = score
        #     # print(f'score {score}')
        #     scores.append(score)
        #     # length.append(len(result[item]))
        # scheme = schemes[scores.index(max(scores))]

        # Get best schemes
        schemes = MLST.get_best_scheme(result)

        # Get Sequence Type
        df_ST = pd.DataFrame()  # init a empty dataframe
        for scheme in schemes:
            col, genotype = MLST.build_genotype(scheme)
            # print(genotype) # genotype为dict {sig:st}
            loci = len(result[scheme])
            print(result[scheme])
            sig = ''
            if loci < genotype[scheme]['nloci']:
                st = 'NA'
                df_tmp = pd.DataFrame.from_dict(
                    result[scheme], orient='index').T
                df_tmp['Note'] = f'Only found {loci} loci in genome, could not determine ST'

            else:
                alleles = []
                for i in col:
                    alleles.append(result[scheme][i])
                alleles_str = '-'.join(alleles)
                if re.search('[\?~]', alleles_str):
                    st = '-'
                else:
                    if alleles_str in genotype[scheme]['profiles']:
                        st = genotype[scheme]['profiles'][alleles_str]
                    else:
                        st = "NewST"
                df_tmp = pd.DataFrame.from_dict(
                    dict(zip(col, alleles)), orient='index').T
            df_tmp['ST'] = st
            df_tmp['Scheme'] = scheme
            df_ST = pd.concat([df_ST, df_tmp])
        return df_ST
