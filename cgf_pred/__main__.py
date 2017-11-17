import sys
import argparse
from cgf_pred import CGFPrediction
import pkg_resources
from glob import glob
import subprocess
from pathlib import Path

def arguments():
    parser = argparse.ArgumentParser()

    # parser.add_argument('genomes',
    #                     help="a path to a directory containing fasta files for all query genomes",
    #                     type=str)

    parser.add_argument('out_file',
                        help="a path to a csv file where eCGF results should be printed.",
                        type=str)

    return parser.parse_args()

def main():

    DATA_PATH = pkg_resources.resource_filename('cgf_pred', 'fastas/')
    paths = {fasta.stem: pkg_resources.resource_filename('cgf_pred', 'fastas/' + str(fasta.stem) + '.fasta')
             for fasta in Path(DATA_PATH).glob('*.fasta')}

    DB_PATH = pkg_resources.resource_filename('cgf_pred', 'csvs/')
    db_paths = {txt.stem: pkg_resources.resource_filename('cgf_pred', 'csvs/' + str(txt.stem) + '.txt')
                for txt in Path(DB_PATH).glob('*.txt')}

    args = arguments()
    CGFPrediction.main(args.genomes, args.out_file, paths['f_primers'], paths['r_primers'], paths['amp_seq'],
                       paths['cj0181_f_primer'], paths['cj0181_r_primer'],
                       db_paths['error_rate'], db_paths['cgf_types_fprints'])

if __name__ == '__main__':
    main()
