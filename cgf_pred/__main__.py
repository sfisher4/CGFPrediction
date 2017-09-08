import sys
import argparse
from cgf_pred import CGFPrediction
import pkg_resources
from glob import glob
import subprocess
from pathlib import Path

def main(args=None):

    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # fasta file with primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # fasta file with primer id's and primer sequences
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    all_genomes = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"

    DATA_PATH = pkg_resources.resource_filename('cgf_pred', 'fastas/')
    # DB_FILE_AMP = pkg_resources.resource_filename('cgf_pred', 'fastas/amp_seq.fasta')
    # DB_FILE_FP = pkg_resources.resource_filename('cgf_pred', 'fastas/f_primers.fasta')
    # DB_FILE_RP = pkg_resources.resource_filename('cgf_pred', 'fastas/r_primers.fasta')
    print(type(DATA_PATH))
    # all_fasta_files = glob(DATA_PATH + '*.fasta')
    paths = {fasta.stem: pkg_resources.resource_filename('cgf_pred', fasta)
             for fasta in Path(DATA_PATH).glob('*.fasta')}
    print(paths)
    for path in paths:
        print(path)

    if args is None:
        args = sys.argv[1: ]
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes",
                        help="a path to a directory that contains fasta files for all genomes of interest",
                        type=str)
    args=parser.parse_args()
    CGFPrediction.main(args.genomes, paths['f_primers'], paths['r_primers'], paths['amp_seq'])
    # CGFPrediction.main(args.genomes, all_fasta_files[1], all_fasta_files[2], all_fasta_files[0])

if __name__ == '__main__':
    main()
