import sys
import argparse
from cgf_pred import CGFPrediction
import os


def main(args=None):
    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # fasta file with primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # fasta file with primer id's and primer sequences
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    all_genomes = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"

    if args is None:
        args = sys.argv[1: ]
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes", help="a path to a directory that contains fasta files for all genomes of interest", type=str)
    args=parser.parse_args()
    CGFPrediction.main(args.genomes, forward_primers, reverse_primers, amplicon_sequences)

if __name__ == '__main__':
    main()
