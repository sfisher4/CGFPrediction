import sys
import argparse
from cgf_pred import CGFPrediction
import pkg_resources
from glob import glob
import subprocess
from pathlib import Path

def main(args=None):

    DATA_PATH = pkg_resources.resource_filename('cgf_pred', 'fastas/')
    paths = {fasta.stem: pkg_resources.resource_filename('cgf_pred', 'fastas/' + str(fasta.stem) + '.fasta')
             for fasta in Path(DATA_PATH).glob('*.fasta')}

    if args is None:
        args = sys.argv[1: ]
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes",
                        help="a path to a directory that contains fasta files for all genomes of interest",
                        type=str)
    parser.add_argument("out_file",
                       help = "a path to a file that the eCGF results should be printed.",
                       type=str)
    #TODO: add optional argument for validation
    args=parser.parse_args()
    CGFPrediction.main(args.genomes, paths['f_primers'], paths['r_primers'], paths['amp_seq'], args.out_file)

if __name__ == '__main__':
    main()
