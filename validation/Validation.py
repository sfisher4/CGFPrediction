import gc
import os
from collections import defaultdict
from memory_profiler import profile
from collections import ChainMap
import pandas as pd

from cgf_pred import CGFPrediction

# @profile
def create_binary_table(forward_primers, reverse_primers, amplicon_sequences, db_directory, table_file, lab_binary_results, ehyb_only):
    # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")

    gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307', 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728',
                 'cj0733', 'cj0736', 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294',
                 'cj1324', 'cj1329','cj1334','cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552',
                 'cj1585', 'cj1679', 'cj1721','cj1727c']

    file = open(lab_binary_results, "r")
    lines = file.readlines()
    file_gene_dict = {}
    for line in lines:
        genes_found_list = []
        count = -2
        print(len(line.split()))
        for word in line.split():
            print(word)
            if word == '1':
                genes_found_list.append(gene_list[count])
            count += 1
        # print(genes_found_list)
        file_gene_dict[line.split(None, 1)[0]] = genes_found_list
        # print(file_gene_dict)
    file.close()

    files = (file for file in os.listdir(db_directory) if file.endswith('.fasta'))
    files_paths = []
    for file in files:
        files_paths.append(os.path.abspath(db_directory) + '/' + file)

    for file_path in files_paths:
        file_name = file_path.partition(db_directory + "/")[2]
        print(file_name)
        cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
                                             amplicon_sequences, file_gene_dict)
        cgf_predictions_dict = cgf_predictions[0]
        pred_exception_dict = cgf_predictions[2]
        assert file_name in file_gene_dict, "error: The db directory contains a file that does not have validation results in 'file'"
        genes_expected = file_gene_dict[file_name]
        genes_found_without_exceptions = [key[6:] for key in cgf_predictions_dict.keys()]

        #for testing
        testing = cgf_predictions_dict['11168_cj1134']
        for lo_hsp in testing:
            for hsp in lo_hsp:
                print('!!! BSR', hsp.name, hsp.bsr)

        gene_exceptions_found = [key[6:] for key in pred_exception_dict.keys()]
        genes_found = ChainMap({}, genes_found_without_exceptions, gene_exceptions_found)
        print(genes_found)

        false_positive = set(genes_found) - set(genes_expected)
        false_negative = set(genes_expected) - set(genes_found)
        print('genes found', genes_found)
        print('genes expected', genes_expected)

        if len(false_positive) > 0 or len(false_negative) > 0:
            print(file_name)
        if len(false_positive) > 0:
            print('false positive', false_positive)
        if len(false_negative) > 0:
            print('false negative', false_negative)

        #PANDAS
        df_dict = {}
        for gene in gene_list:
            if gene in false_positive:
                df_dict[gene] = "+1"
            elif gene in false_negative:
                df_dict[gene] = "-1"
            else:
                df_dict[gene] = "0"
        df = pd.DataFrame(df_dict, index=[file_name[:-6]])
        df.to_csv(table_file, sep=' ', mode='a', header=None)

    #F+/F- Calculation
    t_file = open(table_file)
    t_lo_f_pos = {}
    t_lo_f_neg = {}
    t_lines = t_file.readlines()
    for line in t_lines:
        count = 0
        for word in line.split():
            # print(word)
            if word == "+1":
                try:
                    t_lo_f_pos[gene_list[count]] += 1
                except:
                    t_lo_f_pos[gene_list[count]] = 1
                count += 1
            elif word == "-1":
                try:
                    t_lo_f_neg[gene_list[count]] += 1
                except:
                    t_lo_f_neg[gene_list[count]] = 1
                count += 1
            elif word == "0":
                if gene_list[count] not in t_lo_f_neg:
                    t_lo_f_neg[gene_list[count]] = 0
                if gene_list[count] not in t_lo_f_pos:
                    t_lo_f_pos[gene_list[count]] = 0
                count += 1
    t_file.close()

    with open(table_file, "a") as myfile:
        myfile.write('\n' + "\n" + "# of F+: " + str(t_lo_f_pos) + "\n")
        myfile.write('\n' + "\n" + "# of F-: " + str(t_lo_f_neg) + "\n")





if __name__ == "__main__":

    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"
    # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
    lab_binary_results = "/home/sfisher/Sequences/11168_test_files/cgf40_v2.txt"

    #Output text files
    table_file = "/home/sfisher/Sequences/11168_test_files/tables/oct_30_all_all_singlelen_bsr30.txt"
    # short_file = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/30_ehyb_short_expl.txt"
    # gene_file = '/home/sfisher/Sequences/11168_test_files/eCGF_causation/30_ehyb_expl_per_gene.txt'
    # long_file = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/6_ehyb_full_expl.txt"
    # med_file_f_neg = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/6_ehyb_short_f_neg.txt"
    # med_file_f_pos = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/6_ehyb_short_f_pos.txt"
    # per_gene_f_pos = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/15_ecgf_f+_causes.txt"

    ehyb_only = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307',
                 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
                 'cj0736',
                 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
                 'cj1334',
                 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
                 'cj1727c']

    create_binary_table(forward_primers, reverse_primers, amplicon_sequences, db_directory, table_file, lab_binary_results, ehyb_only)
    # create_short_val_results(forward_primers, reverse_primers, amplicon_sequences, db_directory, short_file, lab_binary_results, ehyb_only)
    # create_gene_results(forward_primers, reverse_primers, amplicon_sequences, db_directory, gene_file, lab_binary_results, ehyb_only)
    # create_long_causation(forward_primers, reverse_primers, amplicon_sequences, db_directory, long_file, med_file_f_neg, med_file_f_pos, lab_binary_results, ehyb_only)
    # per_gene_long_expl_f_pos(forward_primers,reverse_primers,amplicon_sequences, db_directory, per_gene_f_pos, lab_binary_results, ehyb_only)