import unittest
import CGFPrediction
import errno
import os
import gc
import pandas as pd
import numpy as np

MIN_BSR = 0.6

class MyTestCase(unittest.TestCase):
    def test_cgf_prediction(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
        db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes/already_debugged"
        # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
        # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
        file = open("/home/sfisher/Sequences/11168_test_files/cgf_246_results_modified.txt", "r")
        table_file = "/home/sfisher/Sequences/11168_test_files/tables/16_epcr_246-gnomes.txt"

        # bsr
        f_bs_primer_dir = "/home/sfisher/Sequences/BSR/f_primers/"
        r_bs_primer_dir = "/home/sfisher/Sequences/BSR/r_primers/"
        amp_bs_dir = "/home/sfisher/Sequences/BSR/amp_seq/"
        max_f_bits_dict = CGFPrediction.max_bs(f_bs_primer_dir)
        max_r_bits_dict = CGFPrediction.max_bs(r_bs_primer_dir)
        max_amp_bits_dict = CGFPrediction.max_bs(amp_bs_dir)

        lines = file.readlines()

        gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                     'cj0307', 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728',
                     'cj0733', 'cj0736', 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294',
                     'cj1324', 'cj1329','cj1334','cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552',
                     'cj1585', 'cj1679', 'cj1721','cj1727c']

        file_gene_dict = {}
        for line in lines:
            genes_found_list = []
            count = -1
            for word in line.split():
                if word == '1':
                    genes_found_list.append(gene_list[count])
                count += 1
            file_gene_dict[line.split(None, 1)[0]] = genes_found_list

        file.close()

        files = (file for file in os.listdir(db_directory) if file.endswith('.fasta'))
        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        for file_path in files_paths:
            file_name = file_path.partition(db_directory + "/")[2]
            print(file_name)
            cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
                                                 amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)
            cgf_predictions_dict = cgf_predictions[0]
            # print('dict', cgf_predictions_dict)

            genes_expected = file_gene_dict[file_name]
            genes_found = [key for key in cgf_predictions_dict]

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
                print(word)
                if word == "+1":
                    try:
                        t_lo_f_pos[gene_list[count]] += 1
                    except:
                        t_lo_f_pos[gene_list[count]] = 1
                    count += 1
                if word == "-1":
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
            myfile.write("\n" + "# of F+: " + str(t_lo_f_pos) + "\n")
            myfile.write("\n" + "# of F-: " + str(t_lo_f_neg) + "\n")





            # def test_cgf_prediction_debug(self):
    #     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    #     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    #     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    #     db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"
    #     # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    #     # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
    #     file = open("/home/sfisher/Sequences/11168_test_files/cgf_246_results_modified.txt", "r")
    #
    #     # bsr
    #     f_bs_primer_dir = "/home/sfisher/Sequences/BSR/f_primers/"
    #     r_bs_primer_dir = "/home/sfisher/Sequences/BSR/r_primers/"
    #     amp_bs_dir = "/home/sfisher/Sequences/BSR/amp_seq/"
    #     max_f_bits_dict = CGFPrediction.max_bs(f_bs_primer_dir)
    #     max_r_bits_dict = CGFPrediction.max_bs(r_bs_primer_dir)
    #     max_amp_bits_dict = CGFPrediction.max_bs(amp_bs_dir)
    #
    #     lines = file.readlines()
    #
    #     gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
    #                  'cj0307',
    #                  'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
    #                  'cj0736',
    #                  'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
    #                  'cj1334',
    #                  'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
    #                  'cj1727c']
    #
    #     file_gene_dict = {}
    #     for line in lines:
    #         genes_found_list = []
    #         count = -1
    #         for word in line.split():
    #             if word == '1':
    #                 genes_found_list.append(gene_list[count])
    #             count += 1
    #         file_gene_dict[line.split(None, 1)[0]] = genes_found_list
    #
    #     file.close()
    #
    #     files = (file for file in os.listdir(db_directory) if file.endswith('.fasta'))
    #     files_paths = []
    #     for file in files:
    #         files_paths.append(os.path.abspath(db_directory) + '/' + file)
    #
    #     for file_path in files_paths:
    #         file_name = file_path.partition(db_directory + "/")[2]
    #         print(file_name)
    #         cgf_predictions = CGFPrediction.cgf_prediction(forward_primers, reverse_primers, file_path,
    #                                                        amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict, True)
    #         cgf_predictions_dict = cgf_predictions[0]
    #         # print('dict', cgf_predictions_dict)
    #         all_hsp = cgf_predictions[1]
    #         print(type(all_hsp)) #default dictionary
    #
    #         genes_expected = file_gene_dict[file_name]
    #         genes_found = [key for key in cgf_predictions_dict]
    #
    #         false_positive = set(genes_found) - set(genes_expected)
    #         false_negative = set(genes_expected) - set(genes_found)
    #         print('genes found', genes_found)
    #         print('genes expected', genes_expected)
    #
    #         if len(false_positive) > 0 or len(false_negative) > 0:
    #             print(file_name)
    #         if len(false_positive) > 0:
    #             print('false positive', false_positive)
    #         if len(false_negative) > 0:
    #             print('false negative', false_negative)
    #
    #         binary_tree_attrib = ['both_primers_found', 'contig', 'ehybrid', 'ehybrid', 'ehybrid', 'location', None,
    #                        'epcr', None, 'location']
    #         for i in range(10, 16):
    #             binary_tree_attrib.insert(i, None)
    #         binary_tree_attrib.insert(16, 'valid')
    #         for i in range(17, 23):
    #             binary_tree_attrib.insert(i, None)
    #         for i in range(23, 33):
    #             binary_tree_attrib.insert(i, None)
    #         binary_tree_attrib.insert(33, 'snp')
    #         for i in range(34, 68):
    #             binary_tree_attrib.insert(i, None)
    #         binary_tree_attrib.insert(68, 'pcr_distance')
    #         for i in range(69, 139):
    #             binary_tree_attrib.insert(i, None)
    #
    #         for gene_name in false_negative:
    #             for hsp in all_hsp[gene_name]:
    #                 if hsp.name in gene_name:
    #                     result = CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0)
    #                     if hsp.bsr < MIN_BSR:
    #                         continue
    #                         # print(hsp.name, "on contig", hsp.contig_name, "failed because hsp bsr is", hsp.bsr)
    #                     print('bsr', hsp.bsr)
    #                     if result == 2:
    #                         print(hsp.name, "CASE 2")
    #                         print('both primers found:', hsp.both_primers_found)
    #                         print('ehybrid is:', hsp.ehybrid)
    #                         # print('query', hsp.query)
    #                         print('sbjct', hsp.sbjct)
    #                         # print('amp query', hsp.amp_query)
    #                         # print('amp sbjct', hsp.amp_sbjct)
    #                     elif result == 3:
    #                         print(hsp.name, "CASE 3")
    #                         print('both primers found:', hsp.both_primers_found, 'with contig', hsp.contig)
    #                         print('ehybrid is:', hsp.ehybrid)
    #                         print('amp query', hsp.amp_query)
    #                         print('amp sbjct', hsp.amp_sbjct)
    #                     elif result == 4:
    #                         print(hsp.name, "on contig", hsp.contig_name, "CASE 4")
    #                         print('both primers found:', hsp.both_primers_found, 'with contig', hsp.contig)
    #                         print('ehybrid is:', hsp.ehybrid)
    #                         # print('query', hsp.query)
    #                         print('strand', hsp.strand)
    #                         print('sbjct', hsp.sbjct)
    #                         print('partner & partner sbjct', hsp.partner.name, hsp.partner.contig_name, hsp.partner.sbjct)
    #                     elif result == 6 or result == 8 or result == 10:
    #                         print(hsp.name, "on contig", hsp.contig_name, "CASE 6 or 8 or 10")
    #                         print('contig', hsp.contig)
    #                         print('both primers found:', hsp.both_primers_found)
    #                         print('ehybrid is:', hsp.ehybrid, 'with len', hsp.length)
    #                         # print('query', hsp.query)
    #                         print('sbjct', hsp.sbjct)
    #                         print('amp query', hsp.amp_query)
    #                         print('amp sbjct', hsp.amp_sbjct)
    #                     elif result == 11 or result == 15 or result == 137:
    #                         print('FALSE')
    #                     elif result == 12:
    #                         print(hsp.name, "CASE 12")
    #                         print('both primers found:', hsp.both_primers_found, 'with ehybrid', hsp.ehybrid, 'and location', hsp.location)
    #                         if hsp.strand:
    #                             print('leading strand', 'located', hsp.end_dist, 'from end')
    #                         elif not hsp.strand:
    #                             print('lagging strand', 'located,', hsp.end_dist, 'from end')
    #                         # print('query', hsp.query)
    #                         print('sbjct', hsp.sbjct)
    #                         # print('amp query', hsp.amp_query)
    #                         # print('amp sbjct', hsp.amp_sbjct)
    #                         print('Most likely too restrictive to find both primers...')
    #                     elif result == 20:
    #                         print(hsp.name, "on contig", hsp.contig_name, "CASE 20")
    #                         print('both primers found:', hsp.both_primers_found, 'with contig:', hsp.contig, 'and location', hsp.location)
    #                         print('strand', hsp.strand)
    #                         print('sbjct', hsp.sbjct)
    #                         print('partner', hsp.partner.contig_name)
    #                         if hsp.strand:
    #                             print('leading strand', 'located', hsp.end_dist, 'from end')
    #                             # print('primer query:', hsp.query)
    #                         elif not hsp.strand:
    #                             print('lagging strand', 'located,', hsp.end_dist, 'from end')
    #                         print('amp query', hsp.amp_query)
    #                         print('amp sbjct', hsp.amp_sbjct)
    #                         # print('query    ', hsp.query)
    #                     elif result == 34:
    #                         print(hsp.name, "CASE 34")
    #                         print('both primers found:', hsp.both_primers_found, 'with contig:', hsp.contig)
    #                         print('invalid direction to other hsp')
    #                     elif result == 67:
    #                         print(hsp.name, "CASE 67")
    #                         print('false b/c snp:', hsp.snp)
    #                         print('strand', hsp.strand)
    #                         print('sbjct', hsp.sbjct)
    #                         print('contig', hsp.contig_name)
    #                         # print('query', hsp.query)
    #                     elif result == 138:
    #                         print(hsp.name, "CASE 138")
    #                         print('false b/c dist btwn primers')
    #                     else:
    #                         assert False
    #                     print('\n')
    #             del all_hsp
    #             gc.collect()
    #
    #
    #         for gene_name in false_positive:
    #             for lo_hsp in cgf_predictions_dict[gene_name]:
    #                 for hsp in lo_hsp:
    #                     print(hsp.name)
    #                     print(hsp.bsr)
    #                     print('strand', hsp.strand)
    #                     print('sbjct    ', hsp.sbjct)
    #                     print('\n')








if __name__ == '__main__':
    unittest.main()
















