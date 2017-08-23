import unittest
import CGFPrediction
import errno
import os
import gc
import pandas as pd
import numpy as np
from collections import defaultdict
from Bio.Seq import Seq


MIN_BSR = 0.7
MIN_PARTNER_BSR = 0.7

class MyTestCase(unittest.TestCase):
    def test_cgf_prediction(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
        db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes/already_debugged"
        # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
        # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
        file = open("/home/sfisher/Sequences/11168_test_files/cgf40_v2.txt", "r")
        table_file = "/home/sfisher/Sequences/11168_test_files/tables/23_ehyb_all.txt"

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
            count = -2
            for word in line.split():
                # print(word)
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
                                                 amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)
            cgf_predictions_dict = cgf_predictions[0]
            # print('dict', cgf_predictions_dict)
            assert file_name in file_gene_dict, "error: The db directory contains a file that does not have validation results in 'file'"
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
                # print(word)
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
            myfile.write('\n' + "\n" + "# of F+: " + str(t_lo_f_pos) + "\n")
            myfile.write('\n' + "\n" + "# of F-: " + str(t_lo_f_neg) + "\n")





    # def test_cgf_prediction_debug(self):
    #     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    #     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    #     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    #     db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes/debug_f_neg"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
    #     # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"
    #     # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    #     # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
    #     pfile = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/22_246_v2.txt"
    #     short_file_f_neg = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/short_expl_f_neg.txt"
    #     short_file_f_pos = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/short_expl_f_pos.txt"
    #     file = open("/home/sfisher/Sequences/11168_test_files/cgf40_v2.txt", "r")
    #
    #     f_primer_dict = CGFPrediction.create_primer_dict(forward_primers)
    #     r_primer_dict = CGFPrediction.create_primer_dict(reverse_primers)
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
    #         count = -2
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
    #         print_file = open(pfile, "a")
    #         file_name = file_path.partition(db_directory + "/")[2]
    #         cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
    #                                                        amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict, True)
    #         cgf_predictions_dict = cgf_predictions[0]
    #         all_hsp = cgf_predictions[1]
    #
    #         genes_expected = file_gene_dict[file_name]
    #         genes_found = [key for key in cgf_predictions_dict]
    #
    #         false_positive = set(genes_found) - set(genes_expected)
    #         false_negative = set(genes_expected) - set(genes_found)
    #
    #         if len(false_positive) > 0 or len(false_negative) > 0:
    #             # print_file.write('\n' + file_name)
    #             print_file.write('\n' + '\n' + file_name+ '\n')
    #             print_file.write('\n' + 'genes found'+ str(genes_found))
    #             print_file.write('\n' + 'genes expected'+ str(genes_expected))
    #         if len(false_positive) > 0:
    #             print_file.write('\n' + 'false positive'+ str(false_positive))
    #         if len(false_negative) > 0:
    #             print_file.write('\n' + 'false negative'+ str(false_negative) + '\n')
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
    #         ehyb_only = ['cj0421c', 'cj0246c', 'cj0755', 'cj1294', 'cj1324', 'cj1334', 'cj1427c', 'cj1721',
    #                          'cj1329', 'cj1439', 'cj1552', 'cj1551']
    #         for gene_name in false_negative:
    #             cases = defaultdict(int)
    #             hi_bsr = 0
    #             if gene_name in ehyb_only:
    #                 print_file.write('\n' + gene_name + ' not found using ehyb only')
    #                 #TODO!!!
    #             else:
    #                 for hsp in all_hsp[gene_name]:
    #                     if hsp.name in gene_name:
    #                         result = CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0)
    #                         if hsp.bsr < MIN_BSR:
    #                             continue
    #                             # print_file.write('\n' + hsp.name, "on contig", hsp.contig_name, "failed because hsp bsr is", hsp.bsr)
    #                         if hsp.partner != None:
    #                             if hsp.partner.bsr < MIN_PARTNER_BSR:
    #                                 continue
    #                         if hsp.bsr > hi_bsr:
    #                             hi_bsr = hsp.bsr
    #                         print_file.write('\n' + gene_name + ' found with bsr: '+ str(hsp.bsr))
    #                         cases[result] += 1
    #                         if result == 6:
    #                             print_file.write('\n' + 'CASE 6:' + gene_name + ' not found using eCGF with both primers so one primer'
    #                                                     'searched for with EHYBRID ' + '(' + str(hsp.ehybrid) +') not found of long enough length (len=' + str(hsp.length) + ")")
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                         elif result == 10:
    #                             print_file.write('\n' + 'CASE 10: ' + gene_name + ' Found primers on DIFFERENT contigs and '
    #                                                     'EHYBRID: ' + '(' + str(hsp.ehybrid) +')' + ' found but not ' \
    #                                                                                                     'of long enough length')
    #                             print_file.write('\n' + 'Contigs:' + hsp.contig_name + hsp.partner.contig_name)
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 12:
    #                             print_file.write('\n' + 'CASE 12:' + gene_name + ' not found using eCGF with both primers so one primer'
    #                                                     'searched for with EHYBRID found ' + '(' + str(hsp.ehybrid) +')' + 'but not located at the end of a contig. (hsp location is ')
    #                             if hsp.strand:
    #                                 print_file.write(str(hsp.end_dist) + 'from end (on leading strand)')
    #                             elif not hsp.strand:
    #                                 print_file.write(str(hsp.end_dist) + 'from end (on lagging strand)')
    #                             print_file.write('\n' + 'Most likely too restrictive to find both primers...')
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                         elif result == 2:
    #                             print_file.write('\n' + "CASE 2: " + gene_name + ' Not found using eCGF with both primers so one primer'
    #                                                     'searched for with EHYBRID'+ '(' + str(hsp.ehybrid) +')' + 'not found at all!')
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                         elif result == 3:
    #                             print_file.write('\n' + 'CASE 3: '+ gene_name + ' Both primers found on the same contig but EHYBRID'+ '(' + str(hsp.ehybrid) +')' + 'not found at all!')
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 4:
    #                             print_file.write('\n' + "CASE 4: "+ gene_name + ' Found primers on DIFFERENT contigs and'
    #                                                     'EHYBRID: ' + '(' + str(hsp.ehybrid) +')' + 'not found at all!')
    #                             print_file.write('\n' + 'Contigs:' + hsp.contig_name + hsp.partner.contig_name)
    #                             print_file.write('\n' + 'strand'+ str(hsp.strand))
    #                             print_file.write('\n' + 'sbjct'+ hsp.sbjct)
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 8:
    #                             print_file.write('\n' + 'CASE 8: '+ gene_name + ' Both primers found on the same contig '
    #                                                                             'with EHYBRID'+ '(' + str(hsp.ehybrid) +')' +
    #                                              'not found b/c length too short: len=' + str(hsp.length))
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 20:
    #                             print_file.write('\n' + 'CASE 20: '+ gene_name + ' Found primers on DIFFERENT contigs and'
    #                                                     'EHYBRID: ' + '(' + str(hsp.ehybrid) +')' + ' found but at the end of a contig. (hsp location is ')
    #                             if hsp.strand:
    #                                 print_file.write(str(hsp.end_dist) + ' from end (leading strand)')
    #                             elif not hsp.strand:
    #                                 print_file.write(str(hsp.end_dist) + ' from end (lagging strand)')
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 11 or result == 15 or result == 137:
    #                             print_file.write('\n' + 'FALSE')
    #                         elif result == 34:
    #                             print_file.write('\n' + "CASE 34: " + gene_name + " Both primers on the same contig with EHYBRID"
    #                                              + '(' + str(hsp.ehybrid) +')' + "found but epcr not found" + '(' + str(hsp.epcr) +')' +
    #                                              'because strands are not facing each other (valid =' + hsp.valid + ')')
    #                             print_file.write('\n' + 'both primers found:'+ str(hsp.both_primers_found) + 'with contig:'+ str(hsp.contig))
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 138:
    #                             print_file.write('\n' + "CASE 138: " + gene_name + " Both primers on the same contig with EHYBRID"
    #                                              + '(' + str(hsp.ehybrid) +')' + "found but epcr not found" + '(' + str(hsp.epcr) +')' +
    #                                              'b/c primers are not the correct distance from each other')
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #                         elif result == 67:
    #                             print_file.write('\n' + "CASE 67: " + gene_name + " Both primers on the same contig with EHYBRID"
    #                                              + '(' + str(hsp.ehybrid) +')' + "found but epcr not found" + '(' + str(hsp.epcr) +')' +
    #                                              'because of SNP on "3'" end")
    #                             print_file.write('\n' + 'strand '+ str(hsp.strand))
    #                             print_file.write('\n' + 'query '+ hsp.query)
    #                             print_file.write('\n' + 'sbjct '+ hsp.sbjct)
    #                             print_file.write('\n' + 'f_primer seq: ' + str(f_primer_dict[gene_name]))
    #                             print_file.write('\n' + 'r_primer seq: ' + str(r_primer_dict[gene_name]))
    #                             print_file.write('\n' + 'r_comp_f_pri: ' + str(f_primer_dict[gene_name].reverse_complement()))
    #                             print_file.write('\n' + 'r_comp_r_pri: ' + str(r_primer_dict[gene_name].reverse_complement()))
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
    #
    #                         else:
    #                             print_file.write('\n' + str(result))
    #                             assert False
    #                         print_file.write('\n' + '\n')
    #             print_file.write('\n' + '\n' + 'Cases and their occurences for: '+ gene_name + ' ' + str(cases))
    #             print_file.write('\n' + 'Largest bsr: ' + str(hi_bsr) + '\n \n')
    #
    #             short_file = open(short_file_f_neg, "a")
    #             short_file.write('\n' + '\n' + file_name)
    #             short_file.write('\n' + '\n' + gene_name)
    #             short_file.write('\n' + 'F+: ' + str(false_positive))
    #             short_file.write('\n' + 'F-: ' + str(false_negative))
    #             if gene_name in ehyb_only:
    #                 short_file.write('\n' + 'Not found using ehyb only.')
    #                 continue
    #
    #             short_file.write('\n' + 'Num/Case: ' + str(cases))
    #             short_file.write('\n' + 'Highest BSR ' + str(hi_bsr))
    #             #TODO: determine where the highest BSR came from.
    #             if 67 in cases:
    #                 short_file.write('\n' + '67:   SAME CONTIG --> SNP on 3 prime end')
    #             if 138 in cases:
    #                 short_file.write('\n' + '138:  SAME CONTIG --> Dist between primers not correct')
    #             if 34 in cases:
    #                 short_file.write('\n' + '34:   SAME CONTIG --> Primers are not facing each other')
    #             if 3 in cases:
    #                 short_file.write('\n' + '3:    SAME CONTIG --> eHybrid is not found in the same position as the primers found')
    #             if 8 in cases:
    #                 short_file.write('\n' + '8:    SAME CONTIG --> eHybrid was found but not of long enough length')
    #             if 12 in cases:
    #                 short_file.write('\n' + '12:   ONE PRIMER --> eHybrid found but not at end of contig... most likely too restrictive to find both primers.')
    #             if 20 in cases:
    #                 short_file.write('\n' + '20:   DIFF CONTIGS --> eHybrid found but not at end of contig... most likely too restrictive to find both primers')
    #             if 2 in cases:
    #                 short_file.write('\n' + '2:    ONE PRIMER --> eHybrid is not found in the same position as the primers found')
    #             if 4 in cases:
    #                 short_file.write('\n' + '4:    DIFF CONTIG --> eHybrid is not found in the same position as the primers found')
    #             if 6 in cases:
    #                 short_file.write('\n' + '6:    ONE PRIMER --> eHybrid was found but not of long enough length')
    #             if 10 in cases:
    #                 short_file.write('\n' + '10:   DIFF CONTIG --> eHybrid was found but not of long enough length')
    #
    #             short_file.close()
    #
    #
    #         del all_hsp
    #         gc.collect()
    #
    #         #TODO: debug... F+ not showing up!
    #         ehyb_only = ['cj0421c', 'cj0246c', 'cj0755', 'cj1294', 'cj1324', 'cj1334', 'cj1427c', 'cj1721',
    #                          'cj1329', 'cj1439', 'cj1552', 'cj1551']
    #
    #         short_file = open(short_file_f_pos, "a")
    #         short_file.write(file_name)
    #         for gene_name in false_positive:
    #             if gene_name in ehyb_only:
    #                 print_file.write('\n' + gene_name + ' found using ehyb only')
    #             else:
    #                 for lo_hsp in cgf_predictions_dict[gene_name]:
    #                 # print_file.write('\n' + type(lo_hsp))
    #                     if lo_hsp != None:
    #                         for hsp in lo_hsp:
    #                             print_file.write('\n' + hsp.name)
    #                             print_file.write('\n' + str(hsp.bsr))
    #                             print_file.write('\n' + 'strand'+ str(hsp.strand))
    #                             print_file.write('\n' + 'query         '+ hsp.query)
    #                             print_file.write('\n' + 'sbjct         '+ hsp.sbjct)
    #                             print_file.write('\n' + 'f_primer seq: ' + str(f_primer_dict[gene_name]))
    #                             print_file.write('\n' + 'r_primer seq: ' + str(r_primer_dict[gene_name]))
    #                             print_file.write('\n' + 'r_comp_f_pri: ' + str(f_primer_dict[gene_name].reverse_complement()))
    #                             print_file.write('\n' + 'r_comp_r_pri: ' + str(r_primer_dict[gene_name].reverse_complement()))
    #                             print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             print_file.write('\n' + '\n')
    #                             short_file.write('\n' + hsp.name)
    #                             short_file.write('\n' + str(hsp.bsr))
    #                             short_file.write('\n' + 'strand'+ str(hsp.strand))
    #                             short_file.write('\n' + 'query         '+ hsp.query)
    #                             short_file.write('\n' + 'sbjct         '+ hsp.sbjct)
    #                             short_file.write('\n' + 'f_primer seq: ' + str(f_primer_dict[gene_name]))
    #                             short_file.write('\n' + 'r_primer seq: ' + str(r_primer_dict[gene_name]))
    #                             short_file.write('\n' + 'r_comp_f_pri: ' + str(f_primer_dict[gene_name].reverse_complement()))
    #                             short_file.write('\n' + 'r_comp_r_pri: ' + str(r_primer_dict[gene_name].reverse_complement()))
    #                             short_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
    #                             short_file.write('\n' + '\n')
    #                     else:
    #                         print_file.write('\n' + 'no hsp found')
    #         print_file.close()
    #         short_file.close()


if __name__ == '__main__':
    unittest.main()
















