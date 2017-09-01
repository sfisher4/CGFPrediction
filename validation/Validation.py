import gc
import os
from collections import defaultdict

import pandas as pd

from cgf_pred import CGFPrediction

MIN_PARTNER_BSR = 0.7
MIN_BSR = 0.7

def create_binary_table(forward_primers, reverse_primers, amplicon_sequences, db_directory, table_file, lab_binary_results, ehyb_only):
    # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
    file = open(lab_binary_results, "r")

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

def create_short_val_results(forward_primers, reverse_primers, amplicon_sequences, db_directory, result_file, lab_binary_results, ehyb_only):

    file = open(lab_binary_results, "r")

    f_primer_dict = CGFPrediction.create_primer_dict(forward_primers)
    r_primer_dict = CGFPrediction.create_primer_dict(reverse_primers)

    # bsr
    f_bs_primer_dir = "/home/sfisher/Sequences/BSR/f_primers/"
    r_bs_primer_dir = "/home/sfisher/Sequences/BSR/r_primers/"
    amp_bs_dir = "/home/sfisher/Sequences/BSR/amp_seq/"
    max_f_bits_dict = CGFPrediction.max_bs(f_bs_primer_dir)
    max_r_bits_dict = CGFPrediction.max_bs(r_bs_primer_dir)
    max_amp_bits_dict = CGFPrediction.max_bs(amp_bs_dir)

    lines = file.readlines()

    gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307',
                 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
                 'cj0736',
                 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
                 'cj1334',
                 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
                 'cj1727c']

    file_gene_dict = {}
    for line in lines:
        genes_found_list = []
        count = -2
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
        cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
                                             amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict,
                                             True)
        cgf_predictions_dict = cgf_predictions[0]
        all_hsp = cgf_predictions[1]

        genes_expected = file_gene_dict[file_name]
        genes_found = [key for key in cgf_predictions_dict]

        false_positive = set(genes_found) - set(genes_expected)
        false_negative = set(genes_expected) - set(genes_found)

        # if len(false_positive) > 0 or len(false_negative) > 0:
            # print_file.write('\n' + file_name)
            # print_file.write('\n' + '\n' + file_name)

        binary_tree_attrib = ['both_primers_found', 'contig', 'ehybrid', 'ehybrid', 'ehybrid', 'location', None,
                              'epcr', None, 'location']
        for i in range(10, 16):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(16, 'valid')
        for i in range(17, 23):
            binary_tree_attrib.insert(i, None)
        for i in range(23, 33):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(33, 'snp')
        for i in range(34, 68):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(68, 'pcr_distance')
        for i in range(69, 139):
            binary_tree_attrib.insert(i, None)

        short_file = open(result_file, "a")
        if len(false_negative) > 0 or len(false_positive) > 0:
            short_file.write('\n' + '\n' + file_name)
            if len(false_negative) > 0:
                short_file.write('\n' + 'false negative' + str(false_negative))
            if len(false_positive) > 0:
                short_file.write('\n' + 'false positive' + str(false_positive))

        for gene_name in false_negative:
            cases = defaultdict(int)
            hi_bsr = 0
            # if gene_name in ehyb_only:
            #     print_file.write('\n' + ' Not found using ehyb only')
            # else:

            #construct cases!!!
            for hsp in all_hsp[gene_name]:
                if hsp.bsr < MIN_BSR:
                    continue
                if hsp.name in gene_name:
                    result = CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0)
                        # print_file.write('\n' + hsp.name, "on contig", hsp.contig_name, "failed because hsp bsr is", hsp.bsr)
                    try:
                        if hsp.partner != None:
                            if hsp.partner.bsr < MIN_PARTNER_BSR:
                                continue
                    except:
                        continue
                    if hsp.bsr > hi_bsr:
                        hi_bsr = hsp.bsr
                    # print_file.write('\n' + gene_name + ' found with bsr: ' + str(hsp.bsr))
                    cases[result] += 1
                    if 'cj0570' in gene_name:
                        print('result', result)
                        # print_file.write('\n' + '\n')
            # print_file.write('\n' + '\n' + 'Cases and their occurences for: ' + gene_name + ' ' + str(cases))
            # print_file.write('\n' + 'Largest bsr: ' + str(hi_bsr) + '\n \n')

            # short_file.write('\n' + 'F+: ' + str(false_positive))
            # short_file.write('\n' + 'F-: ' + str(false_negative))
            # short_file.write('\n' + 'Num/Case: ' + str(cases))
            # short_file.write('\n' + 'Highest BSR ' + str(hi_bsr))
            short_file.write('\n' + '\n' + 'F-: '+ gene_name)
            if 138 in cases:
                short_file.write('\n' + '138:  SAME CONTIG --> Dist between primers not correct')
            elif 67 in cases:
                short_file.write('\n' + '67:   SAME CONTIG --> SNP on 3 prime end')
            elif 34 in cases:
                short_file.write('\n' + '34:   SAME CONTIG --> Primers are not facing each other')
            elif 3 in cases:
                short_file.write(
                    '\n' + '3:    SAME CONTIG --> eHybrid is not found in the same position as the primers found')
            elif 8 in cases:
                short_file.write('\n' + '8:    SAME CONTIG --> eHybrid was found but not of long enough length')
            elif 12 in cases:
                short_file.write(
                    '\n' + '12:   ONE PRIMER --> eHybrid found but not at end of contig... most likely too restrictive to find both primers.')
            elif 20 in cases:
                short_file.write(
                    '\n' + '20:   DIFF CONTIGS --> eHybrid found but not at end of contig... most likely too restrictive to find both primers')
            elif 2 in cases:
                short_file.write(
                    '\n' + '2:    ONE PRIMER --> eHybrid is not found in the same position as the primers found')
            elif 4 in cases:
                short_file.write(
                    '\n' + '4:    DIFF CONTIG --> eHybrid is not found in the same position as the primers found')
            elif 6 in cases:
                short_file.write('\n' + '6:    ONE PRIMER --> eHybrid was found but not of long enough length')
            elif 10 in cases:
                short_file.write('\n' + '10:   DIFF CONTIG --> eHybrid was found but not of long enough length')
            elif gene_name in ehyb_only:
                short_file.write('\n' + 'Also checked using ehyb only and not found.')
            else:
                short_file.write('\n' + 'Not found with a high enough BSR')
            # short_file.close()

        del all_hsp
        gc.collect()

        # short_file = open(result_file, "a")
        for gene_name in false_positive:
            short_file.write('\n' + '\n' + 'F+: '+ gene_name)
            for lo_hsp in cgf_predictions_dict[gene_name]:
                if lo_hsp != None:
                    for hsp in lo_hsp:
                        if hsp.both_primers_found and hsp.contig and hsp.pcr_distance:
                            short_file.write('\n' + 'BOTH PRIMERS-SAME CONTIG-VERY LIKELY TO BE TRUE POSITIVE')
                        elif hsp.both_primers_found and not hsp.contig and hsp.location:
                            short_file.write('\n' + 'BOTH PRIMERS-DIFF CONTIG-')
                        elif not hsp.both_primers_found and not hsp.contig and hsp.location:
                            short_file.write('\n' + 'ONE PRIMER-')
                        elif gene_name in ehyb_only:
                            short_file.write('\n' + 'EHYB ONLY')
        short_file.close()


def create_gene_results(forward_primers, reverse_primers, amplicon_sequences, db_directory, per_gene_file, lab_binary_results, ehyb_only):

    file = open(lab_binary_results, "r")

    f_primer_dict = CGFPrediction.create_primer_dict(forward_primers)
    r_primer_dict = CGFPrediction.create_primer_dict(reverse_primers)

    # bsr
    f_bs_primer_dir = "/home/sfisher/Sequences/BSR/f_primers/"
    r_bs_primer_dir = "/home/sfisher/Sequences/BSR/r_primers/"
    amp_bs_dir = "/home/sfisher/Sequences/BSR/amp_seq/"
    max_f_bits_dict = CGFPrediction.max_bs(f_bs_primer_dir)
    max_r_bits_dict = CGFPrediction.max_bs(r_bs_primer_dir)
    max_amp_bits_dict = CGFPrediction.max_bs(amp_bs_dir)

    lines = file.readlines()

    gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307',
                 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
                 'cj0736',
                 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
                 'cj1334',
                 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
                 'cj1727c']

    file_gene_dict = {}
    for line in lines:
        genes_found_list = []
        count = -2
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

    cases_per_gene_dict = defaultdict(list)
    num_f_pos = defaultdict(int)
    num_f_neg = defaultdict(int)
    for file_path in files_paths:
        file_name = file_path.partition(db_directory + "/")[2]
        cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
                                             amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict,
                                             True)
        cgf_predictions_dict = cgf_predictions[0]
        all_hsp = cgf_predictions[1]

        genes_expected = file_gene_dict[file_name]
        genes_found = [key for key in cgf_predictions_dict]

        false_positive = set(genes_found) - set(genes_expected)
        false_negative = set(genes_expected) - set(genes_found)

        # if len(false_positive) > 0 or len(false_negative) > 0:
            # print_file.write('\n' + file_name)
            # print_file.write('\n' + '\n' + file_name)

        binary_tree_attrib = ['both_primers_found', 'contig', 'ehybrid', 'ehybrid', 'ehybrid', 'location', None,
                              'epcr', None, 'location']
        for i in range(10, 16):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(16, 'valid')
        for i in range(17, 23):
            binary_tree_attrib.insert(i, None)
        for i in range(23, 33):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(33, 'snp')
        for i in range(34, 68):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(68, 'pcr_distance')
        for i in range(69, 139):
            binary_tree_attrib.insert(i, None)

        for gene_name in false_negative:
            num_f_neg[gene_name] = num_f_neg[gene_name] + 1
            cases = defaultdict(int)
            hi_bsr = 0
            for hsp in all_hsp[gene_name]:
                if hsp.bsr < MIN_BSR:
                    continue
                if hsp.name in gene_name:
                    result = CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0)
                        # print_file.write('\n' + hsp.name, "on contig", hsp.contig_name, "failed because hsp bsr is", hsp.bsr)
                    try:
                        if hsp.partner != None:
                            if hsp.partner.bsr < MIN_PARTNER_BSR:
                                continue
                    except:
                        continue
                    if hsp.bsr > hi_bsr:
                        hi_bsr = hsp.bsr
                    # print_file.write('\n' + gene_name + ' found with bsr: ' + str(hsp.bsr))
                    cases[result] += 1
                        # print_file.write('\n' + '\n')
            # print_file.write('\n' + '\n' + 'Cases and their occurences for: ' + gene_name + ' ' + str(cases))
            # print_file.write('\n' + 'Largest bsr: ' + str(hi_bsr) + '\n \n')

            # cases_per_gene_dict[gene_name].append('F+: ' + str(false_positive))
            # cases_per_gene_dict[gene_name].append('F-: ' + str(false_negative))
            # cases_per_gene_dict[gene_name].append('Num/Case: ' + str(cases))
            # cases_per_gene_dict[gene_name].append('Highest BSR ' + str(hi_bsr))
            if 67 in cases:
                cases_per_gene_dict[gene_name].append('F-  67:   SAME CONTIG --> SNP on 3 prime end')
            elif 138 in cases:
                cases_per_gene_dict[gene_name].append('F-  138:  SAME CONTIG --> Dist between primers not correct')
            elif 34 in cases:
                cases_per_gene_dict[gene_name].append('F-  34:   SAME CONTIG --> Primers are not facing each other')
            elif 3 in cases:
                cases_per_gene_dict[gene_name].append('F-  3:    SAME CONTIG --> eHybrid is not found in the same position as the primers found')
            elif 8 in cases:
                cases_per_gene_dict[gene_name].append('F-  8:    SAME CONTIG --> eHybrid was found but not of long enough length')
            elif 12 in cases:
                cases_per_gene_dict[gene_name].append('F-  12:   ONE PRIMER --> eHybrid found but not at end of contig... most likely too restrictive to find both primers.')
            elif 20 in cases:
                cases_per_gene_dict[gene_name].append('F-  20:   DIFF CONTIGS --> eHybrid found but not at end of contig... most likely too restrictive to find both primers')
            elif 2 in cases:
                cases_per_gene_dict[gene_name].append('F-  2:    ONE PRIMER --> eHybrid is not found in the same position as the primers found')
            elif 4 in cases:
                cases_per_gene_dict[gene_name].append('F-  4:    DIFF CONTIG --> eHybrid is not found in the same position as the primers found')
            elif 6 in cases:
                cases_per_gene_dict[gene_name].append('F-  6:    ONE PRIMER --> eHybrid was found but not of long enough length')
            elif 10 in cases:
                cases_per_gene_dict[gene_name].append('F-  10:   DIFF CONTIG --> eHybrid was found but not of long enough length')
            elif gene_name in ehyb_only:
                cases_per_gene_dict[gene_name].append('F-  : Also checked using ehyb only and not found.')
            else:
                cases_per_gene_dict[gene_name].append('F-  : Not found with a high enough BSR')
            # short_file.close()

        del all_hsp
        gc.collect()

        # short_file = open(result_file, "a")
        for gene_name in false_positive:
            num_f_pos[gene_name] = num_f_pos[gene_name] + 1
            for lo_hsp in cgf_predictions_dict[gene_name]:
                if lo_hsp != None:
                    for hsp in lo_hsp:
                        if hsp.both_primers_found and hsp.contig and hsp.pcr_distance:
                            cases_per_gene_dict[gene_name].append('F+  BOTH PRIMERS-SAME CONTIG-VERY LIKELY TO BE TRUE POSITIVE')
                        elif hsp.both_primers_found and not hsp.contig and hsp.location:
                            cases_per_gene_dict[gene_name].append('F+  BOTH PRIMERS-DIFF CONTIG-')
                        elif not hsp.both_primers_found and not hsp.contig and hsp.location:
                            cases_per_gene_dict[gene_name].append('F+  ONE PRIMER-')
                        elif gene_name in ehyb_only:
                            cases_per_gene_dict[gene_name].append('F+  EHYB ONLY')

    gene_file = open(per_gene_file, "a")
    for key in cases_per_gene_dict.keys():
        gene_file.write('\n \n' + key)
        gene_file.write('\n' + '# F+'+ str(num_f_pos[key]))
        gene_file.write('\n' + '# F-' + str(num_f_neg[key]))
        for val in cases_per_gene_dict[key]:
            gene_file.write('\n' + val)
    gene_file.close()

def create_long_causation(forward_primers, reverse_primers, amplicon_sequences, db_directory, long_file, med_file_f_neg, med_file_f_pos, lab_binary_results, ehyb_only):
    # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
    file = open(lab_binary_results, "r")

    f_primer_dict = CGFPrediction.create_primer_dict(forward_primers)
    r_primer_dict = CGFPrediction.create_primer_dict(reverse_primers)

    # bsr
    f_bs_primer_dir = "/home/sfisher/Sequences/BSR/f_primers/"
    r_bs_primer_dir = "/home/sfisher/Sequences/BSR/r_primers/"
    amp_bs_dir = "/home/sfisher/Sequences/BSR/amp_seq/"
    max_f_bits_dict = CGFPrediction.max_bs(f_bs_primer_dir)
    max_r_bits_dict = CGFPrediction.max_bs(r_bs_primer_dir)
    max_amp_bits_dict = CGFPrediction.max_bs(amp_bs_dir)

    lines = file.readlines()

    gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307',
                 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
                 'cj0736',
                 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
                 'cj1334',
                 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
                 'cj1727c']

    file_gene_dict = {}
    for line in lines:
        genes_found_list = []
        count = -2
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
        print_file = open(long_file, "a")
        file_name = file_path.partition(db_directory + "/")[2]
        cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
                                             amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict, True)
        cgf_predictions_dict = cgf_predictions[0]
        all_hsp = cgf_predictions[1]

        genes_expected = file_gene_dict[file_name]
        genes_found = [key for key in cgf_predictions_dict]

        false_positive = set(genes_found) - set(genes_expected)
        false_negative = set(genes_expected) - set(genes_found)

        if len(false_positive) > 0 or len(false_negative) > 0:
            # print_file.write('\n' + file_name)
            print_file.write('\n' + '\n' + file_name+ '\n')
            # print_file.write('\n' + 'genes found'+ str(genes_found))
            # print_file.write('\n' + 'genes expected'+ str(genes_expected))
        if len(false_positive) > 0:
            print_file.write('\n' + 'false positive'+ str(false_positive))
        if len(false_negative) > 0:
            print_file.write('\n' + 'false negative'+ str(false_negative) + '\n')

        binary_tree_attrib = ['both_primers_found', 'contig', 'ehybrid', 'ehybrid', 'ehybrid', 'location', None,
                       'epcr', None, 'location']
        for i in range(10, 16):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(16, 'valid')
        for i in range(17, 23):
            binary_tree_attrib.insert(i, None)
        for i in range(23, 33):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(33, 'snp')
        for i in range(34, 68):
            binary_tree_attrib.insert(i, None)
        binary_tree_attrib.insert(68, 'pcr_distance')
        for i in range(69, 139):
            binary_tree_attrib.insert(i, None)

        for gene_name in false_negative:
            cases = defaultdict(int)
            hi_bsr = 0
            for hsp in all_hsp[gene_name]:
                if hsp.name in gene_name:
                    result = CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0)
                    if hsp.bsr < MIN_BSR:
                        continue
                        # print_file.write('\n' + hsp.name, "on contig", hsp.contig_name, "failed because hsp bsr is", hsp.bsr)
                    try:
                        if hsp.partner != None:
                            if hsp.partner.bsr < MIN_PARTNER_BSR:
                                continue
                    except:
                        continue
                    if hsp.bsr > hi_bsr:
                        hi_bsr = hsp.bsr
                    print_file.write('\n' + gene_name + ' found with bsr: '+ str(hsp.bsr))
                    cases[result] += 1
                    if result == 6:
                        print_file.write('\n' + 'CASE 6: ' + gene_name + ' not found using eCGF with both primers so one primer'
                                                'searched for with EHYBRID ' + '(' + str(hsp.ehybrid) +') not found of long enough length (len=' + str(hsp.length) + ")")
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                    elif result == 10:
                        print_file.write('\n' + 'CASE 10: ' + gene_name + ' Found primers on DIFFERENT contigs and '
                                                'EHYBRID: ' + '(' + str(hsp.ehybrid) +')' + ' found but not ' \
                                                                                                'of long enough length')
                        print_file.write('\n' + 'Contigs:' + hsp.contig_name + hsp.partner.contig_name)
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 12:
                        print_file.write('\n' + 'CASE 12:' + gene_name + ' not found using eCGF with both primers so one primer'
                                                'searched for with EHYBRID found ' + '(' + str(hsp.ehybrid) +')' + 'but not located at the end of a contig. (hsp location is ')
                        if hsp.strand:
                            print_file.write(str(hsp.end_dist) + 'from end (on leading strand)')
                        elif not hsp.strand:
                            print_file.write(str(hsp.end_dist) + 'from end (on lagging strand)')
                        print_file.write('\n' + 'Most likely too restrictive to find both primers...')
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                    elif result == 2:
                        print_file.write('\n' + "CASE 2: " + gene_name + ' Not found using eCGF with both primers so one primer'
                                                'searched for with EHYBRID'+ '(' + str(hsp.ehybrid) +')' + 'not found at all!')
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                    elif result == 3:
                        print_file.write('\n' + 'CASE 3: '+ gene_name + ' Both primers found on the same contig but EHYBRID'+ '(' + str(hsp.ehybrid) +')' + 'not found at all!')
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 4:
                        print_file.write('\n' + "CASE 4: "+ gene_name + ' Found primers on DIFFERENT contigs and'
                                                'EHYBRID: ' + '(' + str(hsp.ehybrid) +')' + 'not found at all!')
                        print_file.write('\n' + 'Contigs:' + hsp.contig_name + hsp.partner.contig_name)
                        print_file.write('\n' + 'strand'+ str(hsp.strand))
                        print_file.write('\n' + 'sbjct'+ hsp.sbjct)
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 8:
                        print_file.write('\n' + 'CASE 8: '+ gene_name + ' Both primers found on the same contig '
                                                                        'with EHYBRID'+ '(' + str(hsp.ehybrid) +')' +
                                         'not found b/c length too short: len=' + str(hsp.length))
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 20:
                        print_file.write('\n' + 'CASE 20: '+ gene_name + ' Found primers on DIFFERENT contigs and'
                                                'EHYBRID: ' + '(' + str(hsp.ehybrid) +')' + ' found but at the end of a contig. (hsp location is ')
                        if hsp.strand:
                            print_file.write(str(hsp.end_dist) + ' from end (leading strand)')
                        elif not hsp.strand:
                            print_file.write(str(hsp.end_dist) + ' from end (lagging strand)')
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 11 or result == 15 or result == 137:
                        print_file.write('\n' + 'FALSE')
                    elif result == 34:
                        print_file.write('\n' + "CASE 34: " + gene_name + " Both primers on the same contig with EHYBRID"
                                         + '(' + str(hsp.ehybrid) +')' + "found but epcr not found" + '(' + str(hsp.epcr) +')' +
                                         'because strands are not facing each other (valid =' + hsp.valid + ')')
                        print_file.write('\n' + 'both primers found:'+ str(hsp.both_primers_found) + 'with contig:'+ str(hsp.contig))
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 138:
                        print_file.write('\n' + "CASE 138: " + gene_name + " Both primers on the same contig with EHYBRID"
                                         + '(' + str(hsp.ehybrid) +')' + "found but epcr not found" + '(' + str(hsp.epcr) +')' +
                                         'b/c primers are not the correct distance from each other')
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    elif result == 67:
                        print_file.write('\n' + "CASE 67: " + gene_name + " Both primers on the same contig with EHYBRID"
                                         + '(' + str(hsp.ehybrid) +')' + "found but epcr not found" + '(' + str(hsp.epcr) +')' +
                                         'because of SNP on "3'" end")
                        print_file.write('\n' + 'strand '+ str(hsp.strand))
                        print_file.write('\n' + 'query '+ hsp.query)
                        print_file.write('\n' + 'sbjct '+ hsp.sbjct)
                        print_file.write('\n' + 'f_primer seq: ' + str(f_primer_dict[gene_name]))
                        print_file.write('\n' + 'r_primer seq: ' + str(r_primer_dict[gene_name]))
                        print_file.write('\n' + 'r_comp_f_pri: ' + str(f_primer_dict[gene_name].reverse_complement()))
                        print_file.write('\n' + 'r_comp_r_pri: ' + str(r_primer_dict[gene_name].reverse_complement()))
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + 'HSP partners attr:' + str(hsp.partner.__dict__))
                    if gene_name in ehyb_only:
                        print_file.write('\n' + gene_name + '--- also not found using ehyb only')
                    # else:
                    #     print_file.write('\n' + str(result))
                    #     assert False
                    print_file.write('\n' + '\n')
            print_file.write('\n' + '\n' + 'Cases and their occurences for: '+ gene_name + ' ' + str(cases))
            print_file.write('\n' + 'Largest bsr: ' + str(hi_bsr) + '\n \n')

            short_file = open(med_file_f_neg, "a")
            short_file.write('\n' + '\n' + file_name)
            short_file.write('\n' + '\n' + gene_name)
            short_file.write('\n' + 'F+: ' + str(false_positive))
            short_file.write('\n' + 'F-: ' + str(false_negative))

            short_file.write('\n' + 'Num/Case: ' + str(cases))
            short_file.write('\n' + 'Highest BSR ' + str(hi_bsr))
            #TODO: determine where the highest BSR came from.
            if 67 in cases:
                short_file.write('\n' + '67:   SAME CONTIG --> SNP on 3 prime end')
            if 138 in cases:
                short_file.write('\n' + '138:  SAME CONTIG --> Dist between primers not correct')
            if 34 in cases:
                short_file.write('\n' + '34:   SAME CONTIG --> Primers are not facing each other')
            if 3 in cases:
                short_file.write('\n' + '3:    SAME CONTIG --> eHybrid is not found in the same position as the primers found')
            if 8 in cases:
                short_file.write('\n' + '8:    SAME CONTIG --> eHybrid was found but not of long enough length')
            if 12 in cases:
                short_file.write('\n' + '12:   ONE PRIMER --> eHybrid found but not at end of contig... most likely too restrictive to find both primers.')
            if 20 in cases:
                short_file.write('\n' + '20:   DIFF CONTIGS --> eHybrid found but not at end of contig... most likely too restrictive to find both primers')
            if 2 in cases:
                short_file.write('\n' + '2:    ONE PRIMER --> eHybrid is not found in the same position as the primers found')
            if 4 in cases:
                short_file.write('\n' + '4:    DIFF CONTIG --> eHybrid is not found in the same position as the primers found')
            if 6 in cases:
                short_file.write('\n' + '6:    ONE PRIMER --> eHybrid was found but not of long enough length')
            if 10 in cases:
                short_file.write('\n' + '10:   DIFF CONTIG --> eHybrid was found but not of long enough length')
            if gene_name in ehyb_only:
                short_file.write('\n' + '--    Also not found using ehyb only.')
            short_file.close()

        del all_hsp
        gc.collect()

        #TODO: debug... F+ not showing up!

        short_file = open(med_file_f_pos, "a")
        for gene_name in false_positive:
            short_file.write('\n' + '\n' + file_name)
            # if gene_name in ehyb_only:
            #     print_file.write('\n' + gene_name + ' found using ehyb only')
            # else:
            for lo_hsp in cgf_predictions_dict[gene_name]:
            # print_file.write('\n' + type(lo_hsp))
                if lo_hsp != None:
                    if gene_name in ehyb_only:
                        print_file.write('\n' + gene_name + ' found using ehyb only')
                    for hsp in lo_hsp:
                        print_file.write('\n' + hsp.name)
                        if hsp.both_primers_found and hsp.contig and hsp.pcr_distance:
                            short_file.write('\n' + 'BOTH PRIMERS-SAME CONTIG-VERY LIKELY TO BE TRUE POSITIVE')
                            print_file.write('\n' + 'BOTH PRIMERS-SAME CONTIG-VERY LIKELY TO BE TRUE POSITIVE')
                        elif hsp.both_primers_found and not hsp.contig and hsp.location:
                            short_file.write('\n' + 'BOTH PRIMERS-DIFF CONTIG-')
                            print_file.write('\n' + 'BOTH PRIMERS-DIFF CONTIG-')
                        elif not hsp.both_primers_found and not hsp.contig and hsp.location:
                            short_file.write('\n' + 'ONLY ONE PRIMER FOUND-')
                        elif gene_name in ehyb_only:
                            short_file.write('\n' + 'ONLY FOUND USING EHYB')
                            short_file.write('\n' + '% id ' + str(hsp.identities / hsp.length))
                            print_file.write('\n' + 'ONLY FOUND USING EHYB')
                            print_file.write('\n' + '% id ' + str(hsp.identities / hsp.length))
                        print_file.write('\n' + 'bsr' + str(hsp.bsr))
                        print_file.write('\n' + 'strand'+ str(hsp.strand))
                        print_file.write('\n' + 'query         '+ hsp.query)
                        print_file.write('\n' + 'sbjct         '+ hsp.sbjct)
                        print_file.write('\n' + 'f_primer seq: ' + str(f_primer_dict[gene_name]))
                        print_file.write('\n' + 'r_primer seq: ' + str(r_primer_dict[gene_name]))
                        print_file.write('\n' + 'r_comp_f_pri: ' + str(f_primer_dict[gene_name].reverse_complement()))
                        print_file.write('\n' + 'r_comp_r_pri: ' + str(r_primer_dict[gene_name].reverse_complement()))
                        print_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        print_file.write('\n' + '\n')
                        short_file.write('\n' + hsp.name)
                        short_file.write('\n' + str(hsp.bsr))
                        short_file.write('\n' + 'strand'+ str(hsp.strand))
                        short_file.write('\n' + 'query         '+ hsp.query)
                        short_file.write('\n' + 'sbjct         '+ hsp.sbjct)
                        short_file.write('\n' + 'f_primer seq: ' + str(f_primer_dict[gene_name]))
                        short_file.write('\n' + 'r_primer seq: ' + str(r_primer_dict[gene_name]))
                        short_file.write('\n' + 'r_comp_f_pri: ' + str(f_primer_dict[gene_name].reverse_complement()))
                        short_file.write('\n' + 'r_comp_r_pri: ' + str(r_primer_dict[gene_name].reverse_complement()))
                        short_file.write('\n' + 'HSP attr: ' + str(hsp.__dict__))
                        short_file.write('\n' + '\n')
                else:
                    print_file.write('\n' + 'no hsp found')
        print_file.close()
        short_file.close()




if __name__ == "__main__":

    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes/debug_f_neg"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"
    # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
    # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
    lab_binary_results = "/home/sfisher/Sequences/11168_test_files/cgf40_v2.txt"

    #Output text files
    table_file = "/home/sfisher/Sequences/11168_test_files/tables/1_ehybpercid-85.txt"
    short_file = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/30_ehyb_short_expl.txt"
    gene_file = '/home/sfisher/Sequences/11168_test_files/eCGF_causation/30_ehyb_expl_per_gene.txt'
    long_file = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/31_ehyb_full_expl_percid85.txt"
    med_file_f_neg = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/31_ehyb_short_expl_f_neg_percid85.txt"
    med_file_f_pos = "/home/sfisher/Sequences/11168_test_files/eCGF_causation/31_ehyb_short_expl_f_pos_percid85.txt"

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