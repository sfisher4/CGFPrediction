import gc
import os
from collections import defaultdict
from memory_profiler import profile
from collections import ChainMap
import pandas as pd

from cgf_pred import CGFPrediction

# @profile
# def create_binary_table(forward_primers, reverse_primers, amplicon_sequences, db_directory, table_file, lab_binary_results, ehyb_only):
#     # file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
#     # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
#
#     gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
#                  'cj0307', 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728',
#                  'cj0733', 'cj0736', 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294',
#                  'cj1324', 'cj1329','cj1334','cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552',
#                  'cj1585', 'cj1679', 'cj1721','cj1727c']
#
#     file = open(lab_binary_results, "r")
#     lines = file.readlines()
#     file_gene_dict = {}
#     for line in lines:
#         genes_found_list = []
#         count = -2
#         print(len(line.split()))
#         for word in line.split():
#             print(word)
#             if word == '1':
#                 genes_found_list.append(gene_list[count])
#             count += 1
#         # print(genes_found_list)
#         file_gene_dict[line.split(None, 1)[0]] = genes_found_list
#         # print(file_gene_dict)
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
#         cgf_predictions = CGFPrediction.ecgf(forward_primers, reverse_primers, file_path,
#                                              amplicon_sequences, file_gene_dict)
#         cgf_predictions_dict = cgf_predictions[0]
#         pred_exception_dict = cgf_predictions[2]
#         single_primer_dict = cgf_predictions[3]
#         assert file_name in file_gene_dict, "error: The db directory contains a file that does not have validation results in 'file'"
#         genes_expected = file_gene_dict[file_name]
#         genes_found_without_exceptions = [key[6:] for key in cgf_predictions_dict.keys()]
#
#         gene_exceptions_found = [key[6:] for key in pred_exception_dict.keys()]
#         single_primers_found = [key[6:] for key in single_primer_dict.keys()]
#         genes_found = ChainMap({}, genes_found_without_exceptions, gene_exceptions_found, single_primers_found)
#         print(genes_found)
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
#         #PANDAS
#         df_dict = {}
#         for gene in gene_list:
#             if gene in false_positive and \
#                 (gene in gene_exceptions_found or
#                          gene in single_primers_found):
#                 df_dict[gene] = "+0.5"
#             elif gene in false_positive:
#                 df_dict[gene] = "+1"
#             elif gene in false_negative:
#                 df_dict[gene] = "-1"
#             else:
#                 df_dict[gene] = "0"
#         df = pd.DataFrame(df_dict, index=[file_name[:-6]])
#         df.to_csv(table_file, sep=' ', mode='a', header=None)
#
#     #F+/F- Calculation
#     t_file = open(table_file)
#     t_lo_f_pos = {}
#     t_lo_f_neg = {}
#     t_lines = t_file.readlines()
#     for line in t_lines:
#         count = 0
#         for word in line.split():
#             # print(word)
#             if word == "+1":
#                 try:
#                     t_lo_f_pos[gene_list[count]] += 1
#                 except:
#                     t_lo_f_pos[gene_list[count]] = 1
#                 count += 1
#             elif word == "-1":
#                 try:
#                     t_lo_f_neg[gene_list[count]] += 1
#                 except:
#                     t_lo_f_neg[gene_list[count]] = 1
#                 count += 1
#             elif word == "0":
#                 if gene_list[count] not in t_lo_f_neg:
#                     t_lo_f_neg[gene_list[count]] = 0
#                 if gene_list[count] not in t_lo_f_pos:
#                     t_lo_f_pos[gene_list[count]] = 0
#                 count += 1
#             elif word == "+0.5":
#                 try:
#                     t_lo_f_pos[gene_list[count]] += 0.5
#                 except:
#                     t_lo_f_pos[gene_list[count]] = 0.5
#                 count += 1
#
#     t_file.close()
#
#     with open(table_file, "a") as myfile:
#         myfile.write('\n' + "\n" + "# of F+: " + str(t_lo_f_pos) + "\n")
#         myfile.write('\n' + "\n" + "# of F-: " + str(t_lo_f_neg) + "\n")


def new_validation(lab_binary_results, ecgf_binary_results, validation_file, val_total_file, lab_result_delimiter):
    gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307', 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728',
                 'cj0733', 'cj0736', 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294',
                 'cj1324', 'cj1329', 'cj1334', 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552',
                 'cj1585', 'cj1679', 'cj1721', 'cj1727c']

    with open(lab_binary_results, "r") as lab_file:
        next(lab_file)
        lines = lab_file.readlines()
        file_gene_exp_dict = {}
        for line in lines:
            genes_expected = []
            for word in line.split(lab_result_delimiter):
                # print('word', word)
                if word == '1' or word == '1\n':
                    genes_expected.append(1)
                elif word == '0' or word == '0\n':
                    genes_expected.append(0)
                else:
                    print('!!!', word)
            genome_name = line.split(',', 1)[0]
            genome_name = genome_name.replace('.fasta', "")
            file_gene_exp_dict[genome_name] = genes_expected
            # print(len(genes_expected))
        print('genes expected', file_gene_exp_dict)

    with open(ecgf_binary_results, "r") as ecgf_file:
        next(ecgf_file)
        ecgf_lines = ecgf_file.readlines()
        file_gene_found_dict = {}
        for ecgf_line in ecgf_lines:
            genes_found = []
            for word in ecgf_line.split():
                # print(word)
                if word == '1' or word == '1,' or word == '[1,' or word == '1]':
                    genes_found.append(1)
                elif word == '0' or word == '0,' or word == '[0,' or word == '0]':
                    genes_found.append(0)
                elif word == '2' or word == '2,' or word == '[2,' or word == '2]':
                    genes_found.append(2)
            # print(len(genes_found))
            file_gene_found_dict[ecgf_line.split(None, 1)[0]] = genes_found
            # print(len(file_gene_found_dict[ecgf_line.split(None, 1)[0]]))
        print('genes found', file_gene_found_dict)

    dict_compare_results = {}
    total_false_positive = {gene: 0 for gene in gene_list}
    total_false_negatives = {gene: 0 for gene in gene_list}
    total_uncertain_f_pos = {gene: 0 for gene in gene_list}
    total_uncertain_f_neg = {gene: 0 for gene in gene_list}

    with open(validation_file, "a") as print_file:
        for genome in file_gene_exp_dict.keys():
            try:
                lo_found = file_gene_found_dict[genome]
            except:
                print('Whoops, lab data has genome:', genome, 'which is not contained in eCGF data')
                continue

            lo_expected = file_gene_exp_dict[genome]
            # print(len(lo_expected))
            count = 0
            lo_compare = []
            for expect in lo_expected:
                found = lo_found[count]
                if expect == found:
                    lo_compare.append('0')
                elif expect == 1 and found == 0:
                    lo_compare.append('-1')
                    total_false_negatives[gene_list[count]] += 1
                elif expect == 1 and found == 2:
                    lo_compare.append('-0.5')
                    total_uncertain_f_neg[gene_list[count]] += 1
                elif expect == 0 and found == 1:
                    lo_compare.append('+1')
                    total_false_positive[gene_list[count]] += 1
                elif expect == 0 and found == 2:
                    lo_compare.append('+0.5')
                    total_uncertain_f_pos[gene_list[count]] += 1
                else:
                    print(expect)
                    print(found)
                    assert False == True
                count += 1

            dict_compare_results[genome] = lo_compare
            print_file.write(genome + str(dict_compare_results[genome]) + '\n')
        print_file.write('\n' + 'total F+ / gene'+ str(total_false_positive))
        print_file.write('\n' + 'total F- / gene'+ str(total_false_negatives))
        print_file.write('\n' + 'total uncertain F+ / gene' + str(total_uncertain_f_pos))
        print_file.write('\n' + 'total uncertain F- / gene' + str(total_uncertain_f_neg))

        df = pd.DataFrame({
            'F+': total_false_positive,
            'F-': total_false_negatives,
            'Total "2"s when lab data DOES NOT find the gene' : total_uncertain_f_pos,
            'Total "2"s when lab data DOES find the gene': total_uncertain_f_neg
        })
        df.to_csv(total_validation_results, mode='a', sep=" ", )

#     # TODO: remove from program or make seperate case for debugging in the future!
# def fourth_case_check(ehyb_pos, result_dict, f_primer_dict, r_primer_dict, amp_dict, genome_name,
#                   lab_results_dict):
#
#     with open("/home/sfisher/eCGF/all_genomes/Results/false_neg_causes_per_genome", "a") as file:
#         ehyb_pos_names = [hsp.name for hsp in ehyb_pos]
#         file.write('\n \n Genome name: ' + str(genome_name))
#         file.write('\n Genes not found using eCGF but found using eHYB ' + str(ehyb_pos_names))
#         file.write('\n Number of genes found in eHYB only ' + str(len(ehyb_pos)))
#         file.write('\n Number of genes found using eCGF (Probability of being TRULY +ve is ~high) ' + str(
#             len(result_dict.keys())))
#         file.write('\n Number of genes not found at all (Probability of being TRULY -ve is high) ' + str(
#             abs(40 - len(result_dict.keys()) - len(ehyb_pos))))
#         file.write(
#             '\n Below is more info on the genes found using eHYB only (was -ve in eCGF but +ve in eHYB): ')
#         if ".fasta" not in genome_name:
#             lab_results = lab_results_dict[genome_name + '.fasta']
#         else:
#             lab_results = lab_results_dict[genome_name]
#
#         print('lab results', lab_results)
#         for hsp in ehyb_pos:
#             hsp_name = hsp.name
#             perc_id = hsp.identities / hsp.length
#             lab_result = "Found" if hsp.name in lab_results else "Not Found"
#             file.write('\n \n' + hsp_name)
#             file.write('\n Lab Result: ' + lab_result)
#             file.write('\n contig: ' + hsp.contig_name)
#             file.write('\n % id ' + str(perc_id))
#             file.write('\n BSR: ' + str(hsp.bsr))
#             file.write('\n End Dist: ' + str(hsp.end_dist))
#             file.write('\n qcov (% of query seqn that overlaps sbjct seqn) ' + str(
#                 len(hsp.query) / (len(amp_dict[hsp.name]) + 1)))
#             file.write('\n Entire Query Seq      ' + str(len(amp_dict[hsp.name]) + 1) + " 1 " + str(
#                 len(amp_dict[hsp.name])) + " " + str(amp_dict[hsp.name]))
#             file.write(
#                 '\n Matching Query Seq    ' + str(len(hsp.query)) + " " + str(hsp.query_start) + " " + str(
#                     hsp.query_end) + " " * (hsp.query_start + (
#                 7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(
#                     str(hsp.query_end)))) + hsp.query)
#             file.write(
#                 '\n Matching Sbjct Seq    ' + str(len(hsp.sbjct)) + " " + str(hsp.query_start) + " " + str(
#                     hsp.query_end) + " " * (hsp.query_start + (
#                 7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(
#                     str(hsp.query_end)))) + hsp.sbjct)
#             file.write('\n Forward primer         ' + str(len(f_primer_dict[hsp_name])) + "       " + str(
#                 f_primer_dict[hsp_name]))
#             file.write('\n Reverse primer         ' + str(len(r_primer_dict[hsp_name])) + " " * (
#             len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]) + 7) + str(
#                 amp_dict[hsp.name][len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]):]))
#
#             # if (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
#             file.write('\n Forward Primer qcov: ' + str(
#                 (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(
#                     f_primer_dict[hsp_name])) + " (CUTOFF: " + str((QCOV_HSP_PERC / 100)) + " )")
#             # if (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
#             file.write('\n Reverse Primer qcov: ' + str(
#                 (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(
#                     r_primer_dict[hsp_name])) + " (CUTOFF:  " + str((QCOV_HSP_PERC / 100)) + " )")
#             # if perc_id * 100 < PERC_ID_CUTOFF:
#             file.write('\n Perc ID for ehyb in eCGF: ' + str(perc_id) + " (CUTOFF: " + str(
#                 (PERC_ID_CUTOFF / 100)) + " )")
#
# def per_genome_false_neg_check(lab_binary_results, ecgf_binary_results, amp_sequences, database):
#     gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
#                  'cj0307', 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728',
#                  'cj0733', 'cj0736', 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294',
#                  'cj1324', 'cj1329', 'cj1334', 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552',
#                  'cj1585', 'cj1679', 'cj1721', 'cj1727c']
#
#     count = 0
#     lo_genes_found = []
#     lo_genes_found_two = []
#     for result in ecgf_binary_results:
#         if result == 1:
#             lo_genes_found.append(gene_list[count])
#             count += 1
#         elif result == 2:
#             lo_genes_found_two.append(gene_list[count])
#             count += 1
#
#
#     ehyb_blast_hsps = CGFPrediction.create_blastn_object(amp_sequences, database, False, 70)
#     lo_hsp_ehybrid = CGFPrediction.ehyb(ehyb_blast_hsps)  # assigns ehybrid attributes to each hsp from amp vs db
#     ehybrid_pass = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == True]
#     ehyb_pos = [hsp for hsp in ehybrid_pass if hsp.name not in lo_genes_found]
#
#     with open(lab_binary_results, "r") as lab_file:
#         next(lab_file)
#         lines = lab_file.readlines()
#         file_gene_dict = {}
#         for line in lines:
#             genes_expected = []
#             for word in line.split(","):
#                 # print('word', word)
#                 if word == '1' or word == '1\n':
#                     genes_expected.append(1)
#                 elif word == '0' or word == '0\n':
#                     genes_expected.append(0)
#                     # else:
#                     # print('!!!', word)
#             genome_name = line.split(',', 1)[0]
#             # print(genome_name)
#             if '.fasta' not in genome_name:
#                 file_gene_dict[genome_name + '.fasta'] = genes_expected
#             else:
#                 file_gene_dict[genome_name] = genes_expected
#
#                 # print(len(genes_expected))
#                 # print('genes expected', file_gene_dict)
#
#     if os.path.basename(database) in file_gene_dict.keys():
#         fourth_case_check(ehyb_pos, result_dict, dict_f_primers, dict_r_primers, dict_amp,
#                           os.path.basename(database), file_gene_dict)


if __name__ == "__main__":

    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
    # db_directory = "/home/sfisher/Sequences/11168_test_files/debug_genes"

    # lab_binary_results = "/home/sfisher/Sequences/11168_test_files/cgf40_v2.txt"
    # ecgf_binary_results = "/home/sfisher/eCGF/Results_246_gnomes/without_exceptions/oct_20_results.txt"
    lab_binary_results = "/home/sfisher/eCGF/all_genomes/Results/nov_1_labcgf_1046genomes.txt"
    ecgf_binary_results = "/home/sfisher/eCGF/all_genomes/Results/nov_1_all_genomes"
    table_file = "/home/sfisher/Sequences/11168_test_files/tables/testing.txt"
    total_validation_results = "/home/sfisher/eCGF/all_genomes/Results/nov_1_false_sums"
    lab_result_delimiter = ','
    validation_file = "/home/sfisher/eCGF/all_genomes/Results/nov_1_validation_results.txt"

    ehyb_only = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264c', 'cj0297c', 'cj0298c',
                 'cj0307',
                 'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
                 'cj0736',
                 'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
                 'cj1334',
                 'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
                 'cj1727c']

    # create_binary_table(forward_primers, reverse_primers, amplicon_sequences, db_directory, table_file, lab_binary_results, ehyb_only)
    # ecgf_validation()
    new_validation(lab_binary_results, ecgf_binary_results, validation_file, total_validation_results, lab_result_delimiter)

