import unittest
import CGFPrediction
import os
import errno

class MyTestCase(unittest.TestCase):
    # def test_pcr_prediction(self):
    #     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    #     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    #     forward_dict = CGFPrediction.create_primer_dict(forward_primers)
    #     reverse_dict = CGFPrediction.create_primer_dict(reverse_primers)
    #     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     test_11168 = "/home/sfisher/Sequences/11168_complete_genome"
    #     result_dict = CGFPrediction.main(test_11168, forward_primers, reverse_primers, amplicon_sequences)
    #     self.assertEqual(len(result_dict), 1)
    #     result = result_dict["11168_complete_genome.fasta"]
    #     for lo_hsp in result:
    #         self.assertEqual(len(lo_hsp), 2)
    #         self.assertEqual(lo_hsp[0].name, lo_hsp[1].name)
    #         self.assertEqual(CGFPrediction.pcr_directly(lo_hsp[0], lo_hsp[1], amplicon_sequences, forward_dict, reverse_dict), True)
    #     self.assertEqual(len(result), 40)
    # #
    # # def test_cgf_prediction_binary(self):
    # #     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    # #     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    # #     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    # #     test_11168_cases = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    # #     cgf_prediction_dict = CGFPrediction.main(test_11168_cases, forward_primers, reverse_primers, amplicon_sequences)
    # #
    # #     files = [file for file in os.listdir(test_11168_cases) if file.endswith('.fasta')]
    # #     file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt")
    # #     lines = file.readlines()
    # #
    # #
    # #     file_list = []
    # #     for line in lines:
    # #         new_line = line.split()
    # #         file_list.append(new_line)
    # #     file.close()
    # #     print(file_list)
    # #
    # #     dict_result = {}
    # #     for result_list in file_list:
    # #         dict_result[result_list[0]] = result_list[1:]
    # #
    # #     for file in files:
    # #         lo_queries = cgf_prediction_dict[file]
    # #
    # #         gene_list = ['cj0008', 'cj0033', 'cj0035', 'cj0057', 'cj0177', 'cj0181', 'cj0264', 'cj0297c', 'cj0298c',
    # #                      'cj0307',
    # #                      'cj0421c', 'cj0483', 'cj0486', 'cj0566', 'cj0569', 'cj0570', 'cj0625', 'cj0728', 'cj0733',
    # #                      'cj0736',
    # #                      'cj0755', 'cj0860', 'cj0967', 'cj1134', 'cj1136', 'cj1141', 'cj1294', 'cj1324', 'cj1329',
    # #                      'cj1334',
    # #                      'cj1427c', 'cj1431c', 'cj1439', 'cj1550c', 'cj1551', 'cj1552', 'cj1585', 'cj1679', 'cj1721',
    # #                      'cj1727c']
    # #         genes_found = []
    # #         for quer in lo_queries:
    # #             genes_found.append(quer[0].name)
    # #             #
    # #             # try:
    # #             #     if quer[0].contig_name != quer[1].contig_name:
    # #             #         print(file, "primers on different contigs for", quer[0].name)
    # #             #         print('contig name of 1st', quer[0].contig_name)
    # #             #         print('contig name of 2nd', quer[1].contig_name)
    # #             # except IndexError:
    # #             #     print("only one primer found for", quer[0].name)
    # #         print('genes found', genes_found)
    # #
    # #         lo_genes_present = []
    # #         for gene in gene_list:
    # #             if gene in genes_found:
    # #                 lo_genes_present.append('1')
    # #             else:
    # #                 lo_genes_present.append('0')
    # #         # lo_tup_diff_contig_names = []
    # #         # for tup in lo_tup_diff_contig:
    # #         #     lo_tup_diff_contig_names.append(tup[0].name)
    # #         # print('lo hsp on diff contigs', lo_tup_diff_contig_names)
    # #         print(file, lo_genes_present)
    # #         print(dict_result[file])
    # #
    # #         for i in range(39):
    # #             if dict_result[file][i] != lo_genes_present[i]:
    # #                 print('Incorrect result:', gene_list[i])
    # #             # self.assertEqual(dict_result[file][i],lo_genes_present[i])
    #
    # def test_cgf_prediction_same_contig(self):
    #     #TODO: compare the expected results with the actual results using the names rather than binary
    #     #TODO: determine if the primers are on different contigs and if that is affecting the results!
    #
    #     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    #     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    #     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    #
    #     file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
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
    #     files = [file for file in os.listdir(db_directory) if file.endswith('.fasta')]
    #     files_paths = []
    #     for file in files:
    #         files_paths.append(os.path.abspath(db_directory) + '/' + file)
    #
    #     # create new folder for out files
    #     try:
    #         os.mkdir(db_directory + "/out_files")
    #     except OSError as exc:
    #         if exc.errno != errno.EEXIST:
    #             raise
    #         pass
    #
    #     for file_path in files_paths:
    #         print('FILE!!!', file_path)
    #         file_name = file_path.partition(db_directory + "/")[2]
    #         print(file_name)
    #         f_out_file_path = db_directory + "/out_files/" + "f_" + file_name.replace("fasta", "xml")
    #         r_out_file_path = db_directory + "/out_files/" + "r_" + file_name.replace("fasta", "xml")
    #         full_out_file_path = db_directory + "/out_files/" + "full_" + file_name.replace("fasta", "xml")
    #
    #         cgf_predictions = CGFPrediction.pcr_prediction(forward_primers, reverse_primers, file_path, f_out_file_path, r_out_file_path,
    #                                 amplicon_sequences, full_out_file_path)
    #         lo_predictions_names = []
    #         for pred in cgf_predictions[0]:
    #             lo_predictions_names.append(pred[0].name)
    #
    #         lo_diff_contig_genes = cgf_predictions[1]
    #         lo_single_primers = cgf_predictions[3]
    #         hybridized_results_same_contig = []
    #         for gene_name in cgf_predictions[2]:
    #             if "11168_" in gene_name:
    #                 tmp = gene_name[6:]
    #                 hybridized_results_same_contig.append(tmp)
    #
    #         genes_expected = file_gene_dict[file_name]
    #         genes_found = lo_predictions_names
    #         print('genes expected', genes_expected)
    #         print('hybridization results', hybridized_results_same_contig)
    #         print('genes found', genes_found)
    #         print('diff contig genes', lo_diff_contig_genes)
    #         false_positive = set(genes_found) - set(genes_expected) - set(lo_diff_contig_genes) - set(lo_single_primers)
    #         false_negative = set(genes_expected) - set(genes_found) - set(lo_diff_contig_genes) - set(lo_single_primers)
    #         hybridized_results_false_positive = set(hybridized_results_same_contig) - set(genes_expected) - set(lo_diff_contig_genes) - set(lo_single_primers)
    #         hybridized_results_false_negative = set(genes_expected) - set(hybridized_results_same_contig) - set(lo_diff_contig_genes) - set(lo_single_primers)
    #         print('false +ve', false_positive)
    #         print('false -ve', false_negative)
    #         print('hybridization false +ve', hybridized_results_false_positive)
    #         print('hybridization false -ve', hybridized_results_false_negative)
    #         print("\n")
    #         print("\n")
    #
    # def test_cgf_pred_diff_contigs(self):
    #
    #     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    #     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    #     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    #
    #     file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
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
    #     files = [file for file in os.listdir(db_directory) if file.endswith('.fasta')]
    #     files_paths = []
    #     for file in files:
    #         files_paths.append(os.path.abspath(db_directory) + '/' + file)
    #
    #     # create new folder for out files
    #     try:
    #         os.mkdir(db_directory + "/out_files")
    #     except OSError as exc:
    #         if exc.errno != errno.EEXIST:
    #             raise
    #         pass
    #
    #     for file_path in files_paths:
    #         print('FILE!!!', file_path)
    #         file_name = file_path.partition(db_directory + "/")[2]
    #         print(file_name)
    #         f_out_file_path = db_directory + "/out_files/" + "f_" + file_name.replace("fasta", "xml")
    #         r_out_file_path = db_directory + "/out_files/" + "r_" + file_name.replace("fasta", "xml")
    #         full_out_file_path = db_directory + "/out_files/" + "full_" + file_name.replace("fasta", "xml")
    #
    #         cgf_predictions = CGFPrediction.pcr_prediction(forward_primers, reverse_primers, file_path, f_out_file_path, r_out_file_path,
    #                                 amplicon_sequences, full_out_file_path)
    #         lo_predictions_names = []
    #         print('cgf predictions[0]', cgf_predictions[0])
    #         for quer in cgf_predictions[0]:
    #             for hsp in quer:
    #                 lo_predictions_names.append(hsp.name)
    #             # lo_predictions_names.append(quer[0].name)
    #
    #         lo_diff_contig_genes = cgf_predictions[1]
    #         lo_single_primers = cgf_predictions[3]
    #         hybridized_results_same_contig = []
    #         for gene_name in cgf_predictions[2]:
    #             if "11168_" in gene_name:
    #                 tmp = gene_name[6:]
    #                 hybridized_results_same_contig.append(tmp)
    #
    #         genes_expected = file_gene_dict[file_name]
    #         genes_found = lo_predictions_names
    #         print('genes expected', genes_expected)
    #         print('genes found', genes_found)
    #         print('hybridization genes found', hybridized_results_same_contig)
    #         print('diff contig genes', lo_diff_contig_genes)
    #         print('single primers found', lo_single_primers)
    #         false_positive = set(genes_found) - set(genes_expected) - set(lo_single_primers)
    #         false_negative = set(genes_expected) - set(genes_found) - set(lo_single_primers)
    #         hybridized_results_false_positive = set(hybridized_results_same_contig) - set(genes_expected) - set(lo_single_primers)
    #         hybridized_results_false_negative = set(genes_expected) - set(hybridized_results_same_contig) - set(lo_single_primers)
    #         print('hybridization false +ve', hybridized_results_false_positive)
    #         print('hybridization false -ve', hybridized_results_false_negative)
    #         print('false +ve', false_positive)
    #         print('false -ve', false_negative)
    #         print("\n")
    #         print("\n")

            # self.assertEqual(set(), false_positive, "false positive")
            # self.assertEqual(set(), false_negative, "false negative")

    def test_cgf_pred_one_primer(self):

        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"

        file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")

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
            count = -1
            for word in line.split():
                if word == '1':
                    genes_found_list.append(gene_list[count])
                count += 1
            file_gene_dict[line.split(None, 1)[0]] = genes_found_list

        file.close()

        files = [file for file in os.listdir(db_directory) if file.endswith('.fasta')]
        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            print('FILE!!!', file_path)
            file_name = file_path.partition(db_directory + "/")[2]
            print(file_name)
            f_out_file_path = db_directory + "/out_files/" + "f_" + file_name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + file_name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + file_name.replace("fasta", "xml")

            cgf_predictions = CGFPrediction.cgf_prediction(forward_primers, reverse_primers, file_path, f_out_file_path,
                                                           r_out_file_path,
                                                           amplicon_sequences, full_out_file_path)
            lo_predictions_names = []
            result = cgf_predictions[0]
            print('cgf predictions[0]', result)
            for quer in result:
                for hsp in quer:
                    lo_predictions_names.append(hsp.name)
                    # lo_predictions_names.append(quer[0].name)

            lo_diff_contig_genes = cgf_predictions[1]
            lo_single_primers = cgf_predictions[3]
            hybridized_results_same_contig = []
            for gene_name in cgf_predictions[2]:
                if "11168_" in gene_name:
                    tmp = gene_name[6:]
                    hybridized_results_same_contig.append(tmp)

            # blast_object = cgf_predictions[4]
            # forward_blast = cgf_predictions[5]
            # reverse_blast = cgf_predictions[6]
            # full_blast_qcov = cgf_predictions[7]

            genes_expected = file_gene_dict[file_name]
            genes_found = lo_predictions_names
            print('genes expected', genes_expected)
            print('genes found', genes_found)
            print('hybridization genes found', hybridized_results_same_contig)
            print('diff contig genes', lo_diff_contig_genes)
            print('single primers found', lo_single_primers)
            false_positive = set(genes_found) - set(genes_expected) - set(lo_diff_contig_genes)
            false_negative = set(genes_expected) - set(genes_found) - set(lo_diff_contig_genes)
            hybridized_results_false_positive = set(hybridized_results_same_contig) - set(genes_expected)
            hybridized_results_false_negative = set(genes_expected) - set(hybridized_results_same_contig)
            print('hybridization false +ve', hybridized_results_false_positive)
            print('hybridization false -ve', hybridized_results_false_negative)
            print('false +ve', false_positive)

            #testing
            # for hsp_name in false_positive:
            #     for hsp in blast_object.hsp_objects:
            #         if hsp_name in hsp.name:
            #             print('query', hsp.name, hsp.identities, hsp.query)
            #             print('sbjct', hsp.name, hsp.length, hsp.sbjct)

            print('false -ve', false_negative)
            lo_hsp_rejects = cgf_predictions[5]
            print('lo hsp rejects', lo_hsp_rejects)
            for hsp in lo_hsp_rejects:
                for hsp_name in false_negative:
                    if hsp.name == hsp_name:
                        print('possible false -ve', hsp.name, hsp.both_primers_found, hsp.contig, hsp.ehybrid, hsp.epcr, hsp.valid,
                              hsp.snp, hsp.pcr_distance)

            # for hsp_name in false_negative:
            #     for hsp in !!!:
            #         if hsp_name in hsp.name:
            #             print('primers found', hsp.both_primers_found)
            #             print('contig', hsp.contig)

            #testing
            # for hsp_name in false_negative:
            #     for hsp in full_blast_qcov.hsp_objects:
            #         if hsp_name in hsp.name:
            #             print('full', hsp.name, hsp.contig_name, hsp.identities / hsp.length, hsp.gaps)
            #             print(hsp.query)
            #             print(hsp.sbjct)
            #     for hsp in forward_blast.hsp_objects:
            #         if hsp_name in hsp.name:
            #             print('forward', hsp.name, hsp.contig_name, hsp.identities / hsp.length, hsp.gaps)
            #             print(hsp.query)
            #             print(hsp.sbjct)
            #     for hsp in reverse_blast.hsp_objects:
            #         if hsp_name in hsp.name:
            #             print('reverse', hsp.name, hsp.contig_name, hsp.identities / hsp.length, hsp.gaps)
            #             print(hsp.query)
            #             print(hsp.sbjct)

#
# # #TODO: you are currently testing the manipulated data you made for 11168 where you split gene cj0008 into contigs!
#     def test_11168_contigs_cj0008_manipulated(self):
#
#         forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
#         reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
#         forward_dict = CGFPrediction.create_primer_dict(forward_primers)
#         reverse_dict = CGFPrediction.create_primer_dict(reverse_primers)
#         amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
#         db_directory = "/home/sfisher/Sequences/11168_complete_genome/modified_genome_files"
#         result_dict = CGFPrediction.main(db_directory, forward_primers, reverse_primers, amplicon_sequences)
#         # self.assertEqual(len(result_dict), 1)
#
#         print("test 11168 contigs cj0008 manipulated")
#
#         files = [file for file in os.listdir(db_directory) if file.endswith('.fasta')]
#         for file in files:
#             print("file", file)
#             if result_dict[file] == None:
#                 self.assertIn("fail", file)
#             if result_dict[file] != None:
#                 result = result_dict[file]
#                 for lo_hsp in result:
#                     self.assertEqual(len(lo_hsp), 2)
#                     self.assertEqual(lo_hsp[0].name, lo_hsp[1].name)
#                     self.assertEqual(CGFPrediction.pcr_directly(lo_hsp[0], lo_hsp[1], amplicon_sequences, forward_dict, reverse_dict), True)
#                 self.assertEqual(len(result), 40)

if __name__ == '__main__':
    unittest.main()
