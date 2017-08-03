import unittest
import CGFPrediction
import errno
import os

MIN_BSR = 0.6

class MyTestCase(unittest.TestCase):
    def test_cgf_prediction_trial(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
        # db_directory = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
        db_directory = "/home/sfisher/Sequences/11168_test_files/memory_trial_CI-5768"
        file = open("/home/sfisher/Sequences/11168_test_files/CI-5768_results.txt", "r")
        # file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt", "r")
        # file = open("/home/sfisher/Sequences/11168_test_files/cgf_246_results_modified.txt", "r")

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
            cgf_predictions = CGFPrediction.cgf_prediction(forward_primers, reverse_primers, file_path,
                                                           amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)
            cgf_predictions_dict = cgf_predictions[0]
            # print('dict', cgf_predictions_dict)
            all_hsp_list = cgf_predictions[1]

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

            #TODO: possibly change all_hsp_list to a generator!!
            for gene_name in false_negative:
                # for hsp in all_hsp_dict[gene_name]:
                for hsp in all_hsp_list:
                    if hsp.name in gene_name:
                        result = CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0)
                        print('bsr', hsp.bsr)
                        if hsp.bsr < MIN_BSR:
                            print(hsp.name, "on contig", hsp.contig_name, "failed because hsp bsr is", hsp.bsr)
                        elif result == 2:
                            print(hsp.name, "CASE 2")
                            print('both primers found:', hsp.both_primers_found)
                            print('ehybrid is:', hsp.ehybrid)
                            # print('query', hsp.query)
                            print('sbjct', hsp.sbjct)
                            # print('amp query', hsp.amp_query)
                            # print('amp sbjct', hsp.amp_sbjct)
                        elif result == 3:
                            print(hsp.name, "CASE 3")
                            print('both primers found:', hsp.both_primers_found, 'with contig', hsp.contig)
                            print('ehybrid is:', hsp.ehybrid)
                        elif result == 4:
                            print(hsp.name, "on contig", hsp.contig_name, "CASE 4")
                            print('both primers found:', hsp.both_primers_found, 'with contig', hsp.contig)
                            print('ehybrid is:', hsp.ehybrid)
                            # print('query', hsp.query)
                            print('sbjct', hsp.sbjct)
                            print('partner & partner sbjct', hsp.partner.name, hsp.partner.contig_name, hsp.partner.sbjct)
                        elif result == 6 or result == 8 or result == 10:
                            print(hsp.name, "on contig", hsp.contig_name, "CASE 6 or 8 or 10")
                            print('both primers found:', hsp.both_primers_found)
                            print('ehybrid is:', hsp.ehybrid, 'with len', hsp.length)
                            # print('query', hsp.query)
                            print('sbjct', hsp.sbjct)
                            # print('amp query', hsp.amp_query)
                            # print('amp sbjct', hsp.amp_sbjct)
                        elif result == 11 or result == 15 or result == 137:
                            print('FALSE')
                        elif result == 12:
                            print(hsp.name, "CASE 12")
                            print('both primers found:', hsp.both_primers_found, 'with ehybrid', hsp.ehybrid, 'and location', hsp.location)
                            if hsp.strand:
                                print('leading strand', 'located', hsp.end_dist, 'from end')
                            elif not hsp.strand:
                                print('lagging strand', 'located,', hsp.end_dist, 'from end')
                            # print('query', hsp.query)
                            print('sbjct', hsp.sbjct)
                            # print('amp query', hsp.amp_query)
                            # print('amp sbjct', hsp.amp_sbjct)
                            print('Most likely too restrictive to find both primers...')
                        elif result == 20:
                            print(hsp.name, "on contig", hsp.contig_name, "CASE 20")
                            print('both primers found:', hsp.both_primers_found, 'with contig:', hsp.contig, 'and location', hsp.location)
                            if hsp.strand:
                                print('leading strand', 'located', hsp.end_dist, 'from end')
                                # print('primer query:', hsp.query)
                            elif not hsp.strand:
                                print('lagging strand', 'located,', hsp.end_dist, 'from end')
                            # print('amp query', hsp.amp_query)
                            # print('amp sbjct', hsp.amp_sbjct)
                            # print('query    ', hsp.query)
                        elif result == 34:
                            print(hsp.name, "CASE 34")
                            print('both primers found:', hsp.both_primers_found, 'with contig:', hsp.contig)
                            print('invalid direction to other hsp')
                        elif result == 67:
                            print(hsp.name, "CASE 67")
                            print('false b/c snp:', hsp.snp)
                            print('sbjct', hsp.sbjct)
                            # print('query', hsp.query)
                        elif result == 138:
                            print(hsp.name, "CASE 138")
                            print('false b/c dist btwn primers')
                        else:
                            assert False
                        print('\n')

            for gene_name in false_positive:
                for hsp in cgf_predictions_dict[gene_name]:
                    print(hsp.name)
                    print('strand', hsp.strand)
                    print('sbjct    ', hsp.sbjct)
                    print('\n')





                    # print('name', hsp.name)
                    # print('contig', hsp.contig_name)
                    # print('strand', hsp.strand)
                    # print('both primers found', hsp.both_primers_found)
                    # print('contig', hsp.contig)
                    # print('ehybrid', hsp.ehybrid)
                    # print('epcr', hsp.epcr)
                    # print('snp', hsp.snp)
                    # print('pcr distance', hsp.pcr_distance)
                    # print('location', hsp.location)
                    # print('valid', hsp.valid)
                    # print('result!!!', CGFPrediction.false_neg_pred(binary_tree_attrib, hsp, 0))
                    # print("\n")




            # print(binary_tree_attrib)
            # for i in range(11, 16):
            #     binary_tree_attrib.insert(i, str(i))
            # binary_tree_attrib.insert(16, 'valid')
            # for i in range(17, 23):
            #     binary_tree_attrib.insert(i, str(i))
            # for i in range(23,33):
            #     binary_tree_attrib.insert(i, None)
            # binary_tree_attrib.insert(33, 'snp')
            # binary_tree_attrib.insert(34, '34')
            # for i in range(35, 65):
            #     binary_tree_attrib.insert(i, None)
            # for i in range(65, 68):
            #     binary_tree_attrib.insert(i, str(i))
            # print(binary_tree_attrib)

            # for gene_name in false_positive:
            #     print(cgf_predictions_dict[gene_name])
            #
            # for gene_name in false_negative:
            #     for hsp in all_hsp_dict[gene_name]:
            #         if hsp.both_primers_found: #case 1
            #             if hsp.contig: #case 2
            #                 if hsp.ehybrid: #3
            #                     if hsp.epcr: #4
            #                         print('ERROR CASE 15: FALSE RESULT: EPCR SHOULD BE TRUE')
            #                     elif not hsp.epcr: #4
            #                         if hsp.valid: #5
            #                             if hsp.snp: #6
            #                                 print(hsp.name, '-ve b/c snp: CASE 65')
            #                                 if hsp.strand:
            #                                     print('leading strand', hsp.sbjct)
            #                                     print('query         ', hsp.query)
            #                                 elif not hsp.strand:
            #                                     print('lagging strand', hsp.sbjct)
            #                                     print('query         ', hsp.query)
            #                             elif not hsp.snp:
            #                                 if hsp.pcr_distance:
            #                                     print('ERROR: CASE 66: EPCR SHOULD HAVE PRODUCED TRUE')
            #                                 elif not hsp.pcr_distance:
            #                                     print(hsp.name, '-ve b/c distance btwn pcr primers: CASE 67')
            #                         elif not hsp.valid:
            #                             print(hsp.name, "-be b/c primers not facing each other: CASE 34")
            #                 elif not hsp.ehybrid:
            #                     print(hsp.name, '-ve b/c ehybrid with len found:', hsp.amp_len, 'CASE 17 (or possibly 8 or 18)')
            #             elif not hsp.contig:
            #                     if hsp.ehybrid:
            #                         if hsp.location:
            #                             print('error CASE 19!!!')
            #                         elif not hsp.location:
            #                             print(hsp.name, '-ve because not facing end or too far from end of contig: CASE 20')
            #                             if hsp.strand:
            #                                 print('leading strand', 'located', hsp.end_dist, 'from end')
            #                             elif not hsp.strand:
            #                                 print('lagging strand', 'located,', hsp.end_dist, 'from end')
            #                     elif not hsp.ehybrid:
            #                         print(hsp.name, '-ve b/c length failed: primers found on diff contigs with len of', hsp.contig_name, '=', hsp.amp_len, 'CASE 21 or 10 or 22')
            #                         print('amp sbjct', hsp.amp_sbjct)
            #                         print('amp query', hsp.amp_query)
            #         elif not hsp.both_primers_found:
            #             if hsp.ehybrid:
            #                 if hsp.location:
            #                     print('error: false result: location true: CASE 21')
            #                     print('amp sbjct', hsp.amp_sbjct)
            #                     print('amp query', hsp.amp_query)
            #                 elif not hsp.location:
            #                     print(hsp.name, '-ve b/c location on contig failed: CASE 12')
            #                     if hsp.strand:
            #                         print('leading strand', 'located', hsp.end_dist, 'from end')
            #                     elif not hsp.strand:
            #                         print('lagging strand', 'located,', hsp.end_dist, 'from end')
            #             elif not hsp.ehybrid:
            #                 print(hsp.name, '-ve b/c length (ehybrid): CASE 12 or 14', hsp.amp_len)
            #     print('\n')
            #     print('\n')







                    # print('name', hsp.name)
                    # print('both primers found', hsp.both_primers_found)
                    # print('contig', hsp.contig)
                    # print('ehybrid', hsp.ehybrid)
                    # print('epcr', hsp.epcr)
                    # print('snp', hsp.snp)
                    # print('pcr distance', hsp.pcr_distance)
                    # print('location', hsp.location)
                    # print('valid', hsp.valid)
                    # print('strand', hsp.strand)
                    # print("\n")





if __name__ == '__main__':
    unittest.main()
















