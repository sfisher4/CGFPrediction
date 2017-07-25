import unittest
import CGFPrediction
import errno
import os

class MyTestCase(unittest.TestCase):
    def test_cgf_prediction_trial(self):
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
            # print('FILE!!!', file_path)
            file_name = file_path.partition(db_directory + "/")[2]
            # print(file_name)
            f_out_file_path = db_directory + "/out_files/" + "f_" + file_name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + file_name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + file_name.replace("fasta", "xml")

            cgf_predictions_dict = CGFPrediction.cgf_prediction_trial(forward_primers, reverse_primers, file_path, f_out_file_path,
                                                           r_out_file_path,
                                                           amplicon_sequences, full_out_file_path)[0]
            # print('dict', cgf_predictions_dict)
            all_hsp_dict = CGFPrediction.cgf_prediction_trial(forward_primers, reverse_primers, file_path, f_out_file_path,
                                                           r_out_file_path,
                                                           amplicon_sequences, full_out_file_path)[1]

            genes_expected = file_gene_dict[file_name]
            genes_found = [key for key in cgf_predictions_dict]

            #testing
            # for gene_name in genes_found:
            #     for hsp in cgf_predictions_dict[gene_name]:
            #         print('name', hsp.name)
            #         print('both primers found', hsp.both_primers_found)
            #         print('contig', hsp.contig)
            #         print('ehybrid', hsp.ehybrid)
            #         print('epcr', hsp.epcr)
            #         print('pcr distance', hsp.pcr_distance)
            #         print('location', hsp.location)
            #         print('valid', hsp.valid)
            #         print('strand', hsp.strand)
            #         print("\n")
            # print('genes expected', genes_expected)
            # print('genes found', genes_found)

            false_positive = set(genes_found) - set(genes_expected)
            false_negative = set(genes_expected) - set(genes_found)

            if len(false_positive) > 0 or len(false_negative) > 0:
                print(file_name)
            if len(false_positive) > 0:
                print('false positive', false_positive)
            if len(false_negative) > 0:
                print('false negative', false_negative)


            for gene_name in false_positive:
                print(cgf_predictions_dict[gene_name])
            for gene_name in false_negative:
                for hsp in all_hsp_dict[gene_name]:
                    if hsp.both_primers_found: #case 1
                        if hsp.contig: #case 2
                            if hsp.ehybrid: #3
                                if hsp.epcr: #4
                                    print('ERROR CASE 26: FALSE RESULT: EPCR SHOULD BE TRUE')
                                elif not hsp.epcr: #4
                                    if hsp.valid: #5
                                        if hsp.snp: #6
                                            print(hsp.name, '-ve b/c snp: CASE 7') #TODO: determine where the snp was
                                        elif not hsp.snp:
                                            if hsp.pcr_distance:
                                                print('ERROR: CASE 24: EPCR SHOULD HAVE PRODUCED TRUE')
                                            elif not hsp.pcr_distance:
                                                print(hsp.name, '-ve b/c distance btwn pcr primers: CASE 25')
                                    elif not hsp.valid:
                                        print(hsp.name, "-be b/c primers not facing each other: CASE 9")
                            elif not hsp.ehybrid:
                                print(hsp.name, '-ve b/c ehybrid with len found:', hsp.amp_len, 'CASE 10')
                        elif not hsp.contig:
                            if hsp.ehybrid:
                                if hsp.location:
                                    print('error CASE 13!!!')
                                elif not hsp.location:
                                    print(hsp.name, '-ve because not facing end or too far from end of contig: CASE 16') #TODO: determine what way it was facing and how far it is from the end
                            elif not hsp.ehybrid:
                                print(hsp.name, '-ve b/c length failed: primers found on diff contigs with len of', hsp.contig_name, '=', hsp.amp_len, 'CASE 17')
                    elif not hsp.both_primers_found:
                        if hsp.ehybrid:
                            if hsp.location:
                                print('error: false result: location true: CASE 21')
                            elif not hsp.location:
                                print(hsp.name, '-ve b/c location on contig failed: CASE 22') #TODO: determine what way it was facing and how far it is from the end
                        elif not hsp.ehybrid:
                            print(hsp.name, '-ve b/c length (ehybrid): CASE 23', hsp.amp_len)
                print('\n')
                print('\n')







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
















