import unittest
from HSP import HSP
import os
import errno

import PCRPrediction

HSP_THRESHOLD = 0.9
E_VALUE_THRESHOLD = 0.04

    #Files available that will probably not be used
    # database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta"  # contains id and a complete genome sequence
    # full_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"
    # amp_vs_amp_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"
    # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)
    # cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta"  # complete


class TestWithAmp(unittest.TestCase):

    # @classmethod
    # def setUpClass(cls):
        #tODO?

    def test_create_blastn_object(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_perfect_hit"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)

            self.assertEqual(len(forward_blast_object.blast_records), 40)
            self.assertEqual(len(reverse_blast_object.blast_records), 40)
            self.assertEqual(len(full_blast_object.blast_records), 40)

    def test_create_hsp_objects(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_perfect_hit"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)

            print(file_path)

            self.assertEqual(len(forward_blast_object.hsp_objects), 1)
            self.assertEqual(len(reverse_blast_object.hsp_objects), 1)
            if "contig" in file_path:
                self.assertEqual(len(full_blast_object.hsp_objects), 2)
            else:
                self.assertEqual(len(full_blast_object.hsp_objects), 1)

            for hsp_object in forward_blast_object.hsp_objects:
                self.assertNotEqual(hsp_object.contig_name, "", 'hsp_object contig name is not initialized')
                self.assertGreater(E_VALUE_THRESHOLD, hsp_object.expect)
                self.assertNotEqual(hsp_object.expect, -1)
                self.assertNotEqual(hsp_object.start, -1)
                self.assertNotEqual(hsp_object.end, -1)
                self.assertNotEqual(hsp_object.start, hsp_object.end)
                if hsp_object.start < hsp_object.end:
                    self.assertTrue(hsp_object.strand)
                if hsp_object.end < hsp_object.start:
                    self.assertFalse(hsp_object.strand)
                self.assertIsNotNone(hsp_object.strand)
                self.assertNotEqual(hsp_object.length, 0, "length of hsp is 0")
                self.assertEqual(hsp_object.length,
                                 abs(hsp_object.end - hsp_object.start) + 1)
                self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)
                self.assertLess(hsp_object.expect, E_VALUE_THRESHOLD)

                # tests db_length
                if "contig" in file_path:
                    if 'NODE_1' in hsp_object.contig_name:
                        self.assertEqual(hsp_object.db_length, 376)
                    elif 'NODE_2' in hsp_object.contig_name:
                        self.assertEqual(hsp_object.db_length, 236)
                        test = True  # make sure both contigs are reached
                    else:
                        self.assertWarns("shouldn't reach here!!!")

    def test_pcr_directly(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_perfect_hit"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)

            f_primer_dict = PCRPrediction.create_primer_dict(forward_primers)
            r_primer_dict = PCRPrediction.create_primer_dict(reverse_primers)

            if "contig" not in file_path:
                if "cj0008" in file_path:
                    # f_hsp_cj0008 = [hsp for hsp in f_hsp_object if "cj0008" in file_path and hsp.start == 1 and hsp.end == 20]
                    # r_hsp_cj0008 = [hsp for hsp in r_hsp_object if "cj0008" in file_path and hsp.start == 1 and hsp.end == 22]
                    # self.assertEqual(len(f_hsp_cj0008), 1)
                    # self.assertEqual(len(r_hsp_cj0008), 1)
                    for f_hsp_object in forward_blast_object.hsp_objects:
                        for r_hsp_object in reverse_blast_object.hsp_objects:
                            result = PCRPrediction.pcr_directly(f_hsp_object, r_hsp_object, amplicon_sequences, f_primer_dict, r_primer_dict)
                            if f_hsp_object.length == 20 and r_hsp_object.length == 22:
                                self.assertTrue(result)
                            else:
                                self.assertFalse(result)
                elif "cj0483_complete" in file_path:
                    for f_hsp_object in forward_blast_object.hsp_objects:
                        for r_hsp_object in reverse_blast_object.hsp_objects:
                            result = PCRPrediction.pcr_directly(f_hsp_object, r_hsp_object, amplicon_sequences, f_primer_dict, r_primer_dict)
                            if f_hsp_object.length == 20 and r_hsp_object.length == 21:
                                self.assertTrue(result)
                            else:
                                self.assertFalse(result)
                else:
                    self.assertTrue(False) #should never reach here

    def test_entire_gene(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_perfect_hit"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]
        f_primer_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primer_dict = PCRPrediction.create_primer_dict(reverse_primers)

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)

            for hsp in full_blast_object.hsp_objects:
                if "NODE" in hsp.contig_name: #Assuming not the case where both primers are located on the same contig
                    results = PCRPrediction.entire_gene(full_blast_object, hsp, f_primer_dict, r_primer_dict)
                    self.assertEqual(len(results), 2)
                    self.assertIn(hsp, results)



class TestGeneAnnotationError(unittest.TestCase):

    # @classmethod
    # def setUpClass(cls):
    #     cls._amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     cls._forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
    #     cls._reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
    #     db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
    #     cls._f_primers_dict = PCRPrediction.create_primer_dict(cls._forward_primers)
    #     cls._r_primers_dict = PCRPrediction.create_primer_dict(cls._reverse_primers)
    #
    #     files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]
    #
    #     files_paths = []
    #     for file in files:
    #         files_paths.append(os.path.abspath(db_directory) + '/' + file)
    #
    #     # create new folder for out files
    #     try:
    #         os.mkdir(db_directory + "/out_files/")
    #     except OSError as exc:
    #         if exc.errno != errno.EEXIST:
    #             raise
    #         pass
    #
    #     cls._file_paths_dict = {}
    #     for file_path in files_paths:
    #         name = file_path.partition(db_directory + "/")[2]
    #         f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
    #         r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
    #         full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
    #         forward_blast_object = PCRPrediction.create_blastn_object(cls._forward_primers, file_path, f_out_file_path)
    #         reverse_blast_object = PCRPrediction.create_blastn_object(cls._reverse_primers, file_path, r_out_file_path)
    #         full_blast_object = PCRPrediction.create_blastn_object(cls._amplicon_sequences, file_path, full_out_file_path)
    #         objects = []
    #         objects.append(forward_blast_object)
    #         cls._file_paths_dict[file_path] = objects

    #TODO: make tests with different names!!!
    #TODO: test as though one is found but the other may be found???
    #TODO: could both be found?!!!
    def test_create_blast_object(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"

        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        file_paths_dict = {}
        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)
            objects = []
            objects.append(forward_blast_object)
            file_paths_dict[file_path] = objects

            self.assertGreaterEqual(len(forward_blast_object.blast_records), 40)
            self.assertGreaterEqual(len(reverse_blast_object.blast_records), 40)
            self.assertGreaterEqual(len(full_blast_object.blast_records), 40)
            self.assertGreaterEqual(len(forward_blast_object.hsp_objects), 1)
            self.assertGreaterEqual(len(reverse_blast_object.hsp_objects), 1)
            self.assertGreaterEqual(len(full_blast_object.hsp_objects), 1)

    def test_pcr_directly(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
        f_primers_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primers_dict = PCRPrediction.create_primer_dict(reverse_primers)

        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        count = 0
        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)

            for f_hsp in forward_blast_object.hsp_objects:
                r_hsps = [hsp for hsp in reverse_blast_object.hsp_objects if hsp.contig_name == f_hsp.contig_name and hsp.name == f_hsp.name]
                for r_hsp in r_hsps:
                    count += 1
                    result = PCRPrediction.pcr_directly(f_hsp, r_hsp, amplicon_sequences, f_primers_dict, r_primers_dict)
                    self.assertEqual(result, True)
        self.assertEqual(count, 2)

    def test_is_snp(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
        f_primers_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primers_dict = PCRPrediction.create_primer_dict(reverse_primers)

        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)

            for hsp in full_blast_object.hsp_objects:
                if "NODE" in hsp.contig_name:
                    results = PCRPrediction.entire_gene(full_blast_object, hsp, f_primers_dict, r_primers_dict)
                    for result in results:
                        self.assertFalse(result.snp)


    def test_entire_gene(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
        f_primers_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primers_dict = PCRPrediction.create_primer_dict(reverse_primers)

        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        count = 0
        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)

            for hsp in full_blast_object.hsp_objects:
                if "NODE" in hsp.contig_name:
                    count += 1
                    lo_result = PCRPrediction.entire_gene(full_blast_object, hsp, f_primers_dict, r_primers_dict)
                    self.assertEqual(len(lo_result), 2)
                    self.assertIn(hsp, lo_result)
        self.assertEqual(count, 4)

    #TODO: more thorough tests here?
    def test_pcr_prediction(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        dict_predictions = PCRPrediction.main("/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error", forward_primers, reverse_primers, amplicon_sequences)

        for pred in dict_predictions:
            self.assertEqual(len(dict_predictions[pred]), 1)
            for matches in dict_predictions[pred]:
                self.assertEqual(len(matches), 2)


class TestSNPPrimers(unittest.TestCase):

    def test_snp(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_snp"
        f_primer_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primer_dict = PCRPrediction.create_primer_dict(reverse_primers)
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

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
            print(file_path)
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)

            self.assertGreaterEqual(len(forward_blast_object.hsp_objects), 1)
            self.assertGreaterEqual(len(reverse_blast_object.hsp_objects), 1)
            if "6" in file_path:
                for hsp in forward_blast_object.hsp_objects:
                    self.assertIsNone(hsp.snp)
                    PCRPrediction.is_snp_primer_search(hsp, f_primer_dict, r_primer_dict)
                    self.assertFalse(hsp.snp)
                for hsp in reverse_blast_object.hsp_objects:
                    self.assertIsNone(hsp.snp)
                    PCRPrediction.is_snp_primer_search(hsp, f_primer_dict, r_primer_dict)
                    self.assertFalse(hsp.snp)
            else:
                for hsp in forward_blast_object.hsp_objects:
                    self.assertIsNone(hsp.snp)
                    PCRPrediction.is_snp_primer_search(hsp, f_primer_dict, r_primer_dict)
                    self.assertTrue(hsp.snp)
                if "fw" in file_path:
                    for hsp in reverse_blast_object.hsp_objects:
                        self.assertIsNone(hsp.snp)
                        PCRPrediction.is_snp_primer_search(hsp, f_primer_dict, r_primer_dict)
                        self.assertFalse(hsp.snp)
                else:
                    for hsp in reverse_blast_object.hsp_objects:
                        self.assertIsNone(hsp.snp)
                        PCRPrediction.is_snp_primer_search(hsp, f_primer_dict, r_primer_dict)
                        self.assertTrue(hsp.snp)

    def test_pcr_prediction_snp(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        dict_predictions = PCRPrediction.main("/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_snp", forward_primers, reverse_primers, amplicon_sequences)


        predictions_ws = [pred for pred in dict_predictions if "_1" in pred or "_5" in pred]
        predictions = [pred for pred in dict_predictions if pred not in predictions_ws]
        for pred in predictions_ws:
            self.assertEqual([], dict_predictions[pred])
        for pred in predictions:
            self.assertNotEqual([], dict_predictions[pred])
            self.assertEqual(len(dict_predictions[pred]), 1)
            for matches in dict_predictions[pred]:
                self.assertEqual(len(matches), 2)

"""

class TestSNPEntireGene(unittest.TestCase):

    #TODO:!!!
    #SNP_THRESHOLD = 5
    def test_entire_gene(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_snp_contigs"
        f_primers_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primers_dict = PCRPrediction.create_primer_dict(reverse_primers)

        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

        files_paths = []
        for file in files:
            files_paths.append(os.path.abspath(db_directory) + '/' + file)

        # create new folder for out files
        try:
            os.mkdir(db_directory + "/out_files/")
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass

        for file_path in files_paths:
            print(file_path)
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)
            full_blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)
            self.assertGreaterEqual(len(full_blast_object.hsp_objects), 2)
            for hsp_object in full_blast_object.hsp_objects:
                # print(hsp_object.contig_name)
                results = PCRPrediction.entire_gene(full_blast_object, hsp_object, f_primers_dict, r_primers_dict)
                if "6" in file_path:
                    self.assertEqual(len(results), 2)
                    self.assertIn(hsp_object, results)
                    self.assertFalse(hsp_object.snp)
                elif "one_contig" in file_path:
                    self.assertEqual(len(results), 1)
                    for result in results:
                        self.assertIn("NODE_2", result.contig_name)
                        self.assertNotIn("NODE_1", result.contig_name)
                        self.assertFalse(result.snp)
                    other_contig = [hsp for hsp in full_blast_object.hsp_objects if hsp not in results]
                    for hsp in other_contig:
                        self.assertTrue(hsp.snp)
                    #fw snp is true
                elif "10+" in file_path:
                    self.assertEqual(len(results), 2)
                    self.assertIn(hsp_object, results)
                    self.assertFalse(hsp_object.snp)
                else:
                    self.assertEqual(len(results), 0)
                    self.assertTrue(hsp_object.snp)
"""




class TestContigTrunc(unittest.TestCase):

    # @classmethod
    # def setUpClass(cls):
    #     #query
    #     cls._amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    #     cls._f_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
    #     cls._r_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
    #
    #     #databases
    #     both_contigs_306_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_306_bp.fasta"
    #     both_contigs_27_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_27_bp.fasta"
    #     both_contigs_60_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_60_bp.fasta"
    #     both_contigs_61_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_61_bp.fasta"
    #     cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
    #     cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
    #     cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
    #     both_contigs_61bp_db_name_mystery = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db_entire_gene_61_bp_each_contig.fasta"
    #     #one primer databases
    #     one_primer_lead_end = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/one_primer_lead_end.fasta"
    #     one_primer_lag_end = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/one_primer_lag_end.fasta"
    #     one_primer_lead_start = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/one_primer_lead_start.fasta"
    #     one_primer_lag_start = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/one_primer_lag_start.fasta"
    #     one_primer_found_all = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/one_primer_found_all.fasta"
    #
    #     #out files
    #     out_file_306 = '/home/sfisher/Sequences/blast_record/test_blast_306.xml'
    #     out_file_27 = '/home/sfisher/Sequences/blast_record/test_blast_27.xml'
    #     out_file_60 = '/home/sfisher/Sequences/blast_record/test_blast_60.xml'
    #     out_file_61 = '/home/sfisher/Sequences/blast_record/test_blast_61.xml'
    #     out_file_contig_full = '/home/sfisher/Sequences/blast_record/test_blast_contig_full.xml'
    #     out_file_contig_trunc = '/home/sfisher/Sequences/blast_record/test_blast_contig_trunc.xml'
    #     out_file_complete = '/home/sfisher/Sequences/blast_record/test_blast_complete.xml'
    #     out_file_mystery = '/home/sfisher/Sequences/blast_record/test_blast_mystery_61bp.xml'
    #     out_file_lead_end = '/home/sfisher/Sequences/blast_record/out_file_lead_end.xml'
    #     out_file_lag_end = '/home/sfisher/Sequences/blast_record/out_file_lag_end.xml'
    #     out_file_lead_start = '/home/sfisher/Sequences/blast_record/out_file_lead_start.xml'
    #     out_file_lag_start = '/home/sfisher/Sequences/blast_record/out_file_lag_start.xml'
    #     out_file_found_all = '/home/sfisher/Sequences/blast_record/out_file_found_all.xml'
    #
    #     #contig full
    #     cls._blast_object_contig_full = PCRPrediction.create_blastn_object(cls._amplicon_sequences, cj0483_contig_full_amp_seq, out_file_contig_full)
    #     # cls._lo_hsp_objects_contig_full = PCRPrediction.create_hsp_objects(cls._blast_object_contig_full)
    #
    #     #contig truncated
    #     cls._blast_object_contig_trunc = PCRPrediction.create_blastn_object(cls._amplicon_sequences, cj0483_contig_trunc_amp_seq, out_file_contig_trunc)
    #     # cls._lo_hsp_objects_contig_trunc = PCRPrediction.create_hsp_objects(cls._blast_object_contig_trunc)
    #
    #     #complete amp (no contigs)
    #     cls._blast_object_complete = PCRPrediction.create_blastn_object(cls._amplicon_sequences, cj0483_complete_amp_seq, out_file_complete)
    #     # cls._lo_hsp_objects_complete = PCRPrediction.create_hsp_objects(cls._blast_object_complete)
    #
    #     #both contigs length 306 bp
    #     cls._blast_obj_306= PCRPrediction.create_blastn_object(cls._amplicon_sequences, both_contigs_306_bp, out_file_306)
    #     # cls._lo_hsp_objects_306 = PCRPrediction.create_hsp_objects(cls._blast_obj_306)
    #
    #     #both contigs length 27 bp
    #     cls._blast_obj_27 = PCRPrediction.create_blastn_object(cls._amplicon_sequences, both_contigs_27_bp, out_file_27)
    #     # cls._lo_hsp_objects_27 = PCRPrediction.create_hsp_objects(cls._blast_obj_27)
    #
    #     #both contigs length 60 bp
    #     cls._blast_obj_60 = PCRPrediction.create_blastn_object(cls._amplicon_sequences, both_contigs_60_bp, out_file_60)
    #     # cls._lo_hsp_objects_60 = PCRPrediction.create_hsp_objects(cls._blast_obj_60)
    #
    #     #both contigs length 61 bp
    #     cls._blast_obj_61 = PCRPrediction.create_blastn_object(cls._amplicon_sequences, both_contigs_61_bp, out_file_61)
    #     # cls._lo_hsp_objects_61 = PCRPrediction.create_hsp_objects(cls._blast_obj_61)
    #
    #     #both contigs length 61 bp with a different db name then amplicon sequences
    #     cls._blast_obj_61_mystery = PCRPrediction.create_blastn_object(cls._amplicon_sequences, both_contigs_61bp_db_name_mystery, out_file_mystery)
    #     # cls._lo_hsp_objects_61_mystery = PCRPrediction.create_hsp_objects(cls._blast_obj_61_mystery)
    #
    #     #one primer found on leading strand at end of seq
    #     cls._blast_obj_lead_end = PCRPrediction.create_blastn_object(cls._amplicon_sequences, one_primer_lead_end, out_file_lead_end)
    #     # cls._lo_hsp_objects_lead_end = PCRPrediction.create_hsp_objects(cls._blast_obj_lead_end)
    #     #one primer found on lagging strand at end of seq
    #     cls._blast_obj_lag_end = PCRPrediction.create_blastn_object(cls._amplicon_sequences, one_primer_lag_end, out_file_lag_end)
    #     # cls._lo_hsp_objects_lag_end = PCRPrediction.create_hsp_objects(cls._blast_obj_lag_end)
    #     #one primer found on leading strand at start of seq
    #     cls._blast_obj_lead_start = PCRPrediction.create_blastn_object(cls._amplicon_sequences, one_primer_lead_start, out_file_lead_start)
    #     # cls._lo_hsp_objects_lead_start = PCRPrediction.create_hsp_objects(cls._blast_obj_lead_start)
    #     #one primer found on lagging strand at start of seq
    #     cls._blast_obj_lag_start = PCRPrediction.create_blastn_object(cls._amplicon_sequences, one_primer_lag_start, out_file_lag_start)
    #     # cls._lo_hsp_objects_lag_start = PCRPrediction.create_hsp_objects(cls._blast_obj_lag_start)
    #     #one primer found that covers the whole seq
    #     cls._blast_obj_found_all = PCRPrediction.create_blastn_object(cls._amplicon_sequences, one_primer_found_all, out_file_found_all)
    #     # cls._lo_hsp_object_found_all = PCRPrediction.create_hsp_objects(cls._blast_obj_found_all)

    #contigs trunc
    #search full gene
    # def test_create_blast_records_full_search_contigs_truncated(self):
    #     self.assertNotEqual([], self._blast_object_contig_trunc.blast_records)
    #     self.assertEqual(len(self._blast_object_contig_trunc.blast_records), 40)

    def test_create_blast_objects_full_search(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_contig_trunc"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

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
            print(file_path)
            name = file_path.partition(db_directory + "/")[2]
            full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
            blast_object = PCRPrediction.create_blastn_object(amplicon_sequences, file_path, full_out_file_path)
            self.assertGreaterEqual(len(blast_object.blast_records), 40)
            self.assertGreaterEqual(len(blast_object.hsp_objects), 2)

            for hsp_object in blast_object.hsp_objects:
                self.assertNotEqual(hsp_object.contig_name, "", 'hsp_object contig name is not initialized')
                self.assertGreater(E_VALUE_THRESHOLD, hsp_object.expect)
                self.assertNotEqual(hsp_object.expect, -1)
                self.assertNotEqual(hsp_object.start, -1)
                self.assertNotEqual(hsp_object.end, -1)
                self.assertNotEqual(hsp_object.start, hsp_object.end)
                if hsp_object.start < hsp_object.end:
                    self.assertTrue(hsp_object.strand)
                if hsp_object.end < hsp_object.start:
                    self.assertFalse(hsp_object.strand)
                self.assertIsNotNone(hsp_object.strand)
                self.assertNotEqual(hsp_object.length, 0, "length of hsp is 0")
                self.assertEqual(hsp_object.length, abs(hsp_object.end - hsp_object.start) + 1)
                self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)
                self.assertLess(hsp_object.expect, E_VALUE_THRESHOLD)

    #TODO: you are here!!!
    def test_pcr_prediction(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        f_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        r_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_contig_trunc"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

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
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            full_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            # forward_blast_object = PCRPrediction.create_blastn_object(f_primers, file_path, f_out_file_path)
            # reverse_blast_object = PCRPrediction.create_blastn_object(r_primers, file_path, r_out_file_path)
            lo_hsp_predictions = PCRPrediction.pcr_prediction(f_primers, r_primers, db_directory, f_out_file_path, r_out_file_path, amplicon_sequences, full_out_file_path)
            # if "61" in file_path:
            #     self.assertEqual(len(lo_hsp_predictions), 2)
            #     for hsp in lo_hsp_predictions:



#     #complete amp
#     #search full gene
#     def test_create_hsp_objects_full_search_contigs_complete(self):
#
#         self.assertGreaterEqual(len(self._blast_object_complete.hsp_objects), 1)
#
#         for hsp in self._blast_object_complete.hsp_objects:
#             self.assertLessEqual(HSP_THRESHOLD, (hsp.identities / 612))
#
#     # # complete amp, random bp removed < % identities, pass when HSP_THRESHOLD >= 0.9
#     # # 306 bp on each contig,
#     # def test_create_hsp_records_low_perc_identities(self):
#     #
#     #
#     # #complete amp, random bp removed > % identities
#     # #306 on both strands, with
#     # def test_create_hsp_records_high_perc_identities(self):
#
# #TODO: combine create_hsp_objects tests
#     # contigs trunc
#     # full gene search
#     def test_create_hsp_objects_full_search_contigs(self):
#         self.assertGreaterEqual(len(self._blast_object_contig_trunc.hsp_objects), 2)
#
#         for hsp in self._blast_object_contig_trunc.hsp_objects:
#                 if "NODE_1" in hsp.contig_name:
#                     self.assertTrue(HSP_THRESHOLD <= hsp.identities / 376)
#                 elif "NODE_2" in hsp.contig_name:
#                     self.assertTrue(HSP_THRESHOLD <= hsp.identities / 228)
#                 else:
#                     self.assertFalse(True, 'should never reach here unless more than two nodes or error')
#
#         test = False
#         for hsp_object in self._blast_object_contig_trunc.hsp_objects:
#             # tests name
#             self.assertIn('cj0483', hsp_object.name)
#             # tests db_length
#             if 'NODE_1' in hsp_object.contig_name:
#                 self.assertEqual(hsp_object.db_length, 376)
#             elif 'NODE_2' in hsp_object.contig_name:
#                 self.assertEqual(hsp_object.db_length, 228)
#                 test = True  # make sure both contigs are reached
#             else:
#                 self.assertWarns("shouldn't reach here!!!")
#
#             self.assertNotEqual(hsp_object.contig_name, "", 'hsp_object contig name is not initialized')
#             self.assertGreater(E_VALUE_THRESHOLD, hsp_object.expect)
#             self.assertNotEqual(hsp_object.expect, -1)
#             self.assertNotEqual(hsp_object.start, -1)
#             self.assertNotEqual(hsp_object.end, -1)
#             self.assertNotEqual(hsp_object.start, hsp_object.end)
#             if hsp_object.start < hsp_object.end:
#                 self.assertTrue(hsp_object.strand)
#             if hsp_object.end < hsp_object.start:
#                 self.assertFalse(hsp_object.strand)
#             self.assertIsNotNone(hsp_object.strand)
#             self.assertNotEqual(hsp_object.length, 0, "length of hsp is 0")
#             self.assertEqual(hsp_object.length, abs(hsp_object.end - hsp_object.start) + 1)
#             self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)
#             self.assertLess(hsp_object.expect, E_VALUE_THRESHOLD)
#
#         self.assertTrue(test)  # makes sure both contigs are reached (db_length test)
#
#     #contigs trunc (61 each contig)
#     #full gene search
#     def test_create_hsp_objects_full_search_contigs_mystery(self):
#         self.assertGreaterEqual(len(self._blast_obj_61_mystery.hsp_objects), 2)
#
#         for hsp in self._blast_obj_61_mystery.hsp_objects:
#             self.assertTrue(HSP_THRESHOLD <= hsp.identities / 61)
#
#         for hsp_object in self._blast_obj_61_mystery.hsp_objects:
#             # tests name
#             self.assertIn('cj0483', hsp_object.name)
#             # tests db_length
#             self.assertEqual(hsp_object.db_length, 61)
#             self.assertNotEqual(hsp_object.contig_name, "", 'hsp_object contig name is not initialized')
#             self.assertGreater(E_VALUE_THRESHOLD, hsp_object.expect)
#             self.assertNotEqual(hsp_object.expect, -1)
#             self.assertNotEqual(hsp_object.start, -1)
#             self.assertNotEqual(hsp_object.end, -1)
#             self.assertNotEqual(hsp_object.start, hsp_object.end)
#             if hsp_object.start < hsp_object.end:
#                 self.assertTrue(hsp_object.strand)
#             if hsp_object.end < hsp_object.start:
#                 self.assertFalse(hsp_object.strand)
#             self.assertIsNotNone(hsp_object.strand)
#             self.assertNotEqual(hsp_object.length, 0, "length of hsp is 0")
#             self.assertEqual(hsp_object.length, abs(hsp_object.end - hsp_object.start) + 1)
#             self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)
#             self.assertLess(hsp_object.expect, E_VALUE_THRESHOLD)
#
#     #one_primer_lead_end
#     def test_create_hsp_objects_one_primer_lead_end(self):
#         self.assertGreaterEqual(len(self._blast_obj_lead_end.blast_records), 40)
#         self.assertGreaterEqual(len(self._blast_obj_lead_end.hsp_objects), 1)
#
#         for hsp in self._blast_obj_61_mystery.hsp_objects:
#             self.assertTrue(HSP_THRESHOLD <= hsp.identities / 61)
#
#         for hsp_object in self._blast_obj_61_mystery.hsp_objects:
#             # tests name
#             self.assertIn('cj0483', hsp_object.name)
#             # tests db_length
#             self.assertEqual(hsp_object.db_length, 61)
#             self.assertNotEqual(hsp_object.contig_name, "", 'hsp_object contig name is not initialized')
#             self.assertGreater(E_VALUE_THRESHOLD, hsp_object.expect)
#             self.assertNotEqual(hsp_object.expect, -1)
#             self.assertNotEqual(hsp_object.start, -1)
#             self.assertNotEqual(hsp_object.end, -1)
#             self.assertNotEqual(hsp_object.start, hsp_object.end)
#             if hsp_object.start < hsp_object.end:
#                 self.assertTrue(hsp_object.strand)
#             if hsp_object.end < hsp_object.start:
#                 self.assertFalse(hsp_object.strand)
#             self.assertIsNotNone(hsp_object.strand)
#             self.assertNotEqual(hsp_object.length, 0, "length of hsp is 0")
#             self.assertEqual(hsp_object.length, abs(hsp_object.end - hsp_object.start) + 1)
#             self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)
#             self.assertLess(hsp_object.expect, E_VALUE_THRESHOLD)

#TODO: look over!!!
    def test_pcr_prediction_contig_trunc(self):
        db = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_contig_trunc"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

        dict_predictions = PCRPrediction.main(db, forward_primers, reverse_primers, amplicon_sequences)

        predictions_pass = [pred for pred in dict_predictions if "61" in pred or "lag_start" in pred or "lead_end" in pred or "306" in pred or "found_all" in pred]
        predictions = [pred for pred in dict_predictions if pred not in predictions_pass]
        for pred in predictions_pass:
            print('pred', pred)
            self.assertNotEqual([], dict_predictions[pred], pred)
            self.assertEqual(len(dict_predictions[pred]), 1)
            for matches in dict_predictions[pred]:
                self.assertEqual(len(matches), 2)
        for pred in predictions:
            self.assertEqual([], dict_predictions[pred])

    def test_entire_gene(self):
        self.assertTrue(False)
        #TODO!!!
    #
    # def test_one_primer_lead_end(self):
    #     name = "cj0483"
    #     cj0483_object = HSP(name)
    #     lo_hsps = PCRPrediction.entire_gene(self._blast_obj_lead_end.hsp_objects, cj0483_object, self._f_primers, self._r_primers)
    #     self.assertGreaterEqual(len(self._blast_obj_lead_end.hsp_objects), 1)
    #     for hsp_object in self._blast_obj_lead_end.hsp_objects:
    #         self.assertEqual(hsp_object.strand, True)
    #         self.assertNotEqual(hsp_object.end, hsp_object.query_end)
    #         self.assertGreaterEqual(hsp_object.length, 61) #CUTOFF_GENE_LENGTH = 60
    #     self.assertEqual(len(lo_hsps), 1)
    #     self.assertIn('cj0483', lo_hsps[0].name)

    #
    #
    # #306 bp on each contig
    # #TODO: finish testing me
    # def test_entire_gene_306(self):
    #     for hsp_object in self._lo_hsp_objects_306:
    #         lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_306, hsp_object, self._f_primers, self._r_primers)
    #         self.assertTrue(len(lo_queries) == 2)
    #         for query in lo_queries:
    #             self.assertEqual(hsp_object.name, query.name)
    #             self.assertIsNone(query.valid)
    #
    # #27 bp on each contig
    # #TODO: finish testing me
    # def test_entire_gene_27(self):
    #     self.assertGreaterEqual(len(self._lo_hsp_objects_27), 2)
    #     for hsp_object in self._lo_hsp_objects_27:
    #         lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_27, hsp_object, self._f_primers, self._r_primers)
    #         self.assertTrue(len(lo_queries) == 0)
    #
    # def test_entire_gene_60(self):
    #     self.assertTrue(len(self._lo_hsp_objects_60) >= 2)
    #     for hsp_object in self._lo_hsp_objects_60:
    #         lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_60, hsp_object, self._f_primers, self._r_primers)
    #         self.assertEqual(len(lo_queries), 0)
    #
    # def test_entire_gene_61(self):
    #     self.assertTrue(len(self._lo_hsp_objects_61) >= 2)
    #     for hsp_object in self._lo_hsp_objects_61:
    #         lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_61, hsp_object, self._f_primers, self._r_primers)
    #         self.assertTrue(len(lo_queries) == 2)
    #         for query in lo_queries:
    #             self.assertEqual(hsp_object.name, query.name)
    #             self.assertIsNone(query.valid)
    #
    # def test_entire_gene_61_db_name_different(self):
    #     self.assertGreaterEqual(len(self._lo_hsp_objects_61_mystery), 2)
    #     for hsp_object in self._lo_hsp_objects_61_mystery:
    #         lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_61_mystery, hsp_object, self._f_primers, self._r_primers)
    #         self.assertTrue(len(lo_queries) == 2)
    #         for query in lo_queries:
    #             self.assertEqual(hsp_object.name, query.name)
    #             self.assertIsNone(query.valid)

    #
    # def test_one_primer_lag_end(self):
    #     name = "cj0483"
    #     cj0483_object = HSP(name)
    #     lo_hsps = PCRPrediction.entire_gene(self._lo_hsp_objects_lag_end, cj0483_object, self._f_primers, self._r_primers)
    #     self.assertEqual(len(lo_hsps), 0)
    #
    # def test_one_primer_lead_start(self):
    #     name = "cj0483"
    #     cj0483_object = HSP(name)
    #     lo_hsps = PCRPrediction.entire_gene(self._lo_hsp_objects_lead_start, cj0483_object, self._f_primers, self._r_primers)
    #     self.assertEqual(len(lo_hsps), 0)
    #
    # def test_one_primer_lag_start(self):
    #     name = "cj0483"
    #     cj0483_object = HSP(name)
    #     lo_hsps = PCRPrediction.entire_gene(self._lo_hsp_objects_lead_end, cj0483_object, self._f_primers, self._r_primers)
    #     self.assertGreaterEqual(len(self._lo_hsp_objects_lead_end), 1)
    #     for hsp_object in self._lo_hsp_objects_lead_end:
    #         self.assertEqual(hsp_object.strand, True)
    #         self.assertNotEqual(hsp_object.end, hsp_object.query_end)
    #         self.assertGreaterEqual(hsp_object.length, 61) #CUTOFF_GENE_LENGTH = 60
    #     self.assertEqual(len(lo_hsps), 1)
    #     self.assertIn('cj0483', lo_hsps[0].name)
    #

class TestOnePrimerFound(unittest.TestCase):
    def test_create_hsp_objects(self):
        self.assertTrue(False)

class TestPCRDirectly(unittest.TestCase):

    def test_valid_strands_f_leading_r_lagging(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = True
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = False

        f_hsp_objects = [f_hsp_object]
        r_hsp_objects = [r_hsp_object]
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object)
        self.assertTrue(f_hsp_object.valid)
        self.assertTrue(r_hsp_object.valid)

    def test_valid_strands_f_laggin_r_leading(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = False
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = True

        f_hsp_objects = [f_hsp_object]
        r_hsp_objects = [r_hsp_object]
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object)
        self.assertTrue(f_hsp_object.valid)
        self.assertTrue(r_hsp_object.valid)

    def test_valid_strands_f_leading_r_leading(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = True
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = True

        f_hsp_objects = [f_hsp_object]
        r_hsp_objects = [r_hsp_object]
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object)
        self.assertFalse(f_hsp_object.valid)
        self.assertFalse(r_hsp_object.valid)

    def test_valid_strands_f_lagging_r_lagging(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = False
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = False

        test_object = HSP('test00') #ensures only one object is removed when it is not valid
        f_hsp_objects = [f_hsp_object, test_object]
        r_hsp_objects = [r_hsp_object, test_object]
        self.assertEqual(len(f_hsp_objects), 2)
        self.assertEqual(len(r_hsp_objects), 2)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object)
        self.assertFalse(f_hsp_object.valid)
        self.assertFalse(r_hsp_object.valid)


    def test_is_distance(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

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
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)

            for f_hsp_object in forward_blast_object.hsp_objects:
                r_objects = [hsp for hsp in reverse_blast_object.hsp_objects if hsp.name in f_hsp_object.name]
                # f_hsp_object and r_hsp_object should be on the same contig
                for object in r_objects:
                    self.assertEqual(f_hsp_object.contig_name, object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, object, amplicon_sequences)
                    if "_51" in file_path:
                        self.assertFalse(distance)
                    else:
                        self.assertTrue(distance)

    def test_is_snp_primer_search(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]
        f_primers_dict = PCRPrediction.create_primer_dict(forward_primers)
        r_primers_dict = PCRPrediction.create_primer_dict(reverse_primers)

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
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)

            for f_hsp_object in forward_blast_object.hsp_objects:
                r_objects = [hsp for hsp in reverse_blast_object.hsp_objects if hsp.name in f_hsp_object.name]
                # f_hsp_object and r_hsp_object should be on the same contig
                for r_object in r_objects:
                    if PCRPrediction.pcr_directly(f_hsp_object, r_object, amplicon_sequences, f_primers_dict, r_primers_dict):
                        self.assertEqual(f_hsp_object.contig_name, r_object.contig_name)
                        f_hsp_object.snp = None
                        r_object.snp = None
                        PCRPrediction.is_snp_primer_search(f_hsp_object, f_primers_dict, r_primers_dict)
                        self.assertIsNone(r_object.snp)
                        self.assertFalse(f_hsp_object.snp)
                        PCRPrediction.is_snp_primer_search(r_object, f_primers_dict, r_primers_dict)
                        self.assertFalse(r_object.snp)

    def test_pcr_directly(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

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

        self.assertGreater(len(files_paths), 0)
        for file_path in files_paths:
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)

            self.assertGreater(len(forward_blast_object.hsp_objects), 0)
            for f_hsp_object in forward_blast_object.hsp_objects:
                r_objects = [hsp for hsp in reverse_blast_object.hsp_objects if hsp.name in f_hsp_object.name]
                self.assertGreater(len(r_objects), 0)
                # f_hsp_object and r_hsp_object should be on the same contig
                for r_object in r_objects:
                    if f_hsp_object.contig_name == r_object.contig_name:
                        # f_hsp_object and r_hsp_object should be on the same contig
                        self.assertEqual(f_hsp_object.contig_name, r_object.contig_name)
                        f_primers_dict = PCRPrediction.create_primer_dict(forward_primers)
                        r_primers_dict = PCRPrediction.create_primer_dict(reverse_primers)
                        val = PCRPrediction.pcr_directly(f_hsp_object, r_object, amplicon_sequences, f_primers_dict, r_primers_dict)
                        if "_51" in file_path:
                            self.assertFalse(val)
                        else:
                            self.assertTrue(val)

    def test_pcr_prediction_test_db(self):
        db = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        dict_predictions = PCRPrediction.main(db, forward_primers, reverse_primers, amplicon_sequences)

        predictions_fail = [pred for pred in dict_predictions if "_51" in pred]
        predictions = [pred for pred in dict_predictions if pred not in predictions_fail]
        for pred in predictions_fail:
            self.assertEqual([], dict_predictions[pred])
        for pred in predictions:
            self.assertNotEqual([], dict_predictions[pred])


class TestPCRPrediction(unittest.TestCase):
    """

    Information about files in test_pcr_prediction folder:
    #complete_primers = complete genome with 100% identity, valid strands and both primers found (pass)
    #complete_primers_max_distance = complete genome with 50bp removed from middle, 100% identity, valid strands and both primers found #Should pass if MAX_DIST_BTWN_PRIMERS <= 50 (pass)
    #complete_primers_short = complete genome with 51bp removed from middle, 100% identity, valid strands and both primers found #Should pass if MAX_DIST_BTWN_PRIMERS <= 50 (fail)
    #complete_primers_invalid_strands = complete genome with 50bp removed from middle, 100% identity, invalid strands and both primers found #Should pass if MAX_DIST_BTWN_PRIMERS <= 50 (fail)
    #one_contigs = 2 contigs with 100% identity on one contig (both primers on NODE_1) (pass)
    #db_name_different = 2 contigs wtih 100% identity on one contig (both primers on NODE_1) and an unknown name (sh0000) for db seq. (pass)
    #two_60bp_contigs = 2 contigs with 60bp, primer on different contigs, 100% identity. When CUTOFF_GENE_LENGTH >= 60, this should pass
    #two_61_bp_contigs = 2 contigs with 61bp, 100% identity. When CUTOFF_GENE_LENGTH < 61, this should pass (pass)
    #two_contigs_invalid_strands = 2 contigs with 61bp, invalid strands, 100% identity, primers found on different contigs (fail) # When CUTOFF_GENE_LENGTH < 61 this should pass
    #db_name_different_entire_gene = Each contig has 61 bp with a different name of db compared to amp seq (pass) #should produce same result as test_pcr_prediction_two_61bp_contigs

    """

    def test_pcr_prediction(self):
        db = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

        dict_predictions = PCRPrediction.main(db, forward_primers, reverse_primers,amplicon_sequences)

        predictions_fail = [pred for pred in dict_predictions if "_fail" in pred]
        predictions = [pred for pred in dict_predictions if pred not in predictions_fail]
        for pred in predictions_fail:
            self.assertEqual([], dict_predictions[pred])
        for pred in predictions:
            self.assertNotEqual([], dict_predictions[pred])
            self.assertEqual(len(dict_predictions[pred]), 1)
            for p in dict_predictions[pred]:
                self.assertEqual(len(p), 2)


class TestBPRemoved(unittest.TestCase):

    def test_create_hsp_objects_bp_removed(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        db_directory = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed"
        files = [file for file in os.listdir(db_directory) if file.endswith("fasta")]

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
            name = file_path.partition(db_directory + "/")[2]
            f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
            r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
            forward_blast_object = PCRPrediction.create_blastn_object(forward_primers, file_path, f_out_file_path)
            reverse_blast_object = PCRPrediction.create_blastn_object(reverse_primers, file_path, r_out_file_path)

            if "_ws" in file_path: #"ws"
                self.assertEqual(len(forward_blast_object.hsp_objects), 0)
                self.assertIsInstance(forward_blast_object.hsp_objects, list)

            elif "_ws" not in file_path:
                self.assertGreaterEqual(len(forward_blast_object.hsp_objects), 1)
                self.assertIsInstance(forward_blast_object.hsp_objects, list)

                if "cj0008" in file_path:
                    for hsp in forward_blast_object.hsp_objects:
                        self.assertEqual(hsp.identities, 18)
                    for hsp in reverse_blast_object.hsp_objects:
                        self.assertEqual(hsp.identities, 20)
                elif "cj0483" in file_path:
                    for hsp in forward_blast_object.hsp_objects:
                        self.assertEqual(hsp.identities, 19)
                    for hsp in reverse_blast_object.hsp_objects:
                        self.assertEqual(hsp.identities, 20)
                else:
                    self.assertFalse(True) #ensures all files are tested

            else:
                self.assertFalse(True) #ensures all files are tested

    def test_pcr_prediction_bp_removed(self):

        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        dict_predictions = PCRPrediction.main("/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed", forward_primers, reverse_primers, amplicon_sequences)

        predictions_ws = [pred for pred in dict_predictions if "_ws" in pred]
        predictions = [pred for pred in dict_predictions if pred not in predictions_ws]
        self.assertGreater(len(predictions_ws), 0)
        self.assertGreater(len(predictions), 0)
        for pred in predictions_ws:
            self.assertEqual([], dict_predictions[pred])
        for pred in predictions:
            self.assertNotEqual([], dict_predictions[pred])


# class TestMultipleHspsFound(unittest.TestCase):
#
#     # @classmethod
#     # def setUpClass(cls):
#
#     #both primer seq'n cj0008 & cj0483 found on one contig, 100% identities
#     def test_both_primer_seq_same_contig(self):
#
#         db_name_mystery = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db_entire_gene_61_bp_each_contig.fasta"
#         amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
#         forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
#         reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
#         forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
#         reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
#         full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'


#TODO: test with the names of the db name being different than amplicon sequences.

if __name__ == '__main__':
    unittest.main()