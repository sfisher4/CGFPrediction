import unittest
from Blastn import Blastn
from HSP import HSP
from Bio import SeqIO

import PCRPrediction

HSP_THRESHOLD = 0.9
E_VALUE_THRESHOLD = 0.04

    #Files available that will probably not be used
    # database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta"  # contains id and a complete genome sequence
    # full_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"
    # amp_vs_amp_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"
    # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)
    # cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta"  # complete


class TestPrimer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

        #test queries
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences

        #test databases
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
        cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
        cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
        db_name_mystery = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db.fasta"
        one_contig = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_primer_same_contig_with_2_contigs.fasta"

        #contig full
        forward_out_file_contig_full = "/home/sfisher/Sequences/blast_record/forward_primers_blast_contig_full.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file_contig_full = "/home/sfisher/Sequences/blast_record/reverse_primers_blast_contig_full.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cls._forward_blast_object_contig_full = PCRPrediction.create_blastn_object(forward_primers, cj0483_contig_full_amp_seq, forward_out_file_contig_full)
        cls._reverse_blast_object_contig_full = PCRPrediction.create_blastn_object(reverse_primers, cj0483_contig_full_amp_seq, reverse_out_file_contig_full)
        cls._lo_forward_hsp_objects_contig_full = PCRPrediction.create_hsp_objects(cls._forward_blast_object_contig_full)
        cls._lo_reverse_hsp_objects_contig_full = PCRPrediction.create_hsp_objects(cls._reverse_blast_object_contig_full)

        #contig truncated
        forward_out_file_contig_trunc = "/home/sfisher/Sequences/blast_record/forward_primers_blast_contig_trunc.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file_contig_trunc = "/home/sfisher/Sequences/blast_record/reverse_primers_blast_contig_trunc.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cls._forward_blast_object_contig_trunc = PCRPrediction.create_blastn_object(forward_primers, cj0483_contig_trunc_amp_seq, forward_out_file_contig_trunc)
        cls._reverse_blast_object_contig_trunc = PCRPrediction.create_blastn_object(reverse_primers, cj0483_contig_trunc_amp_seq, reverse_out_file_contig_trunc)
        cls._lo_forward_hsp_objects_contig_trunc = PCRPrediction.create_hsp_objects(cls._forward_blast_object_contig_trunc)
        cls._lo_reverse_hsp_objects_contig_trunc = PCRPrediction.create_hsp_objects(cls._reverse_blast_object_contig_trunc)

        #complete amp (no contigs)
        forward_out_file_complete = "/home/sfisher/Sequences/blast_record/forward_primers_blast_complete_amp.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file_complete = "/home/sfisher/Sequences/blast_record/reverse_primers_blast_complete_amp.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cls._forward_blast_object_complete = PCRPrediction.create_blastn_object(forward_primers, cj0483_complete_amp_seq, forward_out_file_complete)
        cls._reverse_blast_object_complete = PCRPrediction.create_blastn_object(reverse_primers, cj0483_complete_amp_seq, reverse_out_file_complete)
        cls._lo_forward_hsp_objects_complete = PCRPrediction.create_hsp_objects(cls._forward_blast_object_complete)
        cls._lo_reverse_hsp_objects_complete = PCRPrediction.create_hsp_objects(cls._reverse_blast_object_complete)

        #database name different then query (mystery query match)
        forward_out_file_mystery = '/home/sfisher/Sequences/blast_record/test_forward_mystery.xml'
        reverse_out_file_mystery = '/home/sfisher/Sequences/blast_record/test_reverse_mystery.xml'
        cls._forward_blast_obj_mystery = PCRPrediction.create_blastn_object(forward_primers, db_name_mystery, forward_out_file_mystery)
        cls._reverse_blast_obj_mystery = PCRPrediction.create_blastn_object(reverse_primers, db_name_mystery, reverse_out_file_mystery)
        cls._lo_forward_hsp_objects_mystery = PCRPrediction.create_hsp_objects(cls._forward_blast_obj_mystery)
        cls._lo_reverse_hsp_objects_mystery = PCRPrediction.create_hsp_objects(cls._reverse_blast_obj_mystery)

        #two contigs with both primers on one contig
        forward_out_file_one_contigs = '/home/sfisher/Sequences/blast_record/test_one_contigs_forward.xml'
        reverse_out_file_one_contigs = '/home/sfisher/Sequences/blast_record/test_one_contigs_reverse.xml'
        cls._forward_blast_obj_one_contigs = PCRPrediction.create_blastn_object(forward_primers, one_contig, forward_out_file_one_contigs)
        cls._reverse_blast_obj_one_contigs = PCRPrediction.create_blastn_object(reverse_primers, one_contig, reverse_out_file_one_contigs)
        cls._lo_forward_hsp_objects_one_contigs = PCRPrediction.create_hsp_objects(cls._forward_blast_obj_one_contigs)
        cls._lo_reverse_hsp_objects_one_contigs = PCRPrediction.create_hsp_objects(cls._reverse_blast_obj_one_contigs)

    #contig trunc
    # Tests to make sure the blastn_object is initialized when called with a primer as query and an amp sequence as db
    def test_create_blastn_primer_search_object(self):
        self.assertNotEqual([], self._forward_blast_object_contig_trunc.blast_records)
        self.assertNotEqual({}, self._forward_blast_object_contig_trunc.hsp_records)
        self.assertGreaterEqual(len(self._forward_blast_object_contig_trunc.blast_records), 1)
        self.assertGreaterEqual(len(self._forward_blast_object_contig_trunc.hsp_records), 1)

    #complete amp & full gene search
    #Tests to make sure the blastn object is initialized when called with amplicon sequences as query and an amp sequence as db
    def test_create_blastn_full_search_object(self):
        self.assertNotEqual([], self._forward_blast_object_complete.blast_records)
        self.assertNotEqual({}, self._forward_blast_object_complete.hsp_records)
        self.assertGreaterEqual(len(self._forward_blast_object_complete.blast_records), 1)
        self.assertGreaterEqual(len(self._forward_blast_object_complete.hsp_records), 1)

    def test_create_blastn_mystery_query_match(self):
        self.assertNotEqual([], self._forward_blast_obj_mystery.blast_records, "Blast records are empty")
        self.assertNotEqual({}, self._forward_blast_obj_mystery.hsp_records, "Hsp records are empty")
        self.assertEqual(len(self._forward_blast_obj_mystery.blast_records), 40, "Blast records are not being created for all 40 primers")

    #Contig full
    def test_create_hsp_records_primer_search(self):

        self.assertGreaterEqual(len(self._reverse_blast_object_contig_full.hsp_records.keys()), 1)

        for blast_record in self._reverse_blast_object_contig_full.blast_records:
            for alignment in blast_record.alignments:
                if alignment in self._reverse_blast_object_contig_full.hsp_records: #ensures the hsp is in the alignment
                    self.assertIsInstance(self._reverse_blast_object_contig_full.hsp_records[alignment], list)
                    for hsp in self._reverse_blast_object_contig_full.hsp_records[alignment]:
                        self.assertTrue(HSP_THRESHOLD <= (hsp.identities / 21)) #the length of primer seq in this case is 21
                    self.assertIsInstance(self._reverse_blast_object_contig_full.hsp_records, dict)

    # complete gene on one contig
    # mystery contains 2 contigs with the entire gene on one contig (NODE_1)
    #should produce the same result as one_contigs
    def test_create_hsp_records_mystery(self):
        self.assertEqual(len(self._forward_blast_obj_mystery.blast_records), 40)
        print(self._forward_blast_obj_mystery.hsp_records)

        self.assertGreaterEqual(len(self._reverse_blast_obj_mystery.hsp_records), 1)
        self.assertGreaterEqual(len(self._forward_blast_obj_mystery.hsp_records), 1)

    #complete gene on one contigs
    def test_create_hsp_records_one_contigs(self):
        self.assertEqual(len(self._forward_blast_obj_one_contigs.blast_records), 40)
        print(self._forward_blast_obj_one_contigs.hsp_records)

        self.assertGreaterEqual(len(self._forward_blast_obj_one_contigs.hsp_records), 1)
        self.assertGreaterEqual(len(self._reverse_blast_obj_one_contigs.hsp_records), 1)

    #complete gene on one contig
    #mystery contains 2 contigs with the entire gene on one contig (NODE_1)
    #should produce the same result as one_contigs db
    def test_create_hsp_objects_mystery(self):
        self.assertEqual(len(self._lo_forward_hsp_objects_mystery), 1)
        self.assertEqual(len(self._lo_reverse_hsp_objects_mystery), 1)

        # for hsp_object in self._lo_reverse_hsp_objects_mystery:
        #     self.assertEqual(hsp_object.name, 'cj0483')

    #two contigs with gene on one of the contigs (NODE_1)
    def test_create_one_contigs_hsp_objects(self):
        self.assertEqual(len(self._lo_forward_hsp_objects_one_contigs), 1)
        self.assertEqual(len(self._lo_reverse_hsp_objects_one_contigs), 1)

        # for hsp_object in self._lo_reverse_hsp_objects_one_contigs:
        #     self.assertEqual(hsp_object.name, 'cj0483')

    # contig full
    def test_create_hsp_objects(self):

        self.assertEqual(len(self._lo_reverse_hsp_objects_contig_full), 1)
        self.assertEqual(len(self._lo_reverse_hsp_objects_contig_full), 1)

        # test name
        test = False
        for hsp_object in self._lo_reverse_hsp_objects_contig_full:
            # tests name
            self.assertEqual(hsp_object.name, 'cj0483')
            # tests db_length
            if 'NODE_1' in hsp_object.contig_name:
                self.assertEqual(hsp_object.db_length, 376)
            elif 'NODE_2' in hsp_object.contig_name:
                self.assertEqual(hsp_object.db_length, 236)
                test = True  # make sure both contigs are reached
            else:
                self.assertWarns("shouldn't reach here!!!")

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
                             abs(hsp_object.end - hsp_object.start) + 1)  # TODO: note: changed in code.
            self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)

        self.assertTrue(test)  # makes sure both contigs are reached (db_length test)

    #complete amp
    #Note: MX_DIST_BTWN_PRIEMRS = 50
    def test_is_distance_exact(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete

        for f_hsp_object in self._lo_forward_hsp_objects_complete:
            for r_hsp_object in self._lo_reverse_hsp_objects_complete:
                if f_hsp_object.name == r_hsp_object.name:
                    #f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, amplicon_sequences, cj0483_complete_amp_seq)
                    self.assertTrue(distance)

class TestEntireGene(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        #query
        cls._amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

        #databases
        both_contigs_306_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_306_bp.fasta"
        both_contigs_27_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_27_bp.fasta"
        both_contigs_60_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_60_bp.fasta"
        both_contigs_61_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_entire_gene_search/cj0483_both_contigs_61_bp.fasta"
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
        cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
        cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
        both_contigs_61bp_db_name_mystery = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db_entire_gene_61_bp_each_contig.fasta"

        #out files
        out_file_306 = '/home/sfisher/Sequences/blast_record/test_blast_306.xml'
        out_file_27 = '/home/sfisher/Sequences/blast_record/test_blast_27.xml'
        out_file_60 = '/home/sfisher/Sequences/blast_record/test_blast_60.xml'
        out_file_61 = '/home/sfisher/Sequences/blast_record/test_blast_61.xml'
        out_file_contig_full = '/home/sfisher/Sequences/blast_record/test_blast_contig_full.xml'
        out_file_contig_trunc = '/home/sfisher/Sequences/blast_record/test_blast_contig_trunc.xml'
        out_file_complete = '/home/sfisher/Sequences/blast_record/test_blast_complete.xml'
        out_file_mystery = '/home/sfisher/Sequences/blast_record/test_blast_mystery_61bp.xml'

        #contig full
        cls._blast_object_contig_full = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, cj0483_contig_full_amp_seq, out_file_contig_full)
        cls._lo_hsp_objects_contig_full = PCRPrediction.create_hsp_objects(cls._blast_object_contig_full)

        #contig truncated
        cls._blast_object_contig_trunc = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, cj0483_contig_trunc_amp_seq, out_file_contig_trunc)
        cls._lo_hsp_objects_contig_trunc = PCRPrediction.create_hsp_objects(cls._blast_object_contig_trunc)

        #complete amp (no contigs)
        cls._blast_object_complete = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, cj0483_complete_amp_seq, out_file_complete)
        cls._lo_hsp_objects_complete = PCRPrediction.create_hsp_objects(cls._blast_object_complete)

        #both contigs length 306 bp
        cls._blast_obj_306= PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, both_contigs_306_bp, out_file_306)
        cls._lo_hsp_objects_306 = PCRPrediction.create_hsp_objects(cls._blast_obj_306)

        #both contigs length 27 bp
        cls._blast_obj_27 = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, both_contigs_27_bp, out_file_27)
        cls._lo_hsp_objects_27 = PCRPrediction.create_hsp_objects(cls._blast_obj_27)

        #both contigs length 60 bp
        cls._blast_obj_60 = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, both_contigs_60_bp, out_file_60)
        cls._lo_hsp_objects_60 = PCRPrediction.create_hsp_objects(cls._blast_obj_60)

        #both contigs length 61 bp
        cls._blast_obj_61 = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, both_contigs_61_bp, out_file_61)
        cls._lo_hsp_objects_61 = PCRPrediction.create_hsp_objects(cls._blast_obj_61)

        #both contigs length 61 bp with a different db name then amplicon sequences
        cls._blast_obj_61_mystery = PCRPrediction.create_blastn_full_search_object(cls._amplicon_sequences, both_contigs_61bp_db_name_mystery, out_file_mystery)
        cls._lo_hsp_objects_61_mystery = PCRPrediction.create_hsp_objects(cls._blast_obj_61_mystery)

    #contigs trunc
    #search full gene
    def test_create_blast_records_full_search_contigs_truncated(self):
        self.assertNotEqual([], self._blast_object_contig_trunc.blast_records)
        self.assertEqual(len(self._blast_object_contig_trunc.blast_records), 40)

    #complete amp
    #search full gene
    def test_create_hsp_records_full_search_contigs_complete(self):

        self.assertGreaterEqual(len(self._blast_object_complete.hsp_records), 1)

        for blast_record in self._blast_object_complete.blast_records:
            for alignment in blast_record.alignments:
                if alignment in self._blast_object_complete.hsp_records: #ensures the hsp is in the alignment
                    self.assertIsInstance(self._blast_object_complete.hsp_records[alignment], list)
                    for hsp in self._blast_object_complete.hsp_records[alignment]:
                        self.assertTrue(HSP_THRESHOLD <= (hsp.identities / 612)) #the length of gene seq in this case is 612 (cj0483_complete_amp_seq)
                    self.assertIsInstance(self._blast_object_complete.hsp_records, dict)

    # contigs trunc
    # full gene search
    def test_create_hsp_records_full_search_contigs(self):

        self.assertGreaterEqual(len(self._blast_object_contig_trunc.hsp_records), 1)
        self.assertGreaterEqual(len(self._blast_object_contig_trunc.hsp_records), 2)

        for blast_record in self._blast_object_contig_trunc.blast_records:
            for alignment in blast_record.alignments:
                if alignment in self._blast_object_contig_trunc.hsp_records:  # ensures the hsp is in the alignment
                    self.assertIsInstance(self._blast_object_contig_trunc.hsp_records[alignment], list)
                    for hsp in self._blast_object_contig_trunc.hsp_records[alignment]:
                        if "NODE_1" in alignment.hit_def:
                            self.assertTrue(HSP_THRESHOLD <= hsp.identities / 376)
                        elif "NODE_2" in alignment.hit_def:
                            self.assertTrue(HSP_THRESHOLD <= hsp.identities / 228)
                        else:
                            self.assertFalse(True, 'should never reach here unless more than two nodes or error')
                    self.assertIsInstance(self._blast_object_contig_trunc.hsp_records, dict)

    #contigs trunc (61 each contig)
    #full gene search
    def test_create_hsp_records_full_search_contigs_mystery(self):


        self.assertGreaterEqual(len(self._blast_obj_61_mystery.hsp_records), 1)
        self.assertGreaterEqual(len(self._blast_obj_61_mystery.hsp_records), 2)

        for blast_record in self._blast_obj_61_mystery.blast_records:
            for alignment in blast_record.alignments:
                if alignment in self._blast_obj_61_mystery.hsp_records:  # ensures the hsp is in the alignment
                    self.assertIsInstance(self._blast_obj_61_mystery.hsp_records[alignment], list)
                    for hsp in self._blast_obj_61_mystery.hsp_records[alignment]:
                            self.assertTrue(HSP_THRESHOLD <= hsp.identities / 61)
                    self.assertIsInstance(self._blast_obj_61_mystery.hsp_records, dict)


    #contigs trunc
    #search full gene
    def test_create_hsp_objects_contig_trunc(self):
        self.assertGreaterEqual(len(self._lo_hsp_objects_contig_trunc), 2)
        test = False
        for hsp_object in self._lo_hsp_objects_contig_trunc:
            # tests name
            self.assertIn('cj0483', hsp_object.name)
            # tests db_length
            if 'NODE_1' in hsp_object.contig_name:
                self.assertEqual(hsp_object.db_length, 376)
            elif 'NODE_2' in hsp_object.contig_name:
                self.assertEqual(hsp_object.db_length, 228)
                test = True  # make sure both contigs are reached
            else:
                self.assertWarns("shouldn't reach here!!!")

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
                             abs(hsp_object.end - hsp_object.start) + 1)  # TODO: note: changed in code.
            self.assertIsNone(hsp_object.valid)  # have not yet intialized valid (do in isValid function)

        self.assertTrue(test)  # makes sure both contigs are reached (db_length test)


    #306 bp on each contig
    #TODO: finish testing me
    def test_entire_gene_306(self):
        self.assertGreaterEqual(len(self._lo_hsp_objects_306), 2)
        for hsp_object in self._lo_hsp_objects_306:
            lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_306, hsp_object)
            self.assertTrue(len(lo_queries) == 2)
            for query in lo_queries:
                self.assertEqual(hsp_object.name, query.name)
                self.assertTrue(query.valid)

    #27 bp on each contig
    #TODO: finish testing me
    def test_entire_gene_27(self):
        self.assertGreaterEqual(len(self._lo_hsp_objects_27), 2)
        for hsp_object in self._lo_hsp_objects_27:
            lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_27, hsp_object)
            self.assertTrue(len(lo_queries) == 0)

    def test_entire_gene_60(self):
        self.assertTrue(len(self._lo_hsp_objects_60) >= 2)
        for hsp_object in self._lo_hsp_objects_60:
            lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_60, hsp_object)
            self.assertEqual(len(lo_queries), 0)

    def test_entire_gene_61(self):
        self.assertTrue(len(self._lo_hsp_objects_61) >= 2)
        for hsp_object in self._lo_hsp_objects_61:
            lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_61, hsp_object)
            self.assertTrue(len(lo_queries) == 2)
            for query in lo_queries:
                self.assertEqual(hsp_object.name, query.name)
                self.assertTrue(query.valid)

    def test_entire_gene_61_db_name_different(self):
        self.assertGreaterEqual(len(self._lo_hsp_objects_61_mystery), 2)
        for hsp_object in self._lo_hsp_objects_61_mystery:
            lo_queries = PCRPrediction.entire_gene(self._lo_hsp_objects_61_mystery, hsp_object)
            self.assertTrue(len(lo_queries) == 2)
            for query in lo_queries:
                self.assertEqual(hsp_object.name, query.name)
                self.assertTrue(query.valid)


class TestPCRDirectly(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences

        # 20 mid bp removed
        forward_out_file_20 = "/home/sfisher/Sequences/blast_record/forward_primers_blast_20_rm.xml"
        reverse_out_file_20 = "/home/sfisher/Sequences/blast_record/reverse_primers_blast_20_rm.xml"
        cj0483_20_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_20_middle_bp_removed.fasta"
        cls._forward_blast_object_20 = PCRPrediction.create_blastn_object(forward_primers, cj0483_20_mid_bp_rm_amp_seq, forward_out_file_20)
        cls._reverse_blast_object_20 = PCRPrediction.create_blastn_object(reverse_primers, cj0483_20_mid_bp_rm_amp_seq, reverse_out_file_20)
        cls._lo_forward_hsp_objects_20 = PCRPrediction.create_hsp_objects(cls._forward_blast_object_20)
        cls._lo_reverse_hsp_objects_20 = PCRPrediction.create_hsp_objects(cls._reverse_blast_object_20)

        #50 mid bp removed
        forward_out_file_50 = "/home/sfisher/Sequences/blast_record/forward_primers_blast_50_rm.xml"
        reverse_out_file_50 = "/home/sfisher/Sequences/blast_record/reverse_primers_blast_50_rm.xml"
        cj0483_50_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_50_middle_bp_removed.fasta"
        cls._forward_blast_object_50 = PCRPrediction.create_blastn_object(forward_primers, cj0483_50_mid_bp_rm_amp_seq, forward_out_file_50)
        cls._reverse_blast_object_50 = PCRPrediction.create_blastn_object(reverse_primers, cj0483_50_mid_bp_rm_amp_seq, reverse_out_file_50)
        cls._lo_forward_hsp_objects_50 = PCRPrediction.create_hsp_objects(cls._forward_blast_object_50)
        cls._lo_reverse_hsp_objects_50 = PCRPrediction.create_hsp_objects(cls._reverse_blast_object_50)

        #51 mid bp removed
        forward_out_file_51 = "/home/sfisher/Sequences/blast_record/forward_primers_blast_51_rm.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file_51 = "/home/sfisher/Sequences/blast_record/reverse_primers_blast_51_rm.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_51_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_51_middle_bp_removed.fasta"
        cls._forward_blast_51 = PCRPrediction.create_blastn_object(forward_primers, cj0483_51_mid_bp_rm_amp_seq, forward_out_file_51)
        cls._reverse_blast_51 = PCRPrediction.create_blastn_object(reverse_primers, cj0483_51_mid_bp_rm_amp_seq, reverse_out_file_51)
        cls._lo_forward_hsp_objects_51 = PCRPrediction.create_hsp_objects(cls._forward_blast_51)
        cls._lo_reverse_hsp_objects_51 = PCRPrediction.create_hsp_objects(cls._reverse_blast_51)

    def test_valid_strands_f_leading_r_lagging(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = True
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = False

        f_hsp_objects = [f_hsp_object]
        r_hsp_objects = [r_hsp_object]
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object, f_hsp_objects, r_hsp_objects)
        self.assertTrue(f_hsp_object.valid)
        self.assertTrue(r_hsp_object.valid)
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

    def test_valid_strands_f_laggin_r_leading(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = False
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = True

        f_hsp_objects = [f_hsp_object]
        r_hsp_objects = [r_hsp_object]
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object, f_hsp_objects, r_hsp_objects)
        self.assertTrue(f_hsp_object.valid)
        self.assertTrue(r_hsp_object.valid)
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

    def test_valid_strands_f_leading_r_leading(self):
        f_hsp_object = HSP('sh0000')
        f_hsp_object.strand = True
        r_hsp_object = HSP('sh0000')
        r_hsp_object.strand = True

        f_hsp_objects = [f_hsp_object]
        r_hsp_objects = [r_hsp_object]
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object, f_hsp_objects, r_hsp_objects)
        self.assertFalse(f_hsp_object.valid)
        self.assertFalse(r_hsp_object.valid)
        self.assertEqual(len(f_hsp_objects), 0)
        self.assertEqual(len(r_hsp_objects), 0)

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

        PCRPrediction.valid_strands(f_hsp_object, r_hsp_object, f_hsp_objects, r_hsp_objects)
        self.assertFalse(f_hsp_object.valid)
        self.assertFalse(r_hsp_object.valid)
        self.assertEqual(len(f_hsp_objects), 1)
        self.assertEqual(len(r_hsp_objects), 1)

    #20 mid bp removed
    #Note: MX_DIST_BTWN_PRIEMRS = 50
    def test_is_distance_20_bp_rm(self):
        cj0483_20_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_20_middle_bp_removed.fasta"  # complete

        for f_hsp_object in self._lo_forward_hsp_objects_20:
            for r_hsp_object in self._lo_reverse_hsp_objects_20:
                if f_hsp_object.name == r_hsp_object.name:
                    #f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, self._amplicon_sequences, cj0483_20_mid_bp_rm_amp_seq)
                    self.assertTrue(distance)

    #50 mid bp removed
    #Note: MAX_DIST_BTWN_PRIMERS = 50
    def test_is_distance_50_bp_rm(self):
        cj0483_50_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_50_middle_bp_removed.fasta"  # complete

        for f_hsp_object in self._lo_forward_hsp_objects_50:
            for r_hsp_object in self._lo_reverse_hsp_objects_50:
                if f_hsp_object.name == r_hsp_object.name:
                    # f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, self._amplicon_sequences, cj0483_50_mid_bp_rm_amp_seq)
                    self.assertTrue(distance)

    #51 mid bp removed
    #Note: MAX_DIST_BTWN_PRIMERS = 50
    def test_is_distance_51_bp_rm(self):
        cj0483_51_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_51_middle_bp_removed.fasta"

        for f_hsp_object in self._lo_forward_hsp_objects_51:
            for r_hsp_object in self._lo_reverse_hsp_objects_51:
                if f_hsp_object.name == r_hsp_object.name:
                    # f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, self._amplicon_sequences, cj0483_51_mid_bp_rm_amp_seq)
                    self.assertFalse(distance)

    #50 mid bp removed
    def test_pcr_directly_dist_50(self):
        cj0483_50_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_50_middle_bp_removed.fasta"

        for f_hsp_object in self._lo_forward_hsp_objects_50:
            for r_hsp_object in self._lo_reverse_hsp_objects_50:
                if f_hsp_object.name == r_hsp_object.name:
                    if f_hsp_object.contig_name == r_hsp_object.contig_name:
                        # f_hsp_object and r_hsp_object should be on the same contig
                        self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                        val = PCRPrediction.pcr_directly(f_hsp_object, r_hsp_object, cj0483_50_mid_bp_rm_amp_seq, self._amplicon_sequences, self._lo_forward_hsp_objects_50, self._lo_reverse_hsp_objects_50)
                        self.assertTrue(val)

    #30 mid bp removed
    def test_pcr_directly_dist_20(self):
        cj0483_20_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_20_middle_bp_removed.fasta"  # complete

        for f_hsp_object in self._lo_forward_hsp_objects_20:
            for r_hsp_object in self._lo_reverse_hsp_objects_20:
                if f_hsp_object.name == r_hsp_object.name:
                    if f_hsp_object.contig_name == r_hsp_object.contig_name:
                        # f_hsp_object and r_hsp_object should be on the same contig
                        self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                        val = PCRPrediction.pcr_directly(f_hsp_object, r_hsp_object, cj0483_20_mid_bp_rm_amp_seq, self._amplicon_sequences, self._lo_forward_hsp_objects_20, self._lo_reverse_hsp_objects_20)
                        self.assertTrue(val)


class TestPCRPrediction(unittest.TestCase):

    def test_pcr_prediction_complete_primers(self):
        complete = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_complete.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
        reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
        full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'

        lo_queries = PCRPrediction.pcr_prediction(forward_primers, reverse_primers, complete, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file)
        self.assertTrue(len(lo_queries) == 1)
        for queries in lo_queries:
            self.assertTrue(len(queries) == 2)
            for query in queries:
                self.assertTrue(query.name == 'cj0483')
                self.assertIn('cj0483', query.contig_name)

    def test_pcr_prediction_one_contigs(self):
        one_contig = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_primer_same_contig_with_2_contigs.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
        reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
        full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'

        lo_queries = PCRPrediction.pcr_prediction(forward_primers, reverse_primers, one_contig, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file)
        self.assertTrue(len(lo_queries) == 1)
        for queries in lo_queries:
            self.assertTrue(len(queries) == 2)
            self.assertEqual(queries[0].contig_name, queries[1].contig_name)
            for query in queries:
                self.assertEqual(query.name, 'cj0483')

    def test_pcr_prediction_db_name_different(self):
        db_name_mystery = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
        reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
        full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'

        lo_queries = PCRPrediction.pcr_prediction(forward_primers, reverse_primers, db_name_mystery, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file)
        self.assertTrue(len(lo_queries) == 1)
        for queries in lo_queries:
            self.assertTrue(len(queries) == 2)
            self.assertEqual(queries[0].contig_name, queries[1].contig_name)
            for query in queries:
                self.assertEqual(query.name, 'cj0483')
                self.assertIn('sh0000', query.contig_name)


    def test_pcr_prediction_two_60bp_contigs(self):
        both_contigs_60_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_both_contigs_60_bp.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
        reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
        full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'

        lo_queries = PCRPrediction.pcr_prediction(forward_primers, reverse_primers, both_contigs_60_bp, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file)
        self.assertTrue(len(lo_queries) <= 0)

    def test_pcr_prediction_two_61bp_contigs(self):
        both_contigs_61_bp = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_both_contigs_61_bp.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
        reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
        full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'

        lo_queries = PCRPrediction.pcr_prediction(forward_primers, reverse_primers, both_contigs_61_bp, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file)
        self.assertTrue(len(lo_queries) == 1)
        for queries in lo_queries:
            self.assertEqual(len(queries), 2)
            self.assertNotEqual(queries[0].contig_name, queries[1].contig_name)
            for query in queries:
                self.assertIn('cj0483', query.name)

    #Each contig has 61 bp with a different name of db compared to amp seq
    #should produce same result as test_pcr_prediction_two_61bp_contigs
    def test_pcr_prediction_db_name_different_entire_gene(self):
        db_name_mystery = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db_entire_gene_61_bp_each_contig.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = '/home/sfisher/Sequences/blast_record/forward_blast.xml'
        reverse_out_file = '/home/sfisher/Sequences/blast_record/reverse_blast.xml'
        full_out_file = '/home/sfisher/Sequences/blast_record/full_out.xml'
        lo_queries = PCRPrediction.pcr_prediction(forward_primers, reverse_primers, db_name_mystery,
                                                  forward_out_file, reverse_out_file, amplicon_sequences,
                                                  full_out_file)
        self.assertEqual(len(lo_queries), 1)
        for queries in lo_queries:
            self.assertEqual(len(queries), 2)
            self.assertNotEqual(queries[0].contig_name, queries[1].contig_name)
            for query in queries:
                self.assertIn('cj0483', query.name)



#TODO: test with the names of the db name being different than amplicon sequences.

if __name__ == '__main__':
    unittest.main()