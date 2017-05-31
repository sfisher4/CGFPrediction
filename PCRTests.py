import unittest
from Blastn import Blastn
from HSP import HSP
from Bio import SeqIO

import PCRPrediction

HSP_THRESHOLD = 0.9
E_VALUE_THRESHOLD = 0.04

class MyTestCase(unittest.TestCase):

    database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta"  # contains id and a complete genome sequence
    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
    forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
    reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
    full_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"
    amp_vs_amp_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"
    test_out_file = '/home/sfisher/Sequences/blast_record/test_blast.xml'

    cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta"  # complete
    cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
    cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
    cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"

    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

    # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)

    lo_query_files = [amplicon_sequences, forward_primers, reverse_primers]
    lo_database_files = [database, cj0483_contig_trunc_amp_seq, cj0483_contig_full_amp_seq, cj0483_complete_amp_seq, cj0008_amp_seq]
    lo_out_files = [forward_out_file, reverse_out_file, full_out_file]

    #Tests to make sure the blastn_object is initialized when called with a primer as query and an amp sequence as db
    def test_create_blastn_primer_search_object(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
        cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
        cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
        test_out_file = '/home/sfisher/Sequences/blast_record/test_blast.xml'
        name = "forward_primers vs complete amp"

        blastn_object = PCRPrediction.create_blastn_object(forward_primers, cj0483_contig_trunc_amp_seq, test_out_file, name)
        self.assertEqual(name, blastn_object.name)
        self.assertNotEqual([], blastn_object.blast_records)
        self.assertNotEqual({}, blastn_object.hsp_records)
        self.assertGreaterEqual(len(blastn_object.blast_records), 1)
        self.assertGreaterEqual(len(blastn_object.hsp_records), 1)

    #Tests to make sure the blastn object is initialized when called with amplicon sequences as query and an amp sequence as db
    def test_create_blastn_full_search_object(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        test_out_file = '/home/sfisher/Sequences/blast_record/test_blast.xml'
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
        name = "amps vs amp"

        blastn_object = PCRPrediction.create_blastn_full_search_object(amplicon_sequences, cj0483_complete_amp_seq, test_out_file, name)
        self.assertEqual(name, blastn_object.name)
        self.assertNotEqual([], blastn_object.blast_records)
        self.assertNotEqual({}, blastn_object.hsp_records)
        self.assertGreaterEqual(len(blastn_object.blast_records), 1)
        self.assertGreaterEqual(len(blastn_object.hsp_records), 1)

    def test_create_hsp_records_primer_search(self):
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
        cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
        cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
        test_out_file = '/home/sfisher/Sequences/blast_record/test_blast.xml'
        name = "primer vs amp"

        blastn_object = PCRPrediction.create_blastn_object(reverse_primers, cj0483_contig_full_amp_seq, test_out_file, name)

        self.assertGreaterEqual(len(blastn_object.hsp_records.keys()), 1)

        for blast_record in blastn_object.blast_records:
            for alignment in blast_record.alignments:
                if alignment in blastn_object.hsp_records: #ensures the hsp is in the alignment
                    self.assertIsInstance(blastn_object.hsp_records[alignment], list)
                    for hsp in blastn_object.hsp_records[alignment]:
                        self.assertTrue(HSP_THRESHOLD <= (hsp.identities / 21)) #the length of primer seq in this case is 21
                    self.assertIsInstance(blastn_object.hsp_records, dict)


    def test_create_hsp_records_full_search(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        test_out_file = '/home/sfisher/Sequences/blast_record/test_blast.xml'
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete
        name = "amps vs amp"

        blastn_object = PCRPrediction.create_blastn_full_search_object(amplicon_sequences, cj0483_complete_amp_seq, test_out_file, name)

        self.assertGreaterEqual(len(blastn_object.hsp_records.keys()), 1)

        for blast_record in blastn_object.blast_records:
            for alignment in blast_record.alignments:
                if alignment in blastn_object.hsp_records: #ensures the hsp is in the alignment
                    self.assertIsInstance(blastn_object.hsp_records[alignment], list)
                    for hsp in blastn_object.hsp_records[alignment]:
                        self.assertTrue(HSP_THRESHOLD <= (hsp.identities / 612)) #the length of gene seq in this case is 612 (cj0483_complete_amp_seq)
                    self.assertIsInstance(blastn_object.hsp_records, dict)

    def test_create_hsp_objects(self):
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
        test_out_file = '/home/sfisher/Sequences/blast_record/test_blast.xml'
        name = "primer vs amp"

        #setup
        blastn_object = PCRPrediction.create_blastn_object(reverse_primers, cj0483_contig_full_amp_seq, test_out_file, name)

        #tests
        hsp_objects = PCRPrediction.create_hsp_objects(blastn_object)
        self.assertGreaterEqual(len(hsp_objects), 1)

        #test name

        test = False
        for hsp_object in hsp_objects:
            #tests name
            self.assertEqual(hsp_object.name, 'cj0483')
            #tests db_length
            if 'NODE_1' in hsp_object.contig_name:
                self.assertEqual(hsp_object.db_length, 376)
            elif 'NODE_2' in hsp_object.contig_name:
                self.assertEqual(hsp_object.db_length, 236)
                test = True #make sure both contigs are reached
            else:
                self.assertWarns("shouldn't reach here!!!")

            self.assertNotEqual(hsp_object.contig_name, "", 'hsp_object contig name is not initialized')
            self.assertGreater(E_VALUE_THRESHOLD, hsp_object.expect)
            self.assertNotEqual(hsp_object.expect, -1)
            self.assertNotEqual(hsp_object.start, -1)
            self.assertNotEqual(hsp_object.end, -1)
            self.assertNotEqual(hsp_object.start, hsp_object.end)
            if hsp_object.start < hsp_object.end :
                self.assertTrue(hsp_object.strand)
            if hsp_object.end < hsp_object.start:
                self.assertFalse(hsp_object.strand)
            self.assertIsNotNone(hsp_object.strand)
            self.assertNotEqual(hsp_object.length, 0, "length of hsp is 0")
            self.assertEqual(hsp_object.length, abs(hsp_object.end - hsp_object.start) + 1) #TODO: note: changed in code.
            self.assertIsNone(hsp_object.valid) #have not yet intialized valid (do in isValid function)

        self.assertTrue(test) #makes sure both contigs are reached (db_length test)
    #
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

    #Note: MX_DIST_BTWN_PRIEMRS = 50
    def test_is_distance_exact(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta"  # complete

        forward_blast = PCRPrediction.create_blastn_object(forward_primers, cj0483_complete_amp_seq, forward_out_file, 'forward')
        reverse_blast = PCRPrediction.create_blastn_object(reverse_primers, cj0483_complete_amp_seq, reverse_out_file, 'reverse')

        lo_forward_hsp_objects = PCRPrediction.create_hsp_objects(forward_blast)
        lo_reverse_hsp_objects = PCRPrediction.create_hsp_objects(reverse_blast)

        for f_hsp_object in lo_forward_hsp_objects:
            for r_hsp_object in lo_reverse_hsp_objects:
                if f_hsp_object.name == r_hsp_object.name:
                    #f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, amplicon_sequences, cj0483_complete_amp_seq)
                    self.assertTrue(distance)

    #Note: MX_DIST_BTWN_PRIEMRS = 50
    def test_is_distance_20_bp_rm(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_20_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_20_middle_bp_removed.fasta"  # complete

        forward_blast = PCRPrediction.create_blastn_object(forward_primers, cj0483_20_mid_bp_rm_amp_seq, forward_out_file, 'forward')
        reverse_blast = PCRPrediction.create_blastn_object(reverse_primers, cj0483_20_mid_bp_rm_amp_seq, reverse_out_file, 'reverse')

        lo_forward_hsp_objects = PCRPrediction.create_hsp_objects(forward_blast)
        lo_reverse_hsp_objects = PCRPrediction.create_hsp_objects(reverse_blast)

        for f_hsp_object in lo_forward_hsp_objects:
            for r_hsp_object in lo_reverse_hsp_objects:
                if f_hsp_object.name == r_hsp_object.name:
                    #f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, amplicon_sequences, cj0483_20_mid_bp_rm_amp_seq)
                    self.assertTrue(distance)

    #Note: MAX_DIST_BTWN_PRIMERS = 50
    def test_is_distance_50_bp_rm(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_50_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_50_middle_bp_removed.fasta"  # complete

        forward_blast = PCRPrediction.create_blastn_object(forward_primers, cj0483_50_mid_bp_rm_amp_seq,
                                                           forward_out_file, 'forward')
        reverse_blast = PCRPrediction.create_blastn_object(reverse_primers, cj0483_50_mid_bp_rm_amp_seq,
                                                           reverse_out_file, 'reverse')

        lo_forward_hsp_objects = PCRPrediction.create_hsp_objects(forward_blast)
        lo_reverse_hsp_objects = PCRPrediction.create_hsp_objects(reverse_blast)

        for f_hsp_object in lo_forward_hsp_objects:
            for r_hsp_object in lo_reverse_hsp_objects:
                if f_hsp_object.name == r_hsp_object.name:
                    # f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, amplicon_sequences, cj0483_50_mid_bp_rm_amp_seq)
                    self.assertTrue(distance)

    #Note: MAX_DIST_BTWN_PRIMERS = 50
    def test_is_distance_51_bp_rm(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_51_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_51_middle_bp_removed.fasta"  # complete

        forward_blast = PCRPrediction.create_blastn_object(forward_primers, cj0483_51_mid_bp_rm_amp_seq, forward_out_file, 'forward')
        reverse_blast = PCRPrediction.create_blastn_object(reverse_primers, cj0483_51_mid_bp_rm_amp_seq, reverse_out_file, 'reverse')

        lo_forward_hsp_objects = PCRPrediction.create_hsp_objects(forward_blast)
        lo_reverse_hsp_objects = PCRPrediction.create_hsp_objects(reverse_blast)

        for f_hsp_object in lo_forward_hsp_objects:
            for r_hsp_object in lo_reverse_hsp_objects:
                if f_hsp_object.name == r_hsp_object.name:
                    # f_hsp_object and r_hsp_object should be on the same contig
                    self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                    distance = PCRPrediction.is_distance(f_hsp_object, r_hsp_object, amplicon_sequences, cj0483_51_mid_bp_rm_amp_seq)
                    self.assertFalse(distance)

    def test_pcr_directly_dist_50(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_50_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_50_middle_bp_removed.fasta"  # complete

        forward_blast = PCRPrediction.create_blastn_object(forward_primers, cj0483_50_mid_bp_rm_amp_seq,
                                                           forward_out_file, 'forward')
        reverse_blast = PCRPrediction.create_blastn_object(reverse_primers, cj0483_50_mid_bp_rm_amp_seq,
                                                           reverse_out_file, 'reverse')

        lo_forward_hsp_objects = PCRPrediction.create_hsp_objects(forward_blast)
        lo_reverse_hsp_objects = PCRPrediction.create_hsp_objects(reverse_blast)

        for f_hsp_object in lo_forward_hsp_objects:
            for r_hsp_object in lo_reverse_hsp_objects:
                if f_hsp_object.name == r_hsp_object.name:
                    if f_hsp_object.contig_name == r_hsp_object.contig_name:
                        # f_hsp_object and r_hsp_object should be on the same contig
                        self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                        val = PCRPrediction.pcr_directly(f_hsp_object, r_hsp_object, cj0483_50_mid_bp_rm_amp_seq, amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects)
                        self.assertTrue(val)

    def test_pcr_directly_dist_20(self):
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"  # contains primer id's and primer sequences
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"  # contains primer id's and primer sequences
        forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml"  # location where the blast record from comparing the forward primers to the db should go
        reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml"  # location where the blast record from comparing the reverse primers to the db should go
        cj0483_20_mid_bp_rm_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_20_middle_bp_removed.fasta"  # complete

        forward_blast = PCRPrediction.create_blastn_object(forward_primers, cj0483_20_mid_bp_rm_amp_seq,
                                                           forward_out_file, 'forward')
        reverse_blast = PCRPrediction.create_blastn_object(reverse_primers, cj0483_20_mid_bp_rm_amp_seq,
                                                           reverse_out_file, 'reverse')

        lo_forward_hsp_objects = PCRPrediction.create_hsp_objects(forward_blast)
        lo_reverse_hsp_objects = PCRPrediction.create_hsp_objects(reverse_blast)

        for f_hsp_object in lo_forward_hsp_objects:
            for r_hsp_object in lo_reverse_hsp_objects:
                if f_hsp_object.name == r_hsp_object.name:
                    if f_hsp_object.contig_name == r_hsp_object.contig_name:
                        # f_hsp_object and r_hsp_object should be on the same contig
                        self.assertEqual(f_hsp_object.contig_name, r_hsp_object.contig_name)
                        val = PCRPrediction.pcr_directly(f_hsp_object, r_hsp_object, cj0483_20_mid_bp_rm_amp_seq, amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects)
                        self.assertTrue(val)


#TODO: test to make sure the correct number of hsp's are created!!!


if __name__ == '__main__':
    unittest.main()



