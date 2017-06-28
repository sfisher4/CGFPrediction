import unittest
import PCRPrediction


class MyTestCase(unittest.TestCase):

    def test_pcr_prediction(self):
        forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
        reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
        amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
        test_11168 = "/home/sfisher/Sequences/11168_complete_genome"
        result_dict = PCRPrediction.main(test_11168, forward_primers, reverse_primers, amplicon_sequences)
        self.assertEqual(len(result_dict), 1)
        result = result_dict["11168_complete_genome.fasta"]
        self.assertEqual(len(result), 40)
        for lo_hsp in result:
            self.assertEqual(len(lo_hsp), 2)
            self.assertEqual(lo_hsp[0].name, lo_hsp[1].name)


if __name__ == '__main__':
    unittest.main()
