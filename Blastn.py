from Bio.Blast import NCBIXML
from Bio import SeqIO
from HSP import HSP

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
HSP_THRESHOLD = 0.9 #todo: determine a HSP_THRESHOLD

class Blastn(object):
    """ A blastn query Object from comparing two nucleotide sequences.

    Attributes:
        blast_records: A list of blast records created from running blastn
        hsp_records: A list of hsp records that contain the blast_records that meet threshold value length
        hsp_objects: A list of hsp objects that contain the blast_records that meet threshold value length
    """

    def __init__(self):
        self.blast_records = []
        #changed!!!
        # self.hsp_records = {}
        self.hsp_objects = []

    def create_blast_records(self, str):
        """ Creates blast records for the Blastn object from xml file str.

        :param str: An xml file containing results from a Blastn query search.
        :return: None
        """
        result_handle = open(str)
        blast_records = NCBIXML.parse(result_handle)  # returns an iterator
        self.blast_records = list(blast_records)
        result_handle.close()

#changed!!!
    # def create_hsp_records(self, query_genes):
    #     """ Creates hsp records for the Blastn object from fasta file genes.
    #
    #     :param query_genes: A fasta file that contains query genes. Should be primers or amplicon sequences.
    #     :return: None
    #     """
    #     dict_hsp = {}
    #
    #     for blast_record in self.blast_records:
    #         for alignment in blast_record.alignments:
    #             for gene in SeqIO.parse(query_genes, "fasta"):
    #                 if blast_record.query in gene.name:
    #                     lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / min(len(gene.seq), alignment.length))] #len(gene.seq) for primers and alignment.length for entire gene search
    #                     if (len(lo_hsp) != 0):  # lo_hsp is not empty
    #                         dict_hsp[alignment] = lo_hsp
    #     for x in dict_hsp:
    #         print('x', x)
    #     print('dict hsp', dict_hsp)
    #     self.hsp_records = dict_hsp

    def create_hsp_objects(self, query_genes):
        """ Creates and initializes all fields of hsp objects from blastn_object input

        :param blastn_object: A Blastn Object
        :return: List of HSP Objects
        """
        #TODO: make a dict?
        hsp_objects = []
        dict_hsp = {}

        #TODO: speed up
        for blast_record in self.blast_records:
            for alignment in blast_record.alignments:
                for gene in SeqIO.parse(query_genes, "fasta"):
                    #TODO: merge if statement with list comprehension
                    if blast_record.query in gene.name:
                        lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / min(len(gene.seq), alignment.length))]
                        dict_hsp[alignment] = lo_hsp

                for hsp in dict_hsp[alignment]:
                    # if hsp.expect < E_VALUE_THRESHOLD: #removed this check b/c passed evalue into blastn query
                    hsp_name = blast_record.query
                    hsp_object = HSP(hsp_name)
                    hsp_object.start = hsp.sbjct_start
                    hsp_object.end = hsp.sbjct_end
                    hsp_object.query_start = hsp.query_start
                    hsp_object.query_end = hsp.query_end
                    hsp_object.alignment = alignment
                    hsp_object.contig_name = alignment.hit_def
                    hsp_object.length = abs(hsp_object.end - hsp_object.start) + 1
                    hsp_object.db_length = alignment.length
                    hsp_object.expect = hsp.expect
                    hsp_object.sbjct = hsp.sbjct
                    hsp_object.query = hsp.query
                    hsp_object.identities = hsp.identities

                    # assuming no contigs (complete genome)
                    if hsp.sbjct_start < hsp.sbjct_end:
                        hsp_object.strand = True
                    elif hsp.sbjct_start > hsp.sbjct_end:
                        hsp_object.strand = False
                    assert hsp.sbjct_start != hsp.sbjct_end

                    hsp_objects.append(hsp_object)
        self.hsp_objects = hsp_objects
