from Bio.Blast import NCBIXML
from Bio import SeqIO

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
HSP_THRESHOLD = 0.9 #todo: determine a HSP_THRESHOLD

class Blastn(object):
    """ A blastn query Object from comparing two nucleotide sequences.

    Attributes:
        blast_records: A list of blast records created from running blastn
        hsp_records: A list of hsp records that contain the blast_records that meet threshold value length
    """

    def __init__(self):
        self.blast_records = []
        self.hsp_records = {}

    def create_blast_records(self, str):
        """ Creates blast records for the Blastn object from xml file str.

        :param str: An xml file containing results from a Blastn query search.
        :return: None
        """
        result_handle = open(str)
        blast_records = NCBIXML.parse(result_handle)  # returns an iterator
        self.blast_records = list(blast_records)
        result_handle.close()

    def create_hsp_records(self, query_genes):
        """ Creates hsp records for the Blastn object from fasta file genes.

        :param query_genes: A fasta file that contains query genes. Should be primers or amplicon sequences.
        :return: None
        """
        dict_hsp = {}
        full_lo_hsp = []
        for blast_record in self.blast_records:
            for gene in SeqIO.parse(query_genes, "fasta"):
                if blast_record.query in gene.name:
                    for alignment in blast_record.alignments:
                        lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / min(len(gene.seq), alignment.length))] #TODO: consider gene 1324 and 1134 #len(gene.seq) for primers and alignment.length for entire gene search
                        if (len(lo_hsp) != 0):  # lo_hsp is not empty
                            dict_hsp[alignment] = lo_hsp
                            full_lo_hsp.append(lo_hsp)
        self.hsp_records = dict_hsp
