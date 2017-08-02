from Bio.Blast import NCBIXML
from Bio import SeqIO
from HSP import HSP
import io

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold

#TODO: use Record.py instead?
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

    def create_blast_records(self, stdout_xml):
        """ Creates blast records for the Blastn object from xml file str.

        :param str: An xml file containing results from a Blastn query search.
        :return: None
        """

        with io.StringIO(stdout_xml) as blast_xml:
            blast_records = NCBIXML.parse(blast_xml)
            self.blast_records = list(blast_records)



    def create_hsp_objects(self, query_genes):
        """ Creates and initializes all fields of hsp objects from blastn_object input

        :param blastn_object: A Blastn Object
        :return: List of HSP Objects
        """
        hsp_objects = []
        dict_hsp = {}

        for blast_record in self.blast_records:
            for alignment in blast_record.alignments:
                for gene in SeqIO.parse(query_genes, "fasta"):
                    if blast_record.query in gene.name:
                        lo_hsp = [hsp for hsp in alignment.hsps]
                        dict_hsp[alignment] = lo_hsp

                for hsp in dict_hsp[alignment]:

                    # if hsp.expect < E_VALUE_CUTOFF: #removed this check b/c passed evalue into blastn query
                    hsp_name = blast_record.query
                    hsp_object = HSP(hsp_name)
                    hsp_object.start = hsp.sbjct_start
                    hsp_object.end = hsp.sbjct_end
                    hsp_object.query_start = hsp.query_start
                    hsp_object.query_end = hsp.query_end
                    hsp_object.alignment = alignment
                    hsp_object.contig_name = alignment.hit_def
                    hsp_object.length = abs(hsp_object.end - hsp_object.start) + 1
                    hsp_object.query_length = abs(hsp.query_end - hsp.query_start) + 1
                    hsp_object.db_length = alignment.length
                    hsp_object.expect = hsp.expect
                    hsp_object.sbjct = hsp.sbjct
                    hsp_object.query = hsp.query
                    hsp_object.identities = hsp.identities
                    hsp_object.gaps = hsp.gaps
                    hsp_object.bits = hsp.bits

                    # assuming no contigs (complete genome)
                    if hsp.sbjct_start < hsp.sbjct_end:
                        hsp_object.strand = True
                    elif hsp.sbjct_start > hsp.sbjct_end:
                        hsp_object.strand = False
                    assert hsp.sbjct_start != hsp.sbjct_end

                    hsp_objects.append(hsp_object)
        self.hsp_objects = hsp_objects
