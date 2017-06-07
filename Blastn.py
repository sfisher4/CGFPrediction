from Bio.Blast import NCBIXML
from Bio import SeqIO

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
HSP_THRESHOLD = 0.9

class Blastn(object):
    """ A blastn query from comparing two nucleotide sequences. Blastn have the following properites:

    Attributes:
        # name_query: A string representing the name of the query run on blastn
        # name_db: A string representing the name of the database run on blastn
        blast_records: A list of blast records created from running blastn
        hsp_records: A list of hsp records that contain the blast_records that meet threshold value length
    """

    def __init__(self):
        """Return a Blastn object with the following attributes from a blastn query:
        - blast records
        - hsp records """

        self.blast_records = []
        self.hsp_records = {}

    def create_blast_records(self, str):
            result_handle = open(str)
            blast_records = NCBIXML.parse(result_handle)  # returns an iterator
            self.blast_records = list(blast_records)
            result_handle.close()

    # genes is a list of primers when called by pcr directly and the db when called by entire_gene
    #TODO: merge with create_hsp_records_entire_gene
    def create_hsp_records(self, genes):
        dict_hsp = {}
        # lo_hsp = []
        full_lo_hsp = []
        # print('br', self.blast_records)
        for blast_record in self.blast_records:
            for gene in SeqIO.parse(genes, "fasta"):
                if blast_record.query in gene.name:
                    for alignment in blast_record.alignments:
                        lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / len(gene.seq))] #TODO: consider gene 1324 and 1134
                        if (len(lo_hsp) != 0):  # lo_hsp is not empty #TODO: correct indentation?
                            dict_hsp[alignment] = lo_hsp
                            full_lo_hsp.append(lo_hsp)
                            #break  # don't need to go through all gene seq'n in primers once the right alignment is found.
        self.hsp_records = dict_hsp

    # genes is a list of primers when called by pcr directly and the db when called by entire_gene
    def create_hsp_records_entire_gene(self, genes):
        dict_hsp = {}
        # lo_hsp = []
        full_lo_hsp = []
        # print('br', self.blast_records)
        for blast_record in self.blast_records:
            for gene in SeqIO.parse(genes, "fasta"):
                if blast_record.query in gene.name:
                    for alignment in blast_record.alignments:
                        # if 'NODE' in alignment.hit_def: #todo: change condition
                        lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / alignment.length)] #TODO: consider gene 1324 and 1134 #divide by the length of the contig
                        # else:
                        #     lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / len(gene.seq))]
                        if (len(lo_hsp) != 0):  # lo_hsp is not empty #TODO: correct indentation?
                            dict_hsp[alignment] = lo_hsp
                            full_lo_hsp.append(lo_hsp)
                            #break  # don't need to go through all gene seq'n in primers once the right alignment is found.
        self.hsp_records = dict_hsp
