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

    def __init__(self, name):
        """Return a Blastn object whose query name is , results from blastn is 
        blast_records and hsp results is hsp_records"""

        # self.name_query = name_query
        self.name = name
        self.blast_records = []
        self.hsp_records = {}

    def create_blast_records(self, str):
            result_handle = open(str)
            blast_records = NCBIXML.parse(result_handle)  # returns an iterator
            self.blast_records = list(blast_records)
            result_handle.close()

#TODO: combine create_hsp_records and create_hsp_records_full_gene_search into one function
    def create_hsp_records(self, genes): #genes is a list of primers when called by pcrdirectly and the db when called by entire_gene
        dict_hsp = {}
        lo_hsp = []

        for blast_record in self.blast_records:
            for alignment in blast_record.alignments:
                for gene in SeqIO.parse(genes, "fasta"):
                    if alignment.title.find(gene.name) != -1:   #ensures the hsp alignment is on the right gene
                                # print('name', self.name)
                                # print(hsp.identities / len(gene.seq))
                                # print('hsp identities', hsp.identities)
                                # print('gene length', len(gene.seq))
                        lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / len(gene.seq))] #TODO: consider gene 1324 and 1134
                    if (len(lo_hsp) != 0):  # lo_hsp is not empty #TODO: correct indentation?
                        # print('ratio', hsp.identities / len(gene.seq))
                        # print('title', alignment.title)
                        dict_hsp[alignment] = lo_hsp
                        #break  # don't need to go through all gene seq'n in primers once the right alignment is found.
        self.hsp_records = dict_hsp





