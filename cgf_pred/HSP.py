
from cgf_pred.Results import Results


#@attr.s
class HSP(Results):
    """A High Sequence Pair object from running blastn on two nucleotide sequences

    Attributes:
        name: The name of the HSP based off of the query name
        expect: The expected frequency of chance occurrence of a HSP
        start: The start element of the hsp in the subject
        end: The end element of the hsp in the subject
        query_start: The start element of the hsp in the query
        query_end: The end element of hsp in the query
        identities:
        strand: True if hsp located on leading strand (+, template) and false if hsp located on lagging strand (-, complement)
        length: length of the hsp
        valid: True if the the leading and lagging strand is opposite (hsp is opposite to the other hsp on the same gene)
        db_length: the length of the db that the high sequence pair is located on
        contig_name: The name of the HSP based off of its contig
        snp: True if there is a 3' SNP within SNP_THRESHOLD distance from the 3' end of the hsp object compared to the primer.
        sbjct: The subject's sequence (db)
        query: The queries sequence
        gaps: Number of gaps in the hsp
    """
    # name = attr.ib()
    # expect = attr.ib(default= -1)
    # start = attr.ib(default= -1)
    # end = attr.ib(default= -1)
    # query_start = attr.ib(default= -1)
    # query_end = attr.ib(default= -1)
    # identities = attr.ib(default= -1)


    def __init__(self, name):
        super(HSP, self).__init__(name)
        self.name = name
        self.expect = -1
        self.start = -1 #not possible to have a start or end that starts at -1
        self.end = -1   #not possible to have a start or end that starts at -1
        self.query_start = -1
        self.query_end = -1
        self.identities = -1
        self.length = 0
        self.db_length = 0
        self.contig_name = ""
        self.sbjct = ""
        self.query = ""
        self.gaps = 0
        self.bsr = -1
        self.duplicate = False

    def set_name(self, name):
        self.name = name

    def __eq__(self, other):
        return self.sbjct == other.sbjct and self.contig_name == other.contig_name




