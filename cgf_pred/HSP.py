
from cgf_pred.Results import Results

class HSP(Results):
    """A High Sequence Pair object from running blastn on two nucleotide sequences

    Attributes:
        name: The name of the HSP based off of the query name
        query_start: The start element of the hsp in the query
        query_end: The end element of hsp in the query
        identities:
        start: The start element of the hsp in the subject
        end: The end element of the hsp in the subject
        length: length of the hsp
        db_length: the length of the db that the high sequence pair is located on
        contig_name: The name of the HSP based off of its contig
        sbjct: The subject's sequence (db)
        query: The queries sequence
        bsr: The Blast Score Ratio HSP
    """


    def __init__(self, name):
        super(HSP, self).__init__(name)
        if '11168_' in name:
            self.name = name
        else:
            self.name = '11168_' + name
        self.query_start = -1
        self.query_end = -1
        self.identities = -1
        self.start = -1
        self.end = -1
        self.length = 0
        self.db_length = 0
        self.contig_name = ""
        self.sbjct = ""
        self.query = ""
        self.bsr = -1

    def set_name(self, name):
        self.name = name


    def __eq__(self, other):
        return self.contig_name == other.contig_name and self.start == other.start and self.end == other.end




