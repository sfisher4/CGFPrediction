
class HSP(object):
    """A High Sequence Pair object from running blastn on two nucleotide sequences

    Attributes:
        name: The name of the HSP based off of the query name
        expect: The expected frequency of chance occurrence of a HSP
        start: The start element of the hsp subject
        end: The end element of the hsp subject
        query_start: The start element of the hsp query
        query_end: The end element of the hsp query
        strand: True if hsp located on leading strand (+, template) and false if hsp located on lagging strand (-, complement)
        length: length of the hsp
        valid: True if the the leading and lagging strand is opposite (hsp is opposite to the other hsp on the same gene)
        db_length: the length of the db that the high sequence pair is located on
        contig_name: The name of the HSP based off of its contig
        snp: True if there is a 3' SNP within SNP_THRESHOLD distance from the 3' end of the hsp object compared to the primer.
        sbjct: The subject's sequence
        query: The queries sequence
    """
    def __init__(self, name):
        self.name = name #The name of the hsp... including which node (contig) located on
        self.expect = -1
        self.start = -1 #not possible to have a start or end that starts at -1
        self.end = -1   #not possible to have a start or end that starts at -1
        self.query_start = -1
        self.query_end = -1
        self.strand = None #assert strand is not None
        self.length = 0
        self.valid = None
        self.db_length = 0
        self.contig_name = "" #ACTUALLY PUT THE CONTIG THE HSP IS ON!
        self.snp = None
        self.sbjct = ""
        self.query = ""



