
class HSP(object):

    """A High Sequence Pair object found from running Blastn that has the following attributes:
        name: The name of the hsp (ran on blastn) based off of the query name
        start: The start element of the hsp subject
        end: The end element of the hsp subject
        strand: True if hsp located on leading strand (+, template) and false if hsp located on lagging strand (-, complement)
        length: length of the hsp
        valid: True if the the leading and lagging strand is opposite (hsp is opposite to the other hsp on the same gene)
        db_length: the length of the db that the high sequence pair is located on 
        blast_record_title: Title of the blast record that the hsp was retrieved from.
        snp: True if there is a 3' SNP within SNP_THRESHOLD distance from the 3' end of the hsp object compared to the primer.
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



