
E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
HSP_THRESHOLD = 1

class Query(object):

    """A Query object that has the following attributes:
        name: The name of the query sequence (ran on blastn)
        start: The start element of the query in the db
        end: The end element of the query in the db
        strand: True if primer located on leading strand (+, template) and false if primer located on lagging strand (-, complement)
    """

    def __init__(self, name):
        self.name = name
        self.start = -1 #not possible to have a start or end that starts at -1
        self.end = -1   #not possible to have a start or end that starts at -1
        self.strand = None #assert strand is not None
        self.length = 0
        self.valid = None
        self.db_length = 0
        self.blast_record_title = ""



