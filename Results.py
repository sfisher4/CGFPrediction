

class Results:
    """
    Attributes:
        both_primers_found: A boolean, True if both primers are found and false otherwise.
        contig: True if located on same contig as another hsp or false otherwise
        ehybrid: True if ehybrid passes and false otherwise
        epcr: True if ehybrid passes and false otherwise
    """
    def __init__(self, name):
        self.both_primers_found = None
        self.contig = None
        self.ehybrid = None
        self.epcr = None
        self.pcr_distance = None
        self.location = None
        self.valid = None
        self.strand = None
        self.snp = None
        self.amp_len = 0
        #TODO: where is the snp?

