
from memory_profiler import profile
# import attr

# @attr.s
class Results():
    """
    Attributes:
        both_primers_found: A boolean, True if both primers are found and false otherwise.
        contig: True if located on same contig as another hsp or false otherwise
        ehybrid: True if ehybrid passes and false otherwise
        epcr: True if ehybrid passes and false otherwise
    """
    #attrs... TODO
    # both_primers_found = attr.ib(default=None)
    # contig = attr.ib(default=None)
    # ehybrid = attr.ib(default=None)
    # epcr = attr.ib(default=None)
    # pcr_distance = attr.ib(default=None)
    # location = attr.ib(default=None)
    # valid = attr.ib(default=None)
    # strand = attr.ib(default=None)
    # snp = attr.ib(default=None)
    # end_dist = attr.ib(default=-1)
    # amp_len = attr.ib(default=-1)
    # amp_query = attr.ib(default="")
    # amp_sbjct = attr.ib(default="")
    # partner = attr.ib(default=None)

    def __init__(self, name):
        self.both_primers_found = None
        self.contig = None
        self.ehybrid = None
        self.epcr = None
        self.pcr_distance = None
        self.location = None
        self.valid = None
        self.strand = None
        self.snp = None #TODO: changed for testing
        self.end_dist = -1
        # self.amp_found = False
        self.amp_len = -1
        self.amp_query = ""
        self.amp_sbjct = ""
        self.partner = None
        self.snp_match = None #True is forward primer, False is reverse primer
    #     #TODO: where is the snp?


