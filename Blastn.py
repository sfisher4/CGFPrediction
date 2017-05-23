class Blastn(object):
    """ A blastn query from comparing two nucleotide sequences. Blastn have the following properites:

    Attributes:
        name_query: A string representing the name of the query run on blastn
        name_db: A string representing the name of the database run on blastn
        blast_records: A list of blast records created from running blastn
        hsp_records: A list of hsp records that contain the blast_records that meet threshold value length
    """

    def __init__(self, blast_records, hsp_records):
        """Return a Blastn object whose query name is *query_name*, db name is *db_name* , results from blastn is 
        blast_records and hsp results is hsp_records"""

        # self.name_query = name_query
        # self.name_db = name_db
        self.hsp_records = hsp_records
        self.blast_records = blast_records
