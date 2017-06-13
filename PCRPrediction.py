from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Blastn import Blastn
from HSP import HSP

#TODO: determine all threshold values
E_VALUE_THRESHOLD = 0.04
MAX_DIST_BTW_CONTIGS = 600
MAX_DIST_BTWN_PRIMERS = 50 #if there is 50 missing, it will work
CUTOFF_GENE_LENGTH = 60
SNP_THRESHOLD = 5 #5 bp must be an exact match with 3' end of primer

def blastn_query(query_genes, database, out_file):
    """ Outputs a blastn query into an xml file.

    :param query_genes: A fasta file that contains the query genes that are being searched in the database
    :param database: A fasta file that contains the database that is being searched against
    :param out_file: A xml file that the blastn query
    :restrictions: Database is formatted using makeblastdb
    :return: None
    """
    #TODO: determine a word_size
    #TODO: restrict the evalue and mismatches
    blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=4, outfmt=5, out=out_file, evalue=E_VALUE_THRESHOLD)
    stdout, stderr = blastn_cline()

def create_blastn_object(query_genes, database, out_file):
    """ Return a blastn object with initialized blast_records and hsp_records

    :param query_genes: A fasta file that contains the query genes that are being searched in the database
    :param database: A fasta file that contains the database that is being searched against
    :param out_file: A xml file that the blastn query
    :restrictions: Database is formatted using makeblastdb
    :return: Blastn object
    """
    blastn_object = Blastn()
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_records(query_genes)
    return blastn_object

def create_hsp_objects(blastn_object):
    """ Returns a hsp object, initializing all fields

    :param blastn_object: A Blastn Object
    :return: HSP Object
    """
    hsp_objects = []

    for blast_record in blastn_object.blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if alignment in blastn_object.hsp_records:
                    if hsp.expect < E_VALUE_THRESHOLD: #TODO: maybe remove this check b/c added evalue in blast query!!!
                        hsp_name = blast_record.query
                        hsp_object = HSP(hsp_name)  #creates a HSP object
                        hsp_object.start = hsp.sbjct_start
                        hsp_object.end = hsp.sbjct_end
                        hsp_object.query_start = hsp.query_start
                        hsp_object.query_end = hsp.query_end
                        hsp_object.contig_name = alignment.hit_def

                        # assuming no contigs (complete genome)
                        if hsp.sbjct_start < hsp.sbjct_end:
                            hsp_object.strand = True
                        elif hsp.sbjct_start > hsp.sbjct_end:
                            hsp_object.strand = False
                        assert hsp.sbjct_start != hsp.sbjct_end

                        hsp_object.length = abs(hsp_object.end - hsp_object.start) + 1
                        hsp_object.db_length = alignment.length
                        hsp_object.expect = hsp.expect
                        hsp_object.sbjct = hsp.sbjct
                        hsp_object.query = hsp.query

                        hsp_objects.append(hsp_object)
    return hsp_objects

def valid_strands(first_hsp_object, second_hsp_object):
    """ Assigns valid attributes to first and second hsp object and modifies lo_first and lo_second hsp objects s.t. they only contain strands with the correct orientation to one another.

    :param first_hsp_object: A HSP object to compare with second_hsp_object
    :param second_hsp_object: A HSP object to compare with first_hsp_object
    :param lo_forward_hsp_objects: A list of hsp objects that was run using forward primer sets
    :param lo_reverse_hsp_objects: A list of hsp objects that was run using reverse primer sets
    :return: None
    """

    if first_hsp_object.name == second_hsp_object.name: #forward and reverse hsp's are from the same query (Note: could still be on different contigs)
        if (first_hsp_object.strand or second_hsp_object.strand) and not (first_hsp_object.strand and second_hsp_object.strand):
            first_hsp_object.valid = True
            second_hsp_object.valid = True
        else:
            first_hsp_object.valid = False
            second_hsp_object.valid = False
            # lo_forward_hsp_objects.remove(first_hsp_object)
            # lo_reverse_hsp_objects.remove(second_hsp_object)

# TODO: faster way?
def create_primer_dict(primers):
    """ Return a dictionary that contains primer sequences: where key is the primer id and value is the primer sequence.

    :param primers: A fasta file containing primer sequences
    :return: A dictionary
    """
    primer_dict = {}
    for primer in SeqIO.parse(primers, "fasta"):
        primer_dict[primer.id] = primer.seq
    return primer_dict

#TODO: only considers if the primers are on the same contig!
#TODO: If this is the way that the program should function, then change this function s.t. the hsp objects have a forward and reverse primer attribute
#TODO: and... make less naive!
#TODO: make sure hsp_object.name is correct!!!
def is_snp(hsp_object, f_primers, r_primers):
    """Return True if a mismatch or removed bp is located in the first SNP_THRESHOLD bp's on the 3' end

    :param hsp_object: A HSP Object.
    :param f_primers: A dictionary containing full forward primer sequences.
    :param r_primers: A dictionary containing full reverse primer sequences.
    :return: A boolean: True if mismatch or removed bp on 3' end, False otherwise
    """
    #TODO: test that this will not work if the primer sequences do not correspond to the hsp_object
    if f_primers[hsp_object.name] and r_primers[hsp_object.name] != "":
        f_primer_end = f_primers[hsp_object.name]
        r_primer_end = r_primers[hsp_object.name]
    # print('primer', f_primer_end[len(f_primer_end) - SNP_THRESHOLD: len(f_primer_end)])
    # print('query', hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD: len(hsp_object.sbjct)])
    # print(f_primer_end)
        if str(f_primer_end[len(f_primer_end) - SNP_THRESHOLD: len(f_primer_end)]) in hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD: len(hsp_object.sbjct)]:
            hsp_object.snp = False
        elif str(r_primer_end[len(r_primer_end) - SNP_THRESHOLD: len(r_primer_end)]) in hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD: len(hsp_object.sbjct)]:
            hsp_object.snp = False
        else:
            hsp_object.snp = True

def is_distance(f_hsp_object, r_hsp_object, amplicon_sequences):
    """ Return True if f_hsp_object and r_hsp_object are within MAX_DIST_BTWN_PRIMERS to each other.

    :param f_hsp_object: A HSP Object on a contig
    :param r_hsp_object: A HSP Object on a contig
    :param amplicon_sequences: A Fasta file that contains amplicon sequences for comparison purposes.
    :restrictions: f_hsp_object and r_hsp_object must be on the same contig.
    :return: A Boolean
    """
    assert f_hsp_object.contig_name == r_hsp_object.contig_name    #located on same contig
    q_distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        if f_hsp_object.name in amplicon.id:
            amplicon_length = len(amplicon.seq)
            return (abs(amplicon_length - q_distance) <= MAX_DIST_BTWN_PRIMERS)

def pcr_directly(forward_hsp_object, reverse_hsp_object, amplicon_sequences):
    """ Determines if forward and reverse hsp object are in the correct orientation and distance from each other to be considered a hit.

    :param forward_hsp_object: A HSP object in question found from lo_forward_hsp_objects
    :param reverse_hsp_object: A HSP object in question found from lo_reverse_hsp_object
    :param amplicon_sequences: A Fasta file that contains amplicon sequences for comparison purposes.
    :return: A Boolean
    """

    valid_strands(forward_hsp_object, reverse_hsp_object)

    #TODO: Possibly add snp condition or return a flag if there is a snp on the 3' end

    assert forward_hsp_object.valid == reverse_hsp_object.valid #should always be True because valid_strands assigns both objects the same valid attribute
    if forward_hsp_object.valid == True: #and forward_hsp_object.snp == False: #foward and reverse hsp objects have the same valid value
        return is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences)
    elif forward_hsp_object.valid == False:
        return False
    # else:
        #!!! forward hsp object is not initialized!!!

# #determines if two strands are facing each other when comparing contigs to amp sequences
# #TODO: !!!
# def is_valid_entire_strands(hsp_1, hsp_2):
#     """ Determine if hsp_1 and hsp_2 are facing each other
#
#     :param hsp_1:
#     :param hsp_2:
#     :return:
#     """
#     if hsp_1.strand == hsp_2.strand:
#         hsp_1.valid = True
#         hsp_2.valid = True
#         return True
#     else:
#         return False

def valid_dir(lo_hsps):
    """ Modify lo_hsps to only contain primer sequences that are facing the end of the contig.

    :param lo_hsps: list of hsp objects
    :return: None
    """
    for hsp in lo_hsps:
        # considers the case where the entire strand is found
        if (hsp.end == hsp.query_end or hsp.start == hsp.query_start) and not (hsp.end == hsp.query_end and hsp.start == hsp.query_start):
            if hsp.strand == False and hsp.end == hsp.query_end:  # TODO: create threshold value for end differences
                print(hsp.end)
                print(hsp.query_end)
                print(hsp.strand)
                lo_hsps.remove(hsp)
            elif hsp.strand == True and hsp.start == hsp.query_start:
                print(hsp.start)
                print(hsp.query_start)
                print(hsp.strand)
                lo_hsps.remove(hsp)

# Given a list of hsp objects and a reference object, finds all of the occurrences of the reference object (name)
# and determines if it is of long enough length and in proper orientation to be considered found.
# Assumes that there are only 2 hsp objects of reference type in lo hsp objects
#NOTE: f_object.name == r_object.name
#only considers the case where one or two primers are found (not > 3)
def entire_gene(lo_hsp_objects, reference_object):
    """ Return the hsp objects from the same blast record query as reference_object that are of proper orientation and long enough length to be considered found.

    :param lo_hsp_objects: A list of hsp objects from a blastn search.
    :param reference_object: A HSP object
    :return: List of hsp objects
    """
    #TODO: consider the case where lo_hsps have 3 elements?
    #It would be highly unlikely that a gene will be found of long enough length on a database more than once?

    lo_hsps = []
    for hsp_object in lo_hsp_objects:
        count = 0
        if reference_object.name in hsp_object.name:  # hsp_object.name may be longer than reference_object.name
            if not hsp_object.length <= CUTOFF_GENE_LENGTH:
                count += 1
                assert count != 3 #Deal with this case later if needed
                lo_hsps.append(hsp_object)
        valid_dir(lo_hsps) #make sure this modifies the list !!!
        if len(lo_hsps) == 1:
            print("only one primer found for", lo_hsps[0].name, "on", lo_hsps[0].contig_name)
    return lo_hsps

#TODO: Change return value to hash?
def pcr_prediction(forward_primers, reverse_primers, database, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file):
    """Return a list of hsp's that would be predicted as hsp's in vitro...

    :param forward_primers: A Fasta file containing forward query primers (ids and sequences)
    :param reverse_primers: A Fasta file containing reverse query primers (ids and sequences)
    :param database: A Fasta file (formatted using makeblastdb) that contains a database (complete or draft)
    :param forward_out_file: An xml file location to store blastn results from forward primers as queries
    :param reverse_out_file: An xml file location to store blastn results from reverse primers as queries
    :param amplicon_sequences: A Fasta file containing amplicon sequences
    :param full_out_file: An xml file location to store blastn results from amplicon sequences as queries
    :return: List of HSP's that were found in database.
    """

    forward_blast = create_blastn_object(forward_primers, database, forward_out_file)  # using amplicon_seq for testing purposes #expect 40 blast records and 1 hsp record
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file)  # using amplicon_seq for testing purposes #expect 40 blast records and 1 hsp record
    blast_object = create_blastn_object(amplicon_sequences, database, full_out_file) #expect 1 blast record and 2 hsp records

    lo_forward_hsp_objects = create_hsp_objects(forward_blast)
    lo_reverse_hsp_objects = create_hsp_objects(reverse_blast)
    lo_hsp_objects = create_hsp_objects(blast_object)

    lo_queries = []

    #TODO: make less naive
    for f_hsp_object in lo_forward_hsp_objects:
        for r_hsp_object in lo_reverse_hsp_objects:
            if f_hsp_object.name == r_hsp_object.name: #if the f_hsp and r_hsp are from the same primer query
                if f_hsp_object.contig_name == r_hsp_object.contig_name: #f_hsp and r_hsp are on the same contig
                    if pcr_directly(f_hsp_object, r_hsp_object, amplicon_sequences):
                        lo_f_r = [f_hsp_object, r_hsp_object]
                        lo_queries.append(lo_f_r)
                else: # primer located on different contigs
                    lo_tmp = entire_gene(lo_hsp_objects, f_hsp_object)
                    if len(lo_tmp) > 0:
                        lo_queries.append(lo_tmp)
    print('lo queries', lo_queries)
    return(lo_queries)

def main():

    database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta" #contains id and a complete genome sequence
    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #fasta file with primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #fasta file with primer id's and primer sequences
    forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml" #location where the blast record from comparing the forward primers to the db should go
    reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml" #location where the blast record from comparing the reverse primers to the db should go
    full_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"

    #Test amp seq
    cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta" #complete
    cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
    cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta" #complete
    cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
    cj1134_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj1134.fasta"
    cj1324_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj1324.fasta"

   # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

    #test data
    test_contig_51_bp_removed = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_51_middle_bp_removed.fasta"
    test_mystery_db_name_same_contig = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db.fasta"
    test_one_contig = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_primer_same_contig_with_2_contigs.fasta"
    test_mystery_db_name_entire_gene = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db_entire_gene_61_bp_each_contig.fasta"
    cj0008_rm_2_ws = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed/cj0008_rm_2_ws.fasta"

    # using an amplicon_sequence as database for testing purposes
    #takes in forward_primers and reverse_primers (queries), database (contigs or complete) [only one database!!!], full amplicon sequences that contain the forward and reverse primers inputted,
    # a forward and reverse out file where the blast results can go, a complete amplicon
    return(pcr_prediction(forward_primers, reverse_primers, cj0483_complete_amp_seq, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file))


if __name__ == "__main__": main()

#TODO: With database input, makeblastdb???
#TODO: Put functions into HSP that might be useful (thought for later)
#TODO: Is createObjects in the right spot?
#TODO: Change main to input files (currently embedded in main)? or somehow create a way to input files?