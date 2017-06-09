from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Blastn import Blastn
from HSP import HSP

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
MAX_DIST_BTW_CONTIGS = 600
MAX_DIST_BTWN_PRIMERS = 50 #if there is 50 missing, it will work
CUTOFF_GENE_LENGTH = 60
SNP_THRESHOLD = 5 #5 bp must be an exact match with 3' end of primer

#create a blastn query and output into an xml file
def blastn_query(query_genes, database, out_file):
    #TODO: Formatted using makeblastdb : query_genes, database (fasta), out_file (xml)

    #TODO: determine a word_size
    #query the database using blastn
    blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=4, outfmt=5, out=out_file) #TODO: restrict the evalue and mismatches
    stdout, stderr = blastn_cline()

# returns a Blastn object with blast_records and hsp_records initialized,
# given two nucl seq to compare (query and db), and outfile
# Note: Currently only works with primers as query_genes
def create_blastn_object(query_genes, database, out_file):
    blastn_object = Blastn()
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_records(query_genes) #todo: changed for testing with full gene
    return blastn_object

# returns a Blastn object with blast_records and hsp_records initialized,
# Note: currently only works with amplicon sequence as query_genes
#TODO: merge with create_blastn_object function!!!
def create_blastn_full_search_object(query_genes, database, out_file):
    blastn_object = Blastn()
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_records_entire_gene(query_genes) #todo: change me!!!
    return blastn_object

#Creates a hsp object, initializing all fields.
def create_hsp_objects(blastn_object):

    hsp_objects = []

    for blast_record in blastn_object.blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if alignment in blastn_object.hsp_records:
                    if hsp.expect < E_VALUE_THRESHOLD:
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

# Given a forward and reverse hsp object, checks if they are in the correct strand orientation to one another
# (ie. one leading (True) and one lagging(False)).
# and assigns their strand attributes, where True is the correct orientation and False is not.
# If False, removes the strands from the forward and reverse object lists.
#TODO: remove lo_forward_hsp_objects??
#TODO: changed to only remove from lo_hsp_objects
def valid_strands(first_hsp_object, second_hsp_object, lo_forward_hsp_objects, lo_reverse_hsp_objects):

    if first_hsp_object.name == second_hsp_object.name: #forward and reverse hsp's are from the same query (Note: could still be on different contigs)
        if (first_hsp_object.strand or second_hsp_object.strand) and not (first_hsp_object.strand and second_hsp_object.strand):
            first_hsp_object.valid = True
            second_hsp_object.valid = True
        else:
            first_hsp_object.valid = False
            second_hsp_object.valid = False
            lo_forward_hsp_objects.remove(first_hsp_object)
            lo_reverse_hsp_objects.remove(second_hsp_object)

# todo: faster way?
def create_primer_dict(primers):
    primer_dict = {}
    for primer in SeqIO.parse(primers, "fasta"):
        primer_dict[primer.id] = primer.seq
        # print('pid', primer.id)
    return primer_dict

#Returns True if there is a mismatch in the first SNP_THRESHOLD base pairs and False otherwise
#f_primers and r_primers are dictionaries containing the primer sequences that was used as the query for blastn search.
#TODO: only considers if the primers are on the same contig!
#TODO: If this is the way that the program should function, then change this function s.t. the hsp objects have a forward and reverse primer attribute
#TODO: and... make less naive!
#TODO: make sure hsp_object.name is correct!!!
def is_snp(hsp_object, f_primers, r_primers):
    f_primer_end = f_primers[hsp_object.name]
    r_primer_end = r_primers[hsp_object.name]
    # print('primer', f_primer_end[len(f_primer_end) - SNP_THRESHOLD: len(f_primer_end)])
    # print('query', hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD: len(hsp_object.sbjct)])
    # print(f_primer_end)
    if str(f_primer_end[len(f_primer_end) - SNP_THRESHOLD: len(f_primer_end)]) in hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD: len(hsp_object.sbjct)]:
        hsp_object.snp = False #Not 3' SNP
    elif str(r_primer_end[len(r_primer_end) - SNP_THRESHOLD: len(r_primer_end)]) in hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD: len(hsp_object.sbjct)]:
        hsp_object.snp = False
    else:
        hsp_object.snp = True #3' SNP

# Returns true if forward_hsp and reverse_hsp have correct distance from each other, compared to the amplicon_sequences (includes primer strand length)
# and False otherwise.
# Requires: f and r hsp objects are on the same contig
def is_distance(f_hsp_object, r_hsp_object, amplicon_sequences, database):
    assert f_hsp_object.contig_name == r_hsp_object.contig_name    #located on same contig
    q_distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
    return is_correct_distance(amplicon_sequences, f_hsp_object, q_distance)

#returns true if the distance between the primers is within MAX_DIST_BTWN_PRIMERS, return false otherwise
#TODO: might want to change MAX_DIST_BTWN_PRIMERS to consider % of distance that is missing
def is_correct_distance(amplicon_sequences, f_hsp_object, q_distance):
    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        #naively assuming starts at 6
        #TODO: change to 'in'?
        amplicon_id = amplicon.id[6:]
        if f_hsp_object.name == amplicon_id:
            amplicon_length = len(amplicon.seq)
            return (abs(amplicon_length - q_distance) <= MAX_DIST_BTWN_PRIMERS) #returns bool

#Given a forward hsp object and a reverse hsp object, checks if they are valid strands to one another and
# if so, returns True if the distance between them are correct, compared to the amplicon sequence and False otherwise
#TODO: possibly change this fcn... too similar to find_distances.
def pcr_directly(forward_hsp_object, reverse_hsp_object, database, amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects):
    valid_strands(forward_hsp_object, reverse_hsp_object, lo_forward_hsp_objects, lo_reverse_hsp_objects)
    #TODO: possibly remove this condition... and return a flag if there is a snp on the 3' end
    if forward_hsp_object.valid == True: #and forward_hsp_object.snp == False: #foward and reverse hsp objects have the same valid value
        return is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences, database)
    else:
        return False

def is_hsp_length(hsp_object):
    if hsp_object.length <= CUTOFF_GENE_LENGTH:
        return False
    else:
        return True

#determines if two strands are facing each other when comparing contigs to amp sequences
def is_valid_entire_strands(hsp_1, hsp_2):
    if hsp_1.strand == hsp_2.strand:
        hsp_1.valid = True
        hsp_2.valid = True
        return True
    else:
        return False

#Given a list of hsp objects and a reference object, finds all of the occurences of the reference object (name)
# and determines if it of long enough length and in proper orientation to be considered found.
# Assumes that there are only 2 hsp objects of reference type in lo hsp objects
def entire_gene(lo_hsp_objects, reference_object):
    lo_queries = []

    for hsp_object in lo_hsp_objects:
        count = 0
        if reference_object.name in hsp_object.name:  # not checking == b/c forward_hsp_object can have more in its name than hsp_object
            #TODO: perhaps consider the case where lo_queries has 3 elements
            count += 1
            assert count != 3
            if is_hsp_length(hsp_object):
                lo_queries.append(hsp_object)
        if len(lo_queries) == 2:
            if not is_valid_entire_strands(lo_queries[0], lo_queries[1]):
                lo_queries = [] #Only contains 2 elements: remove both
    return lo_queries

    #checking all hsp objects (rather than just reference object)
    # lo_queries = entire_gene(lo_hsp_objects, f_hsp_object)
    # if len(lo_queries) == 2:  # this will be the case most of the time
    #     if valid_strands(lo_queries[0],lo_queries[1], lo_hsp_objects):
    #         l_tmp = [lo_queries[0], lo_queries[1]]

    # lo_hsp_objects.sort(key=lambda x: x.name)
    # if len(lo_queries) == 2: #this will be the case most of the time
    #     return valid_strands(lo_queries[0],lo_queries[1], lo_hsp_objects)
    # else:
    #     return False

    # while len(lo_hsp_objects) > 0:
    #     if len(lo_hsp_objects) >= 2:
    #         if lo_hsp_objects[0].name == lo_hsp_objects[1].name:
    #             lo_f_r = [lo_hsp_objects[0], lo_hsp_objects[1]]
    #             lo_queries.append(lo_f_r)
    #             lo_hsp_objects.pop(0)
    #             lo_hsp_objects.pop(0)
    #     else:
    #         lo_hsp_objects.pop(0)
    #return lo_queries

#Returns a list of hsp's that would be predicted as hsp's in vitro...
#TODO: Change return value to hash?
def pcr_prediction(forward_primers, reverse_primers, database, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file):

    forward_blast = create_blastn_object(forward_primers, database, forward_out_file)  # using amplicon_seq for testing purposes #expect 40 blast records and 1 hsp record
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file)  # using amplicon_seq for testing purposes #expect 40 blast records and 1 hsp record
    blast_object = create_blastn_full_search_object(amplicon_sequences, database, full_out_file) #expect 1 blast record and 2 hsp records

    lo_forward_hsp_objects = create_hsp_objects(forward_blast)
    lo_reverse_hsp_objects = create_hsp_objects(reverse_blast)
    lo_hsp_objects = create_hsp_objects(blast_object)

    lo_queries = []

    #TODO: make less naive
    for f_hsp_object in lo_forward_hsp_objects:
        for r_hsp_object in lo_reverse_hsp_objects:
            if f_hsp_object.name == r_hsp_object.name: #if the f_hsp and r_hsp are from the same primer query
                if f_hsp_object.contig_name == r_hsp_object.contig_name: #f_hsp and r_hsp are on the same contig
                    if pcr_directly(f_hsp_object, r_hsp_object, database, amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects):
                        lo_f_r = [f_hsp_object, r_hsp_object]
                        lo_queries.append(lo_f_r)
                else: # primer located on different contigs
                    lo_tmp = entire_gene(lo_hsp_objects, f_hsp_object)
                    if len(lo_tmp) == 2: #don't add empty list (Length 2 from implementation of entire_gene)
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
    return(pcr_prediction(forward_primers, reverse_primers, cj0483_contig_trunc_amp_seq, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file))


if __name__ == "__main__": main()

#TODO: Is it worth having Blastn as a class? (thought for later)
#TODO: Put functions into HSP that might be useful (thought for later)
#TODO: Is createObjects in the right spot?
#TODO: Change main to input files (currently embedded in main)? or somehow create a way to input files?
#TODO: change comments to be in python format (docstring)