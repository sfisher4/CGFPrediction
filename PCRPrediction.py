from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Blastn import Blastn
from HSP import HSP

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
MAX_DIST_BTW_CONTIGS = 600
MAX_DIST_BTWN_PRIMERS = 50

#create a blastn query and output into an xml file
def blastn_query(query_genes, database, out_file):
    #TODO: Formatted using makeblastdb : query_genes, database (fasta), out_file (xml)

    #TODO: determine a word_size
    #query the database using blastn
    blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=4, outfmt=5, out=out_file)
    stdout, stderr = blastn_cline()

# returns a Blastn object with it's name, blast_records and hsp_records initialized,
# given two nucl seq to compare (query and db), outfile and a name to identify the comparison made
# Note: Currently only works with primers as query_genes
def create_blastn_object(query_genes, database, out_file, name):
    blastn_object = Blastn(name)
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_records(query_genes) #todo: changed for testing with full gene

    return blastn_object

# returns a Blastn object with it's name, blast_records and hsp_records initialized,
# Note: currently only works with amplicon sequence as query_genes
def create_blastn_full_search_object(query_genes, database, out_file, name):
    blastn_object = Blastn(name)
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_records(database) #todo: changed for testing with full gene
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
                        # print('hit_def', alignment.hit_def)
                        hsp_object.contig_name = alignment.hit_def #TODO: from here!!!
                        # print(hsp.contig_name)

                        # assuming no contigs (complete genome)
                        if hsp.sbjct_start < hsp.sbjct_end:
                            hsp_object.strand = True
                        elif hsp.sbjct_start > hsp.sbjct_end:
                            hsp_object.strand = False
                        assert hsp.sbjct_start != hsp.sbjct_end

                        hsp_object.length = abs(hsp_object.end - hsp_object.start) + 1
                        hsp_object.db_length = alignment.length
                        hsp_object.expect = hsp.expect

                        hsp_objects.append(hsp_object)
    return hsp_objects

# Given a list of forward and reverse hsp_objects, modifies the lists to include only hsp_objects that are in the correct
# strand orientation (relative to one another) and assigns the objects valid attributes TODO: change my description!!!
#TODO: remove lo_forward_hsp_objects??
# (ie. one leading (True) and one lagging(False)).
def valid_strands(f_hsp_object, r_hsp_object, lo_forward_hsp_objects, lo_reverse_hsp_objects):

    if f_hsp_object.name == r_hsp_object.name: #forward and reverse hsp's are from the same query (Note: could still be on different contigs)
        if (f_hsp_object.strand or r_hsp_object.strand) and not (f_hsp_object.strand and r_hsp_object.strand):
            f_hsp_object.valid = True
            r_hsp_object.valid = True
        else:
            f_hsp_object.valid = False
            r_hsp_object.valid = False
            lo_forward_hsp_objects.remove(f_hsp_object)
            lo_reverse_hsp_objects.remove(r_hsp_object)

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

#returns a list of valid hsp's that have the correct distance from each other, compared to the amplicon_sequences (includes primer strand length)
#Given the hsp_objects (primer hits) that are located on the same contig.
#TODO: possibly change this fcn... too similar to find_distances.
def pcr_directly(forward_hsp_object, reverse_hsp_object, database, amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects):
    valid_strands(forward_hsp_object, reverse_hsp_object, lo_forward_hsp_objects, lo_reverse_hsp_objects)
    if forward_hsp_object.valid == True: #foward and reverse hsp objects have the same valid value
        return is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences, database)
    else:
        return False

#Returns a list of hsp_objects that are within MAX_DIST_BTW_CONTIGS from each other and located in lo_forward (and reverse_ hsp_objects
# Note: given a blast_object that contains hsp' that are on different contigs.
def entire_gene(amplicon_sequences, blast_object, lo_hsp_objects, lo_forward_hsp_objects, lo_reverse_hsp_objects):
    #TODO: check to make sure you are only checking the sequences that starts and end at the primer sequence you were looking for!!!
    #TODO: this is assuming we only found 2 hsp's that are on different contigs!!!
    print('oh no')
    lo_correct_distance_hsp_objects = []

    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        total_len_db = 0
        hsp_object_list = []
        length_amplicon = len(amplicon.seq)
        for f_hsp_object in lo_forward_hsp_objects:
            for hsp_object in lo_hsp_objects:
                if f_hsp_object.name in hsp_object.name: #ensures the primers are facing each other!!! #not checking == b/c forward_hsp_object can have more in its name than hsp_object
                        if hsp_object.name in amplicon.id: #purposely not checking == #TODO: the amplicon is part of the hsp found by the blast record (check name comparison is correct)
                            total_len_db += hsp_object.length
                            hsp_object_list.append(hsp_object)

        if abs(length_amplicon - total_len_db) <= MAX_DIST_BTW_CONTIGS:
            if len(hsp_object_list) != 0:
                lo_correct_distance_hsp_objects.append(hsp_object_list)
                # print(hsp_object_list)
                # print('len', total_len_db)

    print("lo correct distance queries", lo_correct_distance_hsp_objects)
    return lo_correct_distance_hsp_objects

#Returns a list of hsp's that would be predicted as hsp's in vitro...
#TODO: Change return value. (hash?)
def pcr_prediction(forward_primers, reverse_primers, database, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file):

    forward_blast = create_blastn_object(forward_primers, database, forward_out_file, 'forward')  # using amplicon_seq for testing purposes #expect 40 blast records and 1 hsp record
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file, 'reverse')  # using amplicon_seq for testing purposes #expect 40 blast records and 1 hsp record

    lo_forward_hsp_objects = create_hsp_objects(forward_blast)
    lo_reverse_hsp_objects = create_hsp_objects(reverse_blast)

    lo_queries = []

    for f_hsp_object in lo_forward_hsp_objects:
        for r_hsp_object in lo_reverse_hsp_objects:
            if f_hsp_object.name == r_hsp_object.name: #if the f_hsp and r_hsp are from the same primer query
                if f_hsp_object.contig_name == r_hsp_object.contig_name: #f_hsp and r_hsp are on the same contig
                    if pcr_directly(f_hsp_object, r_hsp_object, database, amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects):
                        lo_f_r = [f_hsp_object, r_hsp_object]
                        lo_queries.append(lo_f_r)

                else: # primer located on different contigs
                    # print(f_hsp_object.name, r_hsp_object.name)
                    # print('f contig', f_hsp_object.contig_name)
                    # print('r contig', r_hsp_object.contig_name)
                    blast_object = create_blastn_full_search_object(amplicon_sequences, database, full_out_file, 'full') #expect 1 blast record and 2 hsp records
                    lo_hsp_objects = create_hsp_objects(blast_object)
                    valid_strands(f_hsp_object, r_hsp_object, lo_forward_hsp_objects, lo_reverse_hsp_objects)
                    lo_correct_distance_queries = entire_gene(amplicon_sequences, blast_object, lo_hsp_objects, lo_forward_hsp_objects, lo_reverse_hsp_objects) #TODO: change me to one hsp object
                    lo_queries.append(lo_correct_distance_queries)

    print(lo_queries)

def main():

    database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta" #contains id and a complete genome sequence
    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #contains primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #contains primer id's and primer sequences
    forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml" #location where the blast record from comparing the forward primers to the db should go
    reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml" #location where the blast record from comparing the reverse primers to the db should go
    full_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml"

    cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta" #complete
    cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
    cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta" #complete
    cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"

   # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

    test_contig_51_bp_removed = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_51_middle_bp_removed.fasta"

    # using an amplicon_sequence as database for testing purposes
    #takes in forward_primers and reverse_primers (queries), database (contigs or complete) [only one database!!!], full amplicon sequences that contain the forward and reverse primers inputted,
    # a forward and reverse out file where the blast results can go, a complete amplicon
    return(pcr_prediction(forward_primers, reverse_primers, cj0483_complete_amp_seq, forward_out_file, reverse_out_file, amplicon_sequences, full_out_file)) #TODO: using amplicon sequence as db for testing purposes


if __name__ == "__main__": main()

#TODO: Consider primer cj1134 (actually from gene cj1324)
#TODO: Is it worth having Blastn as a class? (thought for later)
#TODO: Put functions into HSP that might be useful (thought for later)
#TODO: Is createObjects in the right spot?
#TODO: Change main to input files (currently embedded in main)? or somehow create a way to input files?
#TODO: change comments to be in python format (docstring)