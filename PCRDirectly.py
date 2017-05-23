from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Blastn import Blastn

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
HSP_THRESHOLD = 1

#create a blastn query and output into an xml file
def blastn_query(query_genes, database, out_file):
    #TODO: Formatted using makeblastdb : query_genes, database (fasta), out_file (xml)

    #TODO: determine a word_size
    #query the database using blastn
    blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=4, outfmt=5, out=out_file)
    stdout, stderr = blastn_cline()

def create_blast_record(str):
    result_handle = open(str)
    blast_records = NCBIXML.parse(result_handle)  # returns an iterator
    blast_records = list(blast_records)

    return blast_records

#filters blast_records to ensure the hsp have a long enough alignment. Produces a dictionary of hsp that meet threshold value.
#primers is a fasta file that contains the gene sequence to compare the blast_records to.
def create_hsp_records(blast_records, primers):
    dict_hsp = {}

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for primer in SeqIO.parse(primers, "fasta"):
                if alignment.title.find(primer.name) != -1:
                    lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD <= (hsp.identities / len(primer.seq))] #TODO: bug in % identities...

                    if (len(lo_hsp) != 0): #lo_hsp is not empty
                        #print('lo_hsp', lo_hsp)
                        dict_hsp[alignment] = lo_hsp
                    break #don't need to go through all primer seq'n in primers once the right alignment is found.
    #print('length dict_hsp:', len(dict_hsp))
    return dict_hsp

def create_blastn_object(query_genes, database, out_file):
    blastn_query(query_genes, database, out_file)
    blast_records = create_blast_record(out_file)
    hsp_records = create_hsp_records(blast_records, query_genes)
    blastn_object = Blastn(blast_records, hsp_records)
    return blastn_object


#returns a dictionary where key is the query title and value is a boolean that is:
#   True if primer located on leading strand (+, template) and false if primer located on lagging strand (-, complement)
def find_strands(blastn):
    strands = {}

    for blast_record in blastn.blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if alignment in blastn.hsp_records:
                    if hsp in blastn.hsp_records[alignment]:  # high scoring pairs
                        if hsp.expect < E_VALUE_THRESHOLD:
                            query_title = blast_record.query
                            # assuming no contigs (complete genome)
                            if hsp.sbjct_start < hsp.sbjct_end:
                                strands[query_title] = True
                            elif hsp.sbjct_start > hsp.sbjct_end:
                                strands[query_title] = False
                                # hsp.sbjct_end should never = hsp.sbjct_start

    return (strands)


#Given a dictionary of forward and a dict of reverse primers that contains their query title and a bool (true = leading strand, false = lagging strand)
# returns a list of query_title's that contain primer strand pairs that are opposite (one leading, one lagging)
def valid_strands(forward_primers, reverse_primers):
    lo_valid_strand_ids = []
    for primer_id in forward_primers:
        if primer_id in reverse_primers:
            if reverse_primers[primer_id] or forward_primers[primer_id] and not (reverse_primers[primer_id] and forward_primers[primer_id]): #exclusive or... couldn't find non-bitwise xor
                lo_valid_strand_ids.append(primer_id)
    return (lo_valid_strand_ids)


#returns a dictionary of the primers in lo_primers where key is the title of the primer and value is the starting index of the primer.
def get_primer_start(blastn, lo_primers):
    primers = {}

    for blast_record in blastn.blast_records:
        if blast_record.query in lo_primers:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if alignment in blastn.hsp_records:
                        if hsp in blastn.hsp_records[alignment]:  # high scoring pairs
                            if hsp.expect < E_VALUE_THRESHOLD:
                                primer = hsp.sbjct_start
                                primer_title = blast_record.query
                                primers[primer_title] = primer
    return primers

#returns primer strands in lo_valid_strand_direction that have the correct distance from each other, compared to the amplicon_sequences.
def find_distances(lo_valid_strand_direction, forward_blastn, reverse_blastn, amplicon_sequences):

    lo_correct_distance_primers = []
    forward_primers = get_primer_start(forward_blastn, lo_valid_strand_direction) #only checks the strands in opposite direction (valid_strands)
    reverse_primers = get_primer_start(reverse_blastn, lo_valid_strand_direction)

    #start_primers and end_primers should have equal keys.
    for primer_key in forward_primers:
        if primer_key in reverse_primers:    #start_primers and end_primers should have equal keys (should always be true)... delete?
            query_length = abs(forward_primers[primer_key] - reverse_primers[primer_key]) + 1
            for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
                # naively assuming that amplicon_sequences contains 5 int (ex. 11169_) before needed id that can be compared to primer key
                amplicon_id = amplicon.id[6: ]
                if amplicon_id == primer_key:
                    amplicon_length = len(amplicon.seq)
                    if query_length == amplicon_length:
                        lo_correct_distance_primers.append(primer_key)

    print('correct_distance_primers', lo_correct_distance_primers)
    print('lengthof_correct_distance_primers', len(lo_correct_distance_primers))
    return lo_correct_distance_primers


def pcr_directly(forward_primers, reverse_primers, amplicon_sequences, forward_out_file, reverse_out_file):

    forward_blast = create_blastn_object(forward_primers, amplicon_sequences, forward_out_file)
    reverse_blast = create_blastn_object(reverse_primers, amplicon_sequences, reverse_out_file)

    lo_forward_strands = find_strands(forward_blast)
    lo_reverse_strands = find_strands(reverse_blast)

    lo_valid_strand_direction = valid_strands(lo_forward_strands, lo_reverse_strands)
    return find_distances(lo_valid_strand_direction, forward_blast, reverse_blast, amplicon_sequences)


def main():

    #assume database is complete or primers are on the same contigs
    database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta" #contains id and a complete genome sequence
    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #contains primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #contains primer id's and primer sequences
    forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml" #location where the blast record from comparing the forward primers to the db should go
    reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml" #location where the blast record from comparing the reverse primers to the db should go
    amplicon_sequence = "/home/sfisher/Sequences/amplicon_sequences/cj0008.fasta"
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"


    # using amplicon_sequences as database for testing purposes
    pcr_directly(forward_primers, reverse_primers, amplicon_sequences, forward_out_file, reverse_out_file)

    print( "not done")

if __name__ == "__main__": main()

#TODO: Consider primer cj1134 (actually from gene cj1324)
#TODO: memoize?
#TODO: merge SearchGeneDirectly with PCRDirectly