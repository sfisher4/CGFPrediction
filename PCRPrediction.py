from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Blastn import Blastn
from HSP import HSP

E_VALUE_THRESHOLD = 0.04 #TODO: determine an e-value threshold
HSP_THRESHOLD = 1

#create a blastn query and output into an xml file
def blastn_query(query_genes, database, out_file):
    #TODO: Formatted using makeblastdb : query_genes, database (fasta), out_file (xml)

    #TODO: determine a word_size
    #query the database using blastn
    blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=4, outfmt=5, out=out_file)
    stdout, stderr = blastn_cline()

def create_blastn_object(query_genes, database, out_file, name):
    blastn_object = Blastn(name)
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    print('blast created', blastn_object.blast_records)
    blastn_object.create_hsp_records(query_genes)
    return blastn_object

#Creates a hsp object by assigning a name, start, end, length and strand.
def create_hsp_objects(blastn_object):

    hsp_objects = []
    count = 0

    for blast_record in blastn_object.blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if alignment in blastn_object.hsp_records:
                    if hsp.expect < E_VALUE_THRESHOLD:
                        print('hsp expect',hsp.expect)
                        hsp_name = blast_record.query + "_" + str(count)
                        print("name" ,hsp_name)
                        print("hit name", alignment.hit_def)
                        hsp_object = HSP(hsp_name)  #TODO: make sure you can have >1 object with the same name... seems to be fine
                        hsp_object.start = hsp.sbjct_start
                        hsp_object.end = hsp.sbjct_end

                        # assuming no contigs (complete genome)
                        if hsp.sbjct_start < hsp.sbjct_end:
                            hsp_object.strand = True
                        elif hsp.sbjct_start > hsp.sbjct_end:
                            hsp_object.strand = False
                        assert hsp.sbjct_start != hsp.sbjct_end

                        hsp_object.length = (hsp_object.end - hsp_object.start) + 1
                        hsp_object.db_length = alignment.length
                        #print('db length', alignment.length)

                        for description in blast_record.descriptions:
                            #print('blast record title', description.title)
                            hsp_object.blast_record_title = description. title

                        print('add one')
                        hsp_objects.append(hsp_object)
                        count += 1

    print('hsp objects', hsp_objects)
    return hsp_objects

#TODO: test valid_strands using a database that finds >1 primer sequence but doesn't face each other NOTE: checked in assert statement in find_distances
# Given a list of forward and reverse hsp_objects, modifies the lists to include only hsp_objects that are in the correct
# strand orientation and assigns the objects valid attributes (ie. one leading (True) and one lagging(False)).
def valid_strands(lo_forward_hsp_objects, lo_reverse_hsp_objects):

    for f_primer in lo_forward_hsp_objects:
        for r_primer in lo_reverse_hsp_objects:
            if f_primer.name == r_primer.name:
                if f_primer.strand or r_primer.strand and not (f_primer.strand and r_primer.strand):
                    f_primer.valid = True
                    r_primer.valid = True
                else:
                    f_primer.valid = False
                    r_primer.valid = False
                    lo_forward_hsp_objects.remove(f_primer)
                    lo_reverse_hsp_objects.remove(r_primer)

#returns primer strands in lo_valid_strand_direction that have the correct distance from each other, compared to the amplicon_sequences (includes primer strand length)
#TODO: consider with variable length between contigs!
def find_distances(amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects, database, out_file, full_amp):

    lo_correct_distance_primers = []

    for f_hsp_object in lo_forward_hsp_objects:
        assert f_hsp_object.valid == True
        for r_hsp_object in lo_reverse_hsp_objects:
            assert r_hsp_object.valid == True

            #located on same contig:
            if f_hsp_object.blast_record_title == r_hsp_object.blast_record_title:
                q_distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
                if is_correct_distance(amplicon_sequences, f_hsp_object, q_distance):
                    lo_correct_distance_primers.append(f_hsp_object.name)
                    print('lo correct distance primers',lo_correct_distance_primers)  #Note: Should not reach here when inputting with primers on different contigs

            #     search_gene_directly(full_amp, database, out_file)
            else: #located on different contigs directly next to each other!
                print('hsp_object_name', f_hsp_object.name)
                search_gene_directly(full_amp, database, out_file, f_hsp_object.name)
                # #TODO: Consider this case in the case where primers are located on different contigs
                # if f_hsp_object.strand == True:
                #     length_leading_primer = f_hsp_object.end - f_hsp_object.start
                #     length_leading_sequence = f_hsp_object.db_length
                #     length_lagging_primer = r_hsp_object.start - r_hsp_object.end
                #     q_distance = length_leading_primer + (length_leading_sequence - f_hsp_object.end + 1) + length_lagging_primer + r_hsp_object.end
                #
                # else:
                #     assert(f_hsp_object.strand == False)
                #     length_leading_primer = f_hsp_object.start - f_hsp_object.end
                #     length_leading_sequence = r_hsp_object.db_length
                #     length_lagging_primer = r_hsp_object.end - r_hsp_object.start
                #     q_distance = length_leading_primer + (length_leading_sequence - r_hsp_object.end + 1) + length_lagging_primer + f_hsp_object.end


    print('correct_distance_primers', lo_correct_distance_primers)
    print('length of_correct_distance_primers', len(lo_correct_distance_primers))
    return lo_correct_distance_primers

#takes in a hsp object and determines
def is_correct_distance(amplicon_sequence, f_hsp_object, q_distance):
    for amplicon in SeqIO.parse(amplicon_sequence, "fasta"):
        amplicon_id = amplicon.id[6:]
        if f_hsp_object.name == amplicon_id: #TODO: change after changing hsp_object.name
            amplicon_length = len(amplicon.seq)
            print('amp length', amplicon_length)  # TODO: change to take >1 contig into account!
            print('q distance', q_distance)
            if q_distance == amplicon_length:
                return True


def search_gene_directly(full_amp, database, out_file, hsp_primer_object_name):
    blastn_object = create_blastn_object(full_amp, database, out_file, 'full_amp')
    print('len blast records',len(blastn_object.blast_records))
    #TODO: check to make sure you are only checking the sequences that starts and end at the primer sequence you were looking for!!!
    # for blast_record in blastn_object.blast_records:
    #     print('query', blast_record.query)
    hsp_objects = create_hsp_objects(blastn_object)
    for hsp_object in hsp_objects:
        print(hsp_object.length)


def pcr_prediction(forward_primers, reverse_primers, database, amplicon_sequences, forward_out_file, reverse_out_file, full_out_file, full_amp):

    forward_blast = create_blastn_object(forward_primers, database, forward_out_file, 'forward') #using amplicon_seq for testing purposes
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file, 'reverse') #using amplicon_seq for testing purposes

    lo_forward_hsp_objects = create_hsp_objects(forward_blast)
    print('\n')
    lo_reverse_hsp_objects = create_hsp_objects(reverse_blast)

    valid_strands(lo_forward_hsp_objects, lo_reverse_hsp_objects)
    return find_distances(amplicon_sequences, lo_forward_hsp_objects, lo_reverse_hsp_objects, database, full_out_file, full_amp)


def main():

    #assume database is complete or primers are on the same contigs
    database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta" #contains id and a complete genome sequence
    forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #contains primer id's and primer sequences
    reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #contains primer id's and primer sequences
    forward_out_file = "/home/sfisher/Sequences/blast_record/forward_primers_blast.xml" #location where the blast record from comparing the forward primers to the db should go
    reverse_out_file = "/home/sfisher/Sequences/blast_record/reverse_primers_blast.xml" #location where the blast record from comparing the reverse primers to the db should go
    full_out_file = "/home/sfisher/Sequences/blast_record/full_blast.xml" #TODO: create file?

    cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta" #complete
    cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
    cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta" #complete
    cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"

   # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)


    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"


    # using amplicon_sequences as database for testing purposes
    #takes in forward_primers and reverse_primers (queries), database (contigs or complete), full amplicon sequences that contain the forward and reverse primers inputted,
    # a forward and reverse out file where the blast results can go, a complete amplicon
    return(pcr_prediction(forward_primers, reverse_primers, cj0483_contig_full_amp_seq, amplicon_sequences, forward_out_file, reverse_out_file, full_out_file, cj0483_complete_amp_seq)) #TODO: using amplicon sequence as db for testing purposes


if __name__ == "__main__": main()

#TODO: Consider primer cj1134 (actually from gene cj1324)
#TODO: Consider 2 contigs with variable length between
#TODO: possibly rename PCRDirectly
#TODO: Modify SearchGeneDirectly to call according to what I need for getDirection
#TODO: Is it worth having Blastn as a class?
#TODO: Put functions into HSP that might be useful
#TODO: Is createObjects in the right spot?
#TODO: Determine how to input files (currently in main)