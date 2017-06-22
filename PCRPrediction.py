from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Blastn import Blastn
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from HSP import HSP
import os
import errno
from collections import defaultdict

#TODO: determine all threshold values
E_VALUE_THRESHOLD = 0.04
MAX_DIST_BTWN_PRIMERS = 50 #<=
CUTOFF_GENE_LENGTH = 60
SNP_THRESHOLD = 5 #5 bp must be an exact match with 3' end of primer
VALID_DIR_THRESHOLD = 2 # The amount of bp's that can be located on end/start of db primer before reaching the end/start of amp

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
    :param out_file: A string that contains the file location of the xml file being outputted
    :restrictions: Database is formatted using makeblastdb
    :return: Blastn object
    """
    blastn_object = Blastn()
    blastn_query(query_genes, database, out_file)
    blastn_object.create_blast_records(out_file)
    #changed!!!
    blastn_object.create_hsp_objects(query_genes)
    # blastn_object.create_hsp_records(query_genes)
    return blastn_object

#changed!!!
# def create_hsp_objects(blastn_object):
#     """ Creates and initializes all fields of hsp objects from blastn_object input
#
#     :param blastn_object: A Blastn Object
#     :return: HSP Object
#     """
#     hsp_objects = []
#
#     for blast_record in blastn_object.blast_records:
#         for alignment in blast_record.alignments:
#             for hsp in alignment.hsps:
#                 if alignment in blastn_object.hsp_records:
#                     # if hsp.expect < E_VALUE_THRESHOLD: #removed this check b/c passed evalue into blastn query
#                     hsp_name = blast_record.query
#                     hsp_object = HSP(hsp_name)  #creates a HSP object
#                     hsp_object.start = hsp.sbjct_start
#                     hsp_object.end = hsp.sbjct_end
#                     hsp_object.query_start = hsp.query_start
#                     hsp_object.query_end = hsp.query_end
#                     hsp_object.contig_name = alignment.hit_def
#
#                     # assuming no contigs (complete genome)
#                     if hsp.sbjct_start < hsp.sbjct_end:
#                         hsp_object.strand = True
#                     elif hsp.sbjct_start > hsp.sbjct_end:
#                         hsp_object.strand = False
#                     assert hsp.sbjct_start != hsp.sbjct_end
#
#                     hsp_object.length = abs(hsp_object.end - hsp_object.start) + 1
#                     hsp_object.db_length = alignment.length
#                     hsp_object.expect = hsp.expect
#                     hsp_object.sbjct = hsp.sbjct
#                     hsp_object.query = hsp.query
#
#                     hsp_objects.append(hsp_object)
#     print('hsp objects', hsp_objects)
#     return hsp_objects

def valid_strands(first_hsp_object, second_hsp_object):
    """ Assigns valid attributes to first and second hsp object and modifies lo_first and lo_second hsp objects s.t. they only contain strands with the correct orientation to one another.

    :param first_hsp_object: A HSP object to compare with second_hsp_object
    :param second_hsp_object: A HSP object to compare with first_hsp_object
    :param lo_forward_hsp_objects: A list of hsp objects that was run using forward primer sets
    :param lo_reverse_hsp_objects: A list of hsp objects that was run using reverse primer sets
    :return: None
    """

    if first_hsp_object.name == second_hsp_object.name:
        if (first_hsp_object.strand or second_hsp_object.strand) and not (first_hsp_object.strand and second_hsp_object.strand):
            first_hsp_object.valid = True
            second_hsp_object.valid = True
        else:
            first_hsp_object.valid = False
            second_hsp_object.valid = False

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


def is_distance(f_hsp_object, r_hsp_object, amplicon_sequences):
    """ Return True if f_hsp_object and r_hsp_object are within MAX_DIST_BTWN_PRIMERS to each other.

    :param f_hsp_object: A HSP Object on a contig
    :param r_hsp_object: A HSP Object on a contig
    :param amplicon_sequences: A Fasta file that contains amplicon sequences for comparison purposes.
    :restrictions: f_hsp_object and r_hsp_object must be on the same contig.
    :return: A Boolean
    """
    assert f_hsp_object.contig_name == r_hsp_object.contig_name
    distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        if f_hsp_object.name in amplicon.id:
            amplicon_length = len(amplicon.seq)
            return (abs(amplicon_length - distance) <= MAX_DIST_BTWN_PRIMERS)

def pcr_directly(forward_hsp_object, reverse_hsp_object, amplicon_sequences, f_primers, r_primers):
    """ Determines if forward and reverse hsp object are in the correct orientation and distance from each other to be considered a hit.

    :param forward_hsp_object: A HSP object in question found from a forward primer as query
    :param reverse_hsp_object: A HSP object in question found from a reverse primer as query
    :param amplicon_sequences: A Fasta file that contains amplicon sequences for comparison purposes.
    :param f_primers: A dictionary of forward primers with key as primer name and value as primer sequences
    :param r_primers: A dictionary of reverse primers with key as primer name and value as primer sequences
    :return: A Boolean
    """

    valid_strands(forward_hsp_object, reverse_hsp_object)
    is_snp_primer_search(forward_hsp_object, f_primers, r_primers)
    is_snp_primer_search(reverse_hsp_object, f_primers, r_primers)

    assert forward_hsp_object.valid == reverse_hsp_object.valid #should always be True because valid_strands assigns both objects the same valid attribute
    if forward_hsp_object.valid == True and forward_hsp_object.snp == False:
        return is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences)
    elif forward_hsp_object.valid == False or forward_hsp_object.snp == True:
        return False
    # else:
        #!!! forward hsp object is not initialized!!!

def is_snp_primer_search(hsp_object, f_primers, r_primers):
    """

    :param hsp_object:
    :param f_primers:
    :param r_primers:
    :return:
    """

    if hsp_object.name in f_primers and hsp_object.name in r_primers:
        f_primer = f_primers[hsp_object.name]
        r_primer = r_primers[hsp_object.name]
        forward_primer_seq = Seq(str(f_primer[len(f_primer) - SNP_THRESHOLD : len(f_primer)]), generic_dna)
        forward_primer_reverse_complement = forward_primer_seq.reverse_complement()
        reverse_primer_seq = Seq(str(r_primer[len(r_primer) - SNP_THRESHOLD : len(r_primer)]), generic_dna)
        reverse_primer_reverse_complement = reverse_primer_seq.reverse_complement()
        db_end_forward = hsp_object.sbjct[hsp_object.query_end - SNP_THRESHOLD : ]
        db_end_reverse = hsp_object.sbjct[hsp_object.query_end : hsp_object.query_end + SNP_THRESHOLD]

        if str(forward_primer_seq) in db_end_forward\
                or str(forward_primer_reverse_complement) in db_end_forward\
                or str(reverse_primer_seq) in db_end_forward\
                or str(reverse_primer_reverse_complement) in db_end_forward:
            hsp_object.snp = False
        elif str(forward_primer_seq) in db_end_reverse\
                or str(forward_primer_reverse_complement) in db_end_reverse\
                or str(reverse_primer_seq) in db_end_reverse\
                or str(reverse_primer_reverse_complement) in db_end_reverse:
            hsp_object.snp = False
        else:
            hsp_object.snp = True


#TODO: make for entire gene?
#TODO: delete? ... Currently isn't used!!!
def is_snp(hsp_object, f_primers, r_primers):
    #TODO:!!!
    """Return True if a mismatch or removed bp is located in the first SNP_THRESHOLD bp's on the 3' end
    :param hsp_object: A HSP Object that is facing in the right direction (towards contig)
    :param f_primers: A dictionary containing full forward primer sequences.
    :param r_primers: A dictionary containing full reverse primer sequences.
    :return: A boolean: True if mismatch or removed bp on 3' end, False otherwise
    """

    #TODO: Test that this will not work if the primer sequences do not correspond to the hsp_object
    #TODO: "_" only in entire_gene calls?
    #TODO: this only works when there is no deletions in the 5' end of the primer!

    if "_" in hsp_object.name: #handles when entire gene called with amp seq names
        assert("_" in hsp_object.name)
        hsp_object.set_name(hsp_object.name.partition('_')[2])
    if hsp_object.name in f_primers and hsp_object.name in r_primers:
        f_primer = f_primers[hsp_object.name]
        r_primer = r_primers[hsp_object.name]
        forward_primer_seq = Seq(str(f_primer[len(f_primer) - SNP_THRESHOLD : len(f_primer)]), generic_dna)
        forward_primer_reverse_complement = forward_primer_seq.reverse_complement()
        reverse_primer_seq = Seq(str(r_primer[len(r_primer) - SNP_THRESHOLD : len(r_primer)]), generic_dna)
        reverse_primer_reverse_complement = reverse_primer_seq.reverse_complement()

        starting_forward = hsp_object.sbjct[:len(f_primer)] # [len(f_primer) - SNP_THRESHOLD : len(f_primer)]
        starting_reverse = hsp_object.sbjct[:len(r_primer)] # [len(r_primer) - SNP_THRESHOLD : len(r_primer)]
        ending_forward = hsp_object.sbjct[len(hsp_object.sbjct) - len(f_primer) : len(hsp_object.sbjct) - len(f_primer) + SNP_THRESHOLD]
        ending_reverse = hsp_object.sbjct[len(hsp_object.sbjct) - len(r_primer) : len(hsp_object.sbjct) - len(r_primer) + SNP_THRESHOLD]

        # print(hsp_object.sbjct)
        # print(forward_primer_seq)
        # print(forward_primer_reverse_complement)
        # print(reverse_primer_seq)
        # print(reverse_primer_reverse_complement)
        # print('normal forward', starting_forward)
        # print('comp forward', ending_forward)
        # print('normal reverse', ending_reverse)
        # print('reverse comp reverse', starting_reverse)

        #does not consider insertions into the primer seq! If need to consider, make threshold value for # of insertions allowed ie. [:len(primer) + threshold]

        if forward_primer_seq or forward_primer_reverse_complement in starting_forward \
                or forward_primer_reverse_complement or forward_primer_seq in ending_forward\
                and hsp_object.strand == True:
            hsp_object.snp = False
            #TODO: hsp_object.primer = "forward" !!!?
        elif str(reverse_primer_seq) in ending_reverse \
                or str(reverse_primer_reverse_complement) in starting_reverse \
                or str(reverse_primer_reverse_complement) in ending_reverse\
                or str(reverse_primer_seq) in starting_reverse:
            hsp_object.snp = False
        else:
            hsp_object.snp = True

def valid_dir(lo_hsps):
    """ Modify lo_hsps to only contain primer sequences that are facing the end of the contig.

    :param lo_hsps: list of hsp objects
    :return: None
    """

    # valid_lo_hsps = [hsp for hsp in lo_hsps if (hsp.end == hsp.query_end or hsp.start == hsp.query_start)
    #                     and not (hsp.end == hsp.query_end and hsp.start == hsp.query_start)]


    for hsp in lo_hsps:
        # considers the case where the entire strand is found
        if (hsp.end == hsp.query_end or hsp.start == hsp.query_start) and not (hsp.end == hsp.query_end and hsp.start == hsp.query_start):
            if hsp.strand == False and abs(hsp.query_end - hsp.end) <= VALID_DIR_THRESHOLD:
                lo_hsps.remove(hsp)
            elif hsp.strand == True and abs(hsp.start - hsp.query_start) <= VALID_DIR_THRESHOLD:
                lo_hsps.remove(hsp)

#NOTE: f_object.name == r_object.name
#TODO: consider case with cj1134 and cj1324!!!
#TODO: This assumes that it would be highly unlikely that a gene will be found of long enough length on a database more than once!!!
#TODO: possibly just consider if you found it on one primer in the right direction and do not even consider both found??
def entire_gene(blast_object, reference_object, f_primers, r_primers):
    """ Return the hsp objects from the same blast record query as reference_object that are of proper orientation and long enough length to be considered found.

    :param blast_object: A Blastn object with all fields initialized
    :param reference_object: A HSP object
    :param f_primers: A dictionary containing full forward primer sequences.
    :param r_primers: A dictionary containing full reverse primer sequences.
    :restrictions: Does not consider the case where lo_hsps has > 2 elements
    #TODO: It would be highly unlikely that a gene will be found of long enough length on a database more than once?
    :return: List of hsp objects
    """

    reference_list = [hsp for hsp in blast_object.hsp_objects if reference_object.name in hsp.name and not hsp.length <= CUTOFF_GENE_LENGTH]
    #consider % identity?
    assert len(reference_list) < 3 #Deal with this case later if needed
    valid_dir(reference_list)
    if len(reference_list) == 1:
        print("only one primer found for", reference_list[0].name, "on", reference_list[0].contig_name)
    return reference_list

    # changed for list comp.
    # lo_hsps = []
    #
    # for hsp_object in blast_object.hsp_objects:
    #     # print(hsp_object.contig_name)
    #     count = 0
    #     if reference_object.name in hsp_object.name:  # hsp_object.name may be longer than reference_object.name
    #         # is_snp(hsp_object, f_primers, r_primers)
    #         if not hsp_object.length <= CUTOFF_GENE_LENGTH: # and not hsp_object.snp:
    #             count += 1
    #             assert count != 3 #Deal with this case later if needed
    #             lo_hsps.append(hsp_object)
    #
    #     valid_dir(lo_hsps)
    # if len(lo_hsps) == 1:
    #     print("only one primer found for", lo_hsps[0].name, "on", lo_hsps[0].contig_name)
    # # elif len(lo_hsps) == 2:
    # #     print("both primers found for", lo_hsps[0].name, lo_hsps[1].name)
    # return lo_hsps

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

    #TODO: should forward and reverse blast have contents when looking for entire gene?! Write tests for it!
    forward_blast = create_blastn_object(forward_primers, database, forward_out_file)
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file)
    blast_object = create_blastn_object(amplicon_sequences, database, full_out_file)

    lo_queries = []
    dict_f_primers = create_primer_dict(forward_primers)
    dict_r_primers = create_primer_dict(reverse_primers)

    lo_tup_same_queries = [(f_hsp, r_hsp) for f_hsp in forward_blast.hsp_objects for r_hsp in reverse_blast.hsp_objects if f_hsp.name == r_hsp.name]
    lo_tup_pcr_directly = [tup for tup in lo_tup_same_queries if tup[0].contig_name == tup[1].contig_name]
    lo_tup_entire_gene = [tup for tup in lo_tup_same_queries if tup not in lo_tup_pcr_directly] #primers located on different contigs

    for f_r_hsp_object in lo_tup_pcr_directly:
        if pcr_directly(f_r_hsp_object[0], f_r_hsp_object[1], amplicon_sequences, dict_f_primers, dict_r_primers):
            lo_queries.append(list(f_r_hsp_object)) #TODO: change later to keep as tuple and possibly return hash

    for f_r_hsp_object in lo_tup_entire_gene:
        lo_hsps = entire_gene(blast_object, f_r_hsp_object[1], dict_f_primers, dict_r_primers)
        if len(lo_hsps) > 0:
            lo_queries.append(lo_hsps)

    # for f_hsp_object in forward_blast.hsp_objects:
    #     lo_r_same_query = [r_hsp for r_hsp in reverse_blast.hsp_objects if r_hsp.name == f_hsp_object.name]
    #     lo_r_pcr_directly = [r_hsp for r_hsp in lo_r_same_query if f_hsp_object.contig_name == r_hsp.contig_name] #from same query and on same contig
    #     lo_r_entire_gene = [hsp for hsp in lo_r_same_query if hsp not in lo_r_pcr_directly] #from same query and on different contigs
    #     for r_hsp_object in lo_r_pcr_directly:
    #         if pcr_directly(f_hsp_object, r_hsp_object, amplicon_sequences, dict_f_primers, dict_r_primers):
    #             lo_f_r = [f_hsp_object, r_hsp_object]
    #             lo_queries.append(lo_f_r)
    #     for r_hsp_object in lo_r_entire_gene:
    #         lo_tmp = entire_gene(blast_object, r_hsp_object, dict_f_primers, dict_r_primers)
    #         if len(lo_tmp) > 0:
    #             lo_queries.append(lo_tmp)


    return lo_queries

def main(db_directory, forward_primers, reverse_primers, amplicon_sequences):
    """

    :param db_directory: The location of the directory with fasta database files contained (already formatted using makeblastdb)
    :param forward_primers: The fasta file location of the forward primers
    :param reverse_primers: The fasta file location of the reverse primers
    :param amplicon_sequences: The fasta file location of the amplicon sequences
    :return: A dictionary !!!
    """

    files = [file for file in os.listdir(db_directory) if file.endswith('.fasta')]
    files_paths = []
    for file in files:
        files_paths.append(os.path.abspath(db_directory) + '/' + file)

    #create new folder for out files
    try:
        os.mkdir(db_directory + "/out_files")
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    pcr_predictions_dict = {}
    for file_path in files_paths:
        print('FILE!!!')
        name = file_path.partition(db_directory + "/")[2]
        f_out_file_path = db_directory + "/out_files/" + "f_" + name.replace("fasta", "xml")
        r_out_file_path = db_directory + "/out_files/" + "r_" + name.replace("fasta", "xml")
        full_out_file_path = db_directory + "/out_files/" + "full_" + name.replace("fasta", "xml")
        result = pcr_prediction(forward_primers, reverse_primers, file_path, f_out_file_path, r_out_file_path, amplicon_sequences, full_out_file_path)
        pcr_predictions_dict[name] = result

    print('pcr predictions dict', pcr_predictions_dict)
    return(pcr_predictions_dict)


forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #fasta file with primer id's and primer sequences
reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #fasta file with primer id's and primer sequences
amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

if __name__ == "__main__":
    test_pcr_prediction = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction"
    test_bp_removed = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed"
    test_gene_annotation = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
    main(test_pcr_prediction, forward_primers, reverse_primers, amplicon_sequences)






# #Test amp seq
# cj0008_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0008.fasta" #complete
# cj0483_contig_trunc_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_truncation.fasta"
# cj0483_complete_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_complete1.fasta" #complete
# cj0483_contig_full_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full.fasta"
# cj1134_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj1134.fasta"
# cj1324_amp_seq = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj1324.fasta"
#
# # cj0483_contig_full_amp_seq_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/cj0483_contig_full_db" #contigs (no missing contigs inbetween)
#
# #test data
# database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta" #contains id and a complete genome sequence
# test_contig_51_bp_removed = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_db/cj0483_51_middle_bp_removed.fasta"
# test_mystery_db_name_same_contig = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db.fasta"
# test_one_contig = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/cj0483_primer_same_contig_with_2_contigs.fasta"
# test_mystery_db_name_entire_gene = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction/mystery_db_entire_gene_61_bp_each_contig.fasta"
# cj0008_rm_2_ws = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed/cj0008_rm_2_ws.fasta"

#TODO: With database input, makeblastdb or run through command line?