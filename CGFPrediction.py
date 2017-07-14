from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Blastn import Blastn
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from HSP import HSP
import os
import errno
import itertools
from collections import defaultdict

E_VALUE_THRESHOLD = 0.04
MAX_DIST_BTWN_PRIMERS = 50 #<=
CUTOFF_GENE_LENGTH = 70
SNP_THRESHOLD = 5 #5 bp must be an exact match with 3' end of primer
MAX_MM_CONTIG_END = 10 # The amount of bp's that can be located on end/start of db primer before reaching the end/start of amp
PERC_ID_THRESH = 80.0
WORD_SIZE = 7 #switched from 4
QCOV_HSP_PERC = 80.0

def blastn_query(query_genes, database, out_file, type):
    """ Outputs a blastn query into an xml file.

    :param query_genes: A fasta file that contains the query genes that are being searched in the database
    :param database: A fasta file that contains the database that is being searched against
    :param out_file: A xml file that the blastn query
    :restrictions: Database is formatted using makeblastdb
    :return: None
    """

    if type == True:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5, out=out_file,
                                             evalue=E_VALUE_THRESHOLD, perc_identity=PERC_ID_THRESH, qcov_hsp_perc=QCOV_HSP_PERC)
    else:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5,
                                             out=out_file, evalue=E_VALUE_THRESHOLD, perc_identity=PERC_ID_THRESH)
    stdout, stderr = blastn_cline()

def create_blastn_object(query_genes, database, out_file, type=False):
    """ Return a blastn object with initialized blast_records and hsp_records

    :param query_genes: A fasta file that contains the query genes that are being searched in the database
    :param database: A fasta file that contains the database that is being searched against
    :param out_file: A string that contains the file location of the xml file being outputted
    :restrictions: Database is formatted using makeblastdb
    :return: Blastn object
    """
    blastn_object = Blastn()
    blastn_query(query_genes, database, out_file, type)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_objects(query_genes)
    # blastn_object.create_hsp_records(query_genes)
    return blastn_object

def valid_strands(first_hsp_object: HSP, second_hsp_object: HSP) -> None :
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

def pcr_directly(forward_hsp_object:HSP, reverse_hsp_object:HSP, amplicon_sequences, f_primers:dict, r_primers:dict) -> bool:
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

    if forward_hsp_object.snp == True:
        print("There is a 3' SNP on", forward_hsp_object.name)
    if reverse_hsp_object.snp == True:
        print("There is a 3' SNP on", reverse_hsp_object.name)

    assert forward_hsp_object.valid == reverse_hsp_object.valid #should always be True because valid_strands assigns both objects the same valid attribute
    if forward_hsp_object.valid == True : # and forward_hsp_object.snp == False:
        return is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences)
    elif forward_hsp_object.valid == False : # or (forward_hsp_object.snp == True and reverse_hsp_object == True):
        return False
    else:
        assert True == False

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

    #REMEMBER: if the hsp is on the lagging strand, end will be at the beginning of the strand
    # for hsp in lo_hsps:
    #     # considers the case where the entire strand is found
    #     print('contig', hsp.contig_name)
    #     print('end', hsp.end)
    #     print('start', hsp.start)
    #     print('strand', hsp.strand)
    #     print('sbjct', hsp.sbjct)
    #     print('len sbjct', len(hsp.sbjct))
    #     print('db length', hsp.db_length)
    #     print('diff', abs(hsp.end - hsp.db_length))
    #     #not covering the entire gene
    #     if hsp.strand == True:
    #         if not abs(hsp.end - hsp.db_length) <= MAX_MM_CONTIG_END and not hsp.start <= MAX_MM_CONTIG_END:
    #             #TODO: ensure if the sequence is on the lagging strand that this is still the case
    #             #the hsp is not facing the end of the contig
    #             if not abs(hsp.end - hsp.db_length) <= MAX_MM_CONTIG_END:
    #                 print('result diff', abs(hsp.end - hsp.db_length))
    #                 print('valid dir thres', MAX_MM_CONTIG_END)
    #                 print('contig name', hsp.contig_name)
    #                 lo_hsps.remove(hsp)
    #     elif hsp.strand == False:
    #         if not abs(hsp.start - hsp.db_length) <= MAX_MM_CONTIG_END and not hsp.end <= MAX_MM_CONTIG_END:
    #             if not abs(hsp.start - hsp.db_length) <= MAX_MM_CONTIG_END:
    #                 print('result diff', abs(hsp.start - hsp.db_length))
    #                 print('valid dir thres', MAX_MM_CONTIG_END)
    #                 print('contig name', hsp.contig_name)
    #                 lo_hsps.remove(hsp)

    #TODO: FIX ME!!!!
    # REMEMBER: if the hsp is on the lagging strand, end will be at the beginning of the strand
    count = 0
    for hsp in lo_hsps:
        print('contig', hsp.contig_name)
        print('end', hsp.end)
        print('start', hsp.start)
        print('sbjct', hsp.sbjct)
        print('len sbjct', len(hsp.sbjct))
        print('hsp end - length hsp sbjct', abs(hsp.end - len(hsp.sbjct)))

        # considers the case where the entire strand is found

        if not (abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END and hsp.start <= MAX_MM_CONTIG_END):
            # the the hsp is not facing the end of the contig
            # TODO: ensure if the sequence is on the lagging strand that this is still the case
            if abs(hsp.end - len(hsp.sbjct)) > MAX_MM_CONTIG_END:
                lo_hsps.remove(hsp)
            else:
                print('end align', abs(hsp.end - len(hsp.sbjct)))
        else:
            count += 1
            print('the hsp extends along the entire contig! This should happen 40 times')
    print ('count', count)

            # for hsp in lo_hsps:
    #     # considers the case where the entire strand is found
    #     if (hsp.end == hsp.query_end or hsp.start == hsp.query_start) and not (hsp.end == hsp.query_end and hsp.start == hsp.query_start):
    #         #Need to consider bp's that are removed here!!!
    #         if hsp.strand == False and abs(hsp.query_end - hsp.end) <= MAX_MM_CONTIG_END:
    #             lo_hsps.remove(hsp)
    #         elif hsp.strand == True and abs(hsp.start - hsp.query_start) <= MAX_MM_CONTIG_END:
    #             lo_hsps.remove(hsp)

#NOTE: f_object.name == r_object.name...
#TODO: consider case with cj1134 and cj1324!!!
#TODO: This assumes that it would be highly unlikely that a gene will be found of long enough length on a database more than once!!!
def entire_gene(blast_object:Blastn, reference_object:HSP, f_primers_dict:dict, r_primers_dict:dict) -> list:
    """ Return the hsp objects from the same blast record query as reference_object that are of proper orientation and long enough length to be considered found.

    :param blast_object: A Blastn object with all fields initialized
    :param reference_object: A HSP object
    :param f_primers_dict: A dictionary containing full forward primer sequences.
    :param r_primers_dict: A dictionary containing full reverse primer sequences.
    :restrictions: Does not consider the case where lo_hsps has > 2 elements
    :return: List of hsp objects
    """
    # print('hsp objects', blast_object.hsp_objects)
    # for hsp in blast_object.hsp_objects:
    #     if reference_object.name in hsp.name:
    #         print(hsp.name)
    #         print(hsp.length)

    reference_list = [hsp for hsp in blast_object.hsp_objects if reference_object.name in hsp.name and not hsp.length <= CUTOFF_GENE_LENGTH]
    # print('ref list 1', reference_list)
    # for hsp in reference_list:
    #     print(hsp.name)
    #     print(hsp.length)
    #     print(hsp.sbjct)
    #     print(hsp.query)
    valid_dir(reference_list)
    # print('ref list 2', reference_list)
    assert len(reference_list) < 3 #Deal with this case later if needed
    # if len(reference_list) == 1: #and called with primers on different contigs
        # print("only one primer found for", reference_list[0].name, "on", reference_list[0].contig_name)
    return reference_list

def ehybridization(blast_object:Blastn) -> list:
    lo_hsp = [hsp for hsp in blast_object.hsp_objects if hsp.length >= CUTOFF_GENE_LENGTH]
    return lo_hsp

def pcr_diff_contigs():
    """
        1) Primer found (is in f_primer or r_primer hsp objects)
        2) primer facing end of contig?
        3) 3'SNP
        4)
    :return:
    """

def same_contig_pred(blast_object:Blastn, lo_tup_same_contig:list, dict_f_primers:dict, dict_r_primers:dict) -> list:
    """ Finds query in db through ehybridization followed by epcr.
        WORKFLOW: ehybridization --> epcr

    :param blast_object: A Blastn object.
    :param lo_tup_same_contig: A list of tuples with forward and reverse hsp objects
    :param dict_f_primers: A dict containing forward primers
    :param dict_r_primers: A dict containing reverse primers
    :return: A list of lists of forward and reverse HSP objects.
    """

    lo_queries = []
    lo_hybridized_hsp = ehybridization(blast_object)

    #TODO: added for testing!!!
    lo_hybridized_hsp_names = []
    for hsp in lo_hybridized_hsp:
        lo_hybridized_hsp_names.append(hsp.name)
    # print('hybridized', lo_hybridized_hsp_names)


    lo_hsp_start = [hsp.start for hsp in lo_hybridized_hsp]
    lo_hsp_end = [hsp.end for hsp in lo_hybridized_hsp]
    for hsp_tup in lo_tup_same_contig:
        if any((hsp_tup[0].start == hsp_start or hsp_tup[1].start == hsp_start) for hsp_start in lo_hsp_start):
            if pcr_directly(hsp_tup[0], hsp_tup[1], amplicon_sequences, dict_f_primers, dict_r_primers):
                lo_queries.append(list(hsp_tup))
        elif any ((hsp_tup[0].end == hsp_end or hsp_tup[1].end == hsp_end) for hsp_end in lo_hsp_end):
            if pcr_directly(hsp_tup[0], hsp_tup[1], amplicon_sequences, dict_f_primers, dict_r_primers):
                lo_queries.append(list(hsp_tup))

   #TODO: added for testing!!!
    lo_hsp_names = []
    for quer in lo_queries:
        for hsp in quer:
            lo_hsp_names.append(hsp.name)
    # print('lo_queries', lo_hsp_names)
    # print('diff', len(lo_hsp_names) - 2 * len(lo_hybridized_hsp_names))

    #TODO: added hybridized results for testing purposes. delete from return value later
    return (lo_queries, lo_hybridized_hsp_names)

def diff_contig_pred(blast_object:Blastn, lo_tup_diff_contig:list, dict_f_primers:dict, dict_r_primers:dict):

    #f_r_hsp_object should not have duplicates when primer is on both contigs!!! only send the reference object once!!!
    """ Primers located on different contigs.
        WORKFLOW: entire_gene -> epcr
        1) both forward and reverse primers found (lo_tup_diff_contig) found
        2) just forward or just reverse primer found
        3)
    """
    lo_queries = []

    for f_r_hsp_object in lo_tup_diff_contig:
        lo_hsps = entire_gene(blast_object, f_r_hsp_object[1], dict_f_primers, dict_r_primers)
        if len(lo_hsps) > 0:
            lo_queries.append(lo_hsps)
    return lo_queries


#TODO: Change return value to hash?
def pcr_prediction(forward_primers:str, reverse_primers:str, database:str, forward_out_file:str, reverse_out_file:str, amplicon_sequences:str, full_out_file:str):
    """ Return a list of hsp's that would be predicted as hsp's in vitro...

    :param forward_primers: A Fasta file location containing forward query primers (ids and sequences)
    :param reverse_primers: A Fasta file location containing reverse query primers (ids and sequences)
    :param database: A Fasta file path (formatted using makeblastdb) that contains a database (complete or draft)
    :param forward_out_file: An xml file location to store blastn results from forward primers as queries
    :param reverse_out_file: An xml file location to store blastn results from reverse primers as queries
    :param amplicon_sequences: A Fasta file containing amplicon sequences
    :param full_out_file: An xml file location to store blastn results from amplicon sequences as queries
    :return: List of HSP's that were found in database.
    """

    forward_blast = create_blastn_object(forward_primers, database, forward_out_file, True)
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file, True)
    blast_object = create_blastn_object(amplicon_sequences, database, full_out_file)
    full_blast_qcov = create_blastn_object(amplicon_sequences, database, full_out_file, True)

    dict_f_primers = create_primer_dict(forward_primers)
    dict_r_primers = create_primer_dict(reverse_primers)

    lo_tup_same_queries = [(f_hsp, r_hsp) for f_hsp in forward_blast.hsp_objects for r_hsp in reverse_blast.hsp_objects if f_hsp.name == r_hsp.name]
    lo_tup_same_contig = [tup for tup in lo_tup_same_queries if tup[0].contig_name == tup[1].contig_name]
    lo_tup_diff_contig = [tup for tup in lo_tup_same_queries if tup not in lo_tup_same_contig]

    #TODO: added for testing purposes. Can delete later.
    lo_hsp_diff_contigs = [hsp[0].name for hsp in lo_tup_diff_contig]
    print('tuples on different contigs', database, lo_hsp_diff_contigs)

    # TODO: this considers the case where one primer sequence is found. Test me.
    try:
        lo_f_primers, lo_r_primers = zip(*lo_tup_same_queries)
    except ValueError:
        lo_f_primers = []
        lo_r_primers = []

    f_hsp_single_primers = (hsp.name for hsp in forward_blast.hsp_objects if hsp not in lo_f_primers)
    r_hsp_single_primers = (hsp.name for hsp in reverse_blast.hsp_objects if hsp not in lo_r_primers)
    lo_single_primers = list(itertools.chain(f_hsp_single_primers, r_hsp_single_primers))
    print('single primers found', lo_single_primers)
    #predicts when primers are located on the same contig.
    lo_same_contig_queries_tup = same_contig_pred(full_blast_qcov, lo_tup_same_contig, dict_f_primers, dict_r_primers)
    lo_same_contig_queries = lo_same_contig_queries_tup[0]
    lo_hybridized_results = lo_same_contig_queries_tup[1]

    #TODO: uncomment this and delete line: lo_queries = lo_same_contig_queries when ready to test for primers on different contigs
    # lo_diff_contig_queries = diff_contig_pred(blast_object, lo_tup_diff_contig, dict_f_primers, dict_r_primers)
    # lo_queries = list(itertools.chain(lo_same_contig_queries, lo_diff_contig_queries))

    lo_queries = lo_same_contig_queries
    # print(lo_queries)
    # for quer in lo_queries:
    #     print(quer)
    #     for hsp in quer:
    #         print(hsp.name)

    #TODO: added return value lo_hsp_diff_contigs for testing purposes!!! delete later.
    return_list = (lo_queries, lo_hsp_diff_contigs, lo_hybridized_results, lo_single_primers)
    return return_list

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

    cgf_predictions_dict = {}
    for file_path in files_paths:
        print('FILE!!!', file_path)
        file_name = file_path.partition(db_directory + "/")[2]
        print(file_name)
        f_out_file_path = db_directory + "/out_files/" + "f_" + file_name.replace("fasta", "xml")
        r_out_file_path = db_directory + "/out_files/" + "r_" + file_name.replace("fasta", "xml")
        full_out_file_path = db_directory + "/out_files/" + "full_" + file_name.replace("fasta", "xml")

        result = pcr_prediction(forward_primers, reverse_primers, file_path, f_out_file_path, r_out_file_path, amplicon_sequences, full_out_file_path)[0]
        cgf_predictions_dict[file_name] = result

    return(cgf_predictions_dict)


forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #fasta file with primer id's and primer sequences
reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #fasta file with primer id's and primer sequences
amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"

if __name__ == "__main__":
    test_pcr_prediction = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction"
    test_bp_removed = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed"
    test_gene_annotation = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
    test_contig_trunc = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_contig_trunc"
    test_valid_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_valid_dir"
    test_11168 = "/home/sfisher/Sequences/11168_complete_genome"
    test_11168_cases = "/home/sfisher/Sequences/11168_test_files/gnomes_for_shannah"
    main(test_11168_cases, forward_primers, reverse_primers, amplicon_sequences)


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