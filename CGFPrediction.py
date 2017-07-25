from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Blastn import Blastn
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from HSP import HSP
import os
import errno
import itertools
import subprocess
from collections import Counter
from collections import defaultdict

SCORE=135
E_VALUE_CUTOFF = 1.0
PERC_ID_CUTOFF = 70
QCOV_HSP_PERC = 70
WORD_SIZE = 7
MAX_MARGIN_BTWN_PRIMERS = 50 #<= #TODO: change this value to something smaller...
CUTOFF_GENE_LENGTH = 70
SNP_THRESHOLD = 5 #5 bp must be an exact match with 3' end of primer
MAX_MM_CONTIG_END = 10 # The amount of bp's that can be located on end/start of db primer before reaching the end/start of amp
MAX_MARGIN_AMP = 0

def blastn_query(query_genes, database, out_file, qcov):
    """ Outputs a blastn query into an xml file.

    :param query_genes: A fasta file that contains the query genes that are being searched in the database
    :param database: A fasta file that contains the database that is being searched against
    :param out_file: A xml file that the blastn query
    :restrictions: Database is formatted using makeblastdb
    :return: None
    """
    #eval, perc iden, qcov
    if qcov == True:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5, out=out_file,
                                             evalue=E_VALUE_CUTOFF, perc_identity=PERC_ID_CUTOFF, qcov_hsp_perc=QCOV_HSP_PERC)

    #blast score
    # if qcov == True:
    #     blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5, out=out_file)

    else:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5,
                                             out=out_file, evalue=E_VALUE_CUTOFF, perc_identity=PERC_ID_CUTOFF)
    stdout, stderr = blastn_cline()

def create_blastn_object(query_genes, database, out_file, qcov=False):
    """ Return a blastn object with initialized blast_records and hsp_records

    :param query_genes: A fasta file that contains the query genes that are being searched in the database
    :param database: A fasta file that contains the database that is being searched against
    :param out_file: A string that contains the file location of the xml file being outputted
    :restrictions: Database is formatted using makeblastdb
    :return: Blastn object
    """
    blastn_object = Blastn()
    blastn_query(query_genes, database, out_file, qcov)
    blastn_object.create_blast_records(out_file)
    blastn_object.create_hsp_objects(query_genes)
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

def is_distance_trial(f_hsp_object, r_hsp_object, amplicon_sequences):
    assert f_hsp_object.contig_name == r_hsp_object.contig_name
    distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        if f_hsp_object.name in amplicon.id:
            amplicon_length = len(amplicon.seq)
            f_hsp_object.pcr_distance = (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)
            r_hsp_object.pcr_distance = (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)

def is_distance(f_hsp_object, r_hsp_object, amplicon_sequences):
    """ Return True if f_hsp_object and r_hsp_object are within MAX_MARGIN_BTWN_PRIMERS to each other.

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
            return (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)

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

    #testing
    snp_list = []
    if forward_hsp_object.snp == True:
        snp_list.append(forward_hsp_object.name)
        snp_list.append(forward_hsp_object.contig_name)
    if reverse_hsp_object.snp == True:
        snp_list.append(reverse_hsp_object.name)
        snp_list.append(reverse_hsp_object.contig_name)
    # print('snp list', snp_list)

    assert forward_hsp_object.valid == reverse_hsp_object.valid #should always be True because valid_strands assigns both objects the same valid attribute
    if forward_hsp_object.valid == True:
        if (forward_hsp_object.snp == False or reverse_hsp_object.snp == False):
            pcr_distance = is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences)
            #Result
            forward_hsp_object.pcr_distance = pcr_distance
            reverse_hsp_object.pcr_distance = pcr_distance

            return pcr_distance
        elif forward_hsp_object.snp == True and reverse_hsp_object.snp == True:
            return False
    elif forward_hsp_object.valid == False:
        return False
    else:
        assert True == False #valid not initialized

def pcr_directly_trial(forward_hsp_object:HSP, reverse_hsp_object:HSP, amplicon_sequences, f_primers:dict, r_primers:dict) -> bool:

    valid_strands(forward_hsp_object, reverse_hsp_object)
    is_snp_primer_search(forward_hsp_object, f_primers, r_primers)
    is_snp_primer_search(reverse_hsp_object, f_primers, r_primers)
    is_distance_trial(forward_hsp_object, reverse_hsp_object, amplicon_sequences)

    assert forward_hsp_object.valid == reverse_hsp_object.valid
    if forward_hsp_object.valid == True:
        if (forward_hsp_object.snp == False or reverse_hsp_object.snp == False):
            assert forward_hsp_object.pcr_distance == reverse_hsp_object.pcr_distance
            forward_hsp_object.epcr = reverse_hsp_object.pcr_distance
            reverse_hsp_object.epcr = reverse_hsp_object.pcr_distance
        else:
            forward_hsp_object.epcr, reverse_hsp_object.epcr = False, False
    else:
        forward_hsp_object.epcr, reverse_hsp_object.epcr = False, False


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

        # print('forward primer', forward_primer_seq)
        # print('forward primer rcomp', forward_primer_reverse_complement)
        # print('reverse primer', reverse_primer_seq)
        # print('reverse priemr rcomp', reverse_primer_reverse_complement)
        # print('db end forward', db_end_forward)
        # print('db end reverse', db_end_reverse)

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

def valid_dir_trial(hsp: HSP):
    if not (abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END and hsp.start <= MAX_MM_CONTIG_END):
        if hsp.strand == True and abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END:
            hsp.location = True
        elif hsp.strand == False and abs(hsp.end) <= MAX_MM_CONTIG_END:
            hsp.location = True
        else:
            hsp.location = False
    else:
        hsp.location = True

def ehybridization_trial(blast_object: Blastn):
    lo_ehybrid_hsp = [hsp for hsp in blast_object.hsp_objects if hsp.length >= CUTOFF_GENE_LENGTH]
    lo_failures = [hsp for hsp in blast_object.hsp_objects if hsp not in lo_ehybrid_hsp]
    for hsp in lo_ehybrid_hsp:
        hsp.ehybrid = True
    for hsp in lo_failures:
        hsp.ehybrid = False
    lo_ehybrid_results = itertools.chain(lo_ehybrid_hsp, lo_failures)
    return lo_ehybrid_results


def cgf_prediction_trial(forward_primers:str, reverse_primers:str, database:str, forward_out_file:str, reverse_out_file:str, amplicon_sequences:str, full_out_file:str) -> list:

    #TODO: make sure that after adding a hsp to the results list that if the hsp is modified in later fcns, it is not affected.

    forward_blast = create_blastn_object(forward_primers, database, forward_out_file, True)
    reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file, True)
    blast_object = create_blastn_object(amplicon_sequences, database, full_out_file)
    full_blast_qcov = create_blastn_object(amplicon_sequences, database, full_out_file, True)
    dict_f_primers = create_primer_dict(forward_primers)
    dict_r_primers = create_primer_dict(reverse_primers)

    lo_tup_same_queries = [(f_hsp, r_hsp) for f_hsp in forward_blast.hsp_objects for r_hsp in reverse_blast.hsp_objects if f_hsp.name == r_hsp.name]
    lo_tup_same_contig = [tup for tup in lo_tup_same_queries if tup[0].contig_name == tup[1].contig_name]
    lo_tup_diff_contig = [tup for tup in lo_tup_same_queries if tup not in lo_tup_same_contig]
    try:
        lo_f_primers, lo_r_primers = zip(*lo_tup_same_queries)
    except ValueError:
        lo_f_primers = []
        lo_r_primers = []
    f_hsp_single_primers = (hsp for hsp in forward_blast.hsp_objects if hsp not in lo_f_primers)
    r_hsp_single_primers = (hsp for hsp in reverse_blast.hsp_objects if hsp not in lo_r_primers)
    lo_hsp_single_primers = list(itertools.chain(f_hsp_single_primers, r_hsp_single_primers))

    result_dict = defaultdict(list)
    all_hsp = defaultdict(list)

    #same contig
    lo_hsp_ehybrid_qcov = ehybridization_trial(full_blast_qcov) # assigns ehybrid attributes to each hsp from amp vs db
    ehybrid_qcov_hsp = [hsp for hsp in lo_hsp_ehybrid_qcov if hsp.ehybrid == True]
    ehybrid_qcov_fail = [hsp for hsp in lo_hsp_ehybrid_qcov if hsp.ehybrid == False]
    # print('ehybrid qcov', ehybrid_qcov_hsp)
    # print('ehybrid fail qcov', ehybrid_qcov_fail)
    # print('set!!!', set(full_blast_qcov.hsp_objects) - set(ehybrid_qcov_hsp) - set(ehybrid_qcov_failed_hsp) == set())
    #diff contigs and one found
    lo_hsp_ehybrid = ehybridization_trial(blast_object)  # assigns ehybrid attributes to each hsp from amp vs db
    ehybrid_hsp = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == True]
    ehybrid_hsp_fail = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == False]
    # print('ehybrid', ehybrid_hsp)
    # ehybrid_failed_hsp = [hsp for hsp in blast_object.hsp_objects if hsp.ehybrid == False]
    # print('set2!!!', set(blast_object.hsp_objects) - set(ehybrid_hsp) - set(ehybrid_failed_hsp) == set())

    #Same contig prediction
    lo_tup_same_contig_results = []
    for tup in lo_tup_same_contig:
        f_hsp = tup[0]
        r_hsp = tup[1]

        f_hsp.both_primers_found = True
        r_hsp.both_primers_found = True
        f_hsp.contig = True
        r_hsp.contig = True

        pcr_directly_trial(f_hsp, r_hsp, amplicon_sequences, dict_f_primers, dict_r_primers) #assigns valid and snp attributes and pcr_distance and epcr
        assert f_hsp.epcr == r_hsp.epcr
        # if f_hsp.epcr == True and r_hsp.epcr == True:
        #     lo_tup_same_contig_pcr_results.append(tup)

        for hsp in ehybrid_qcov_hsp:
            if f_hsp.start == hsp.start or r_hsp.start == hsp.start or f_hsp.end == hsp.end or r_hsp.end == hsp.end:
                f_hsp.ehybrid, r_hsp_ehybrid = True, True
                f_hsp.amp_len, r_hsp.amp_len = hsp.length, hsp.length
        for hsp in ehybrid_qcov_fail:
            if f_hsp.start == hsp.start or r_hsp.start == hsp.start or f_hsp.end == hsp.end or r_hsp.end == hsp.end:
                f_hsp.ehybrid, r_hsp.ehybrid = False, False
                f_hsp.amp_len, r_hsp.amp_len = hsp.length, hsp.length
        # assert f_hsp.ehybrid == r_hsp.ehybrid

        if f_hsp.ehybrid and f_hsp.epcr:
            # assert result_dict[f_hsp.name] == None
            result_dict[f_hsp.name].append(f_hsp)
            result_dict[r_hsp.name].append(r_hsp)
            # assert len(result_dict[f_hsp.name]) == 2

        all_hsp[f_hsp.name].append(f_hsp)
        all_hsp[r_hsp.name].append(r_hsp)


    #Different contig prediction
    lo_tup_diff_contig_results = []
    for tup in lo_tup_diff_contig:
        f_hsp = tup[0]
        r_hsp = tup[1]

        f_hsp.both_primers_found = True
        r_hsp.both_primers_found = True
        f_hsp.contig = False
        r_hsp.contig = False

        #todo: make into helper fcn
        for hsp in ehybrid_hsp:
            if f_hsp.start == hsp.start or f_hsp.end == hsp.end:
                f_hsp.ehybrid = True
                f_hsp.amp_len = hsp.length
                valid_dir_trial(f_hsp)
            if r_hsp.start == hsp.start or r_hsp.end == hsp.end:
                r_hsp.ehybrid = True
                r_hsp.amp_len = hsp.length
                valid_dir_trial(r_hsp)
            # assert len(result_dict[f_hsp.name]) == 0
            if (f_hsp.ehybrid and r_hsp.ehybrid) and (f_hsp.location and r_hsp.location):
                result_dict[f_hsp.name].append((f_hsp, r_hsp))
            elif f_hsp.ehybrid and f_hsp.location:
                result_dict[f_hsp.name].append((f_hsp, None))
            elif r_hsp.ehybrid and r_hsp.location:
                result_dict[r_hsp.name].append((None, r_hsp))
            # assert len(result_dict[f_hsp.name]) < 3
        for hsp in ehybrid_hsp_fail:
            if f_hsp.start == hsp.start or f_hsp.end == hsp.end:
                f_hsp.ehybrid = False
                f_hsp.amp_len = hsp.length
            if r_hsp.start == hsp.start or r_hsp.end == hsp.end:
                r_hsp.ehybrid = False
                r_hsp.amp_len = hsp.length

        all_hsp[f_hsp.name].append(f_hsp)
        all_hsp[r_hsp.name].append(r_hsp)

    #One primer found prediction
    lo_tup_one_primer_results = []
    for single_hsp in lo_hsp_single_primers:

        single_hsp.both_primers_found = False
        single_hsp.contig = False

        for blast_hsp in ehybrid_hsp:
            if single_hsp.start == blast_hsp.start or single_hsp.end == blast_hsp.end:
                single_hsp.ehybrid = True
                single_hsp.amp_len = blast_hsp.length
                valid_dir_trial(single_hsp)
                if single_hsp.location == True:
                    assert result_dict[single_hsp.name] == None
                    result_dict[single_hsp.name].append(single_hsp)
                    assert len(result_dict[single_hsp.name]) < 2
        for blast_hsp in ehybrid_hsp_fail:
            if single_hsp.start == blast_hsp.start or single_hsp.end == blast_hsp.end:
                single_hsp.ehybrid = False
                single_hsp.amp_len = blast_hsp.length

        all_hsp[single_hsp.name].append(single_hsp)


    results_list = [result_dict, all_hsp]

    return results_list

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

    #TODO: use subprocess to makeblastdb using command line.
    # subprocess.check_call(['makeblastdb', '-in', '/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests/06-1505.fasta', '-dbtype', 'nucl'])


    #create new folder for out files
    try:
        os.mkdir(db_directory + "/out_files")
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    cgf_predictions_dict = {}
    for file_path in files_paths:
        # print('FILE!!!', file_path)
        file_name = file_path.partition(db_directory + "/")[2]
        # print(file_name)
        f_out_file_path = db_directory + "/out_files/" + "f_" + file_name.replace("fasta", "xml")
        r_out_file_path = db_directory + "/out_files/" + "r_" + file_name.replace("fasta", "xml")
        full_out_file_path = db_directory + "/out_files/" + "full_" + file_name.replace("fasta", "xml")

        result = cgf_prediction_trial(forward_primers, reverse_primers, file_path, f_out_file_path, r_out_file_path, amplicon_sequences, full_out_file_path)[0]
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






# #TODO: make for entire gene?
# #TODO: delete? ... Currently isn't used!!!
# def is_snp(hsp_object, f_primers, r_primers):
#     #TODO:!!!
#     """Return True if a mismatch or removed bp is located in the first SNP_THRESHOLD bp's on the 3' end
#     :param hsp_object: A HSP Object that is facing in the right direction (towards contig)
#     :param f_primers: A dictionary containing full forward primer sequences.
#     :param r_primers: A dictionary containing full reverse primer sequences.
#     :return: A boolean: True if mismatch or removed bp on 3' end, False otherwise
#     """
#
#     #TODO: Test that this will not work if the primer sequences do not correspond to the hsp_object
#     #TODO: "_" only in entire_gene calls?
#     #TODO: this only works when there is no deletions in the 5' end of the primer!
#
#     if "_" in hsp_object.name: #handles when entire gene called with amp seq names
#         assert("_" in hsp_object.name)
#         hsp_object.set_name(hsp_object.name.partition('_')[2])
#     if hsp_object.name in f_primers and hsp_object.name in r_primers:
#         f_primer = f_primers[hsp_object.name]
#         r_primer = r_primers[hsp_object.name]
#         forward_primer_seq = Seq(str(f_primer[len(f_primer) - SNP_THRESHOLD : len(f_primer)]), generic_dna)
#         forward_primer_reverse_complement = forward_primer_seq.reverse_complement()
#         reverse_primer_seq = Seq(str(r_primer[len(r_primer) - SNP_THRESHOLD : len(r_primer)]), generic_dna)
#         reverse_primer_reverse_complement = reverse_primer_seq.reverse_complement()
#
#         starting_forward = hsp_object.sbjct[:len(f_primer)]  # [len(f_primer) - SNP_THRESHOLD : len(f_primer)]
#         starting_reverse = hsp_object.sbjct[:len(r_primer)] # [len(r_primer) - SNP_THRESHOLD : len(r_primer)]
#         ending_forward = hsp_object.sbjct[len(hsp_object.sbjct) - len(f_primer) : len(hsp_object.sbjct) - len(f_primer) + SNP_THRESHOLD]
#         ending_reverse = hsp_object.sbjct[len(hsp_object.sbjct) - len(r_primer) : len(hsp_object.sbjct) - len(r_primer) + SNP_THRESHOLD]
#
#         # print(hsp_object.sbjct)
#         # print(forward_primer_seq)
#         # print(forward_primer_reverse_complement)
#         # print(reverse_primer_seq)
#         # print(reverse_primer_reverse_complement)
#         # print('normal forward', starting_forward)
#         # print('comp forward', ending_forward)
#         # print('normal reverse', ending_reverse)
#         # print('reverse comp reverse', starting_reverse)
#
#         #does not consider insertions into the primer seq! If need to consider, make threshold value for # of insertions allowed ie. [:len(primer) + threshold]
#
#         if forward_primer_seq or forward_primer_reverse_complement in starting_forward \
#                 or forward_primer_reverse_complement or forward_primer_seq in ending_forward\
#                 and hsp_object.strand == True:
#             hsp_object.snp = False
#             #TODO: hsp_object.primer = "forward" !!!?
#         elif str(reverse_primer_seq) in ending_reverse \
#                 or str(reverse_primer_reverse_complement) in starting_reverse \
#                 or str(reverse_primer_reverse_complement) in ending_reverse\
#                 or str(reverse_primer_seq) in starting_reverse:
#             hsp_object.snp = False
#         else:
#             hsp_object.snp = True



# def valid_dir(lo_hsp:list) -> list:
#     """ Modify lo_hsps to only contain primer sequences that are facing and located at the end of the contig
#
#     :param lo_hsps: list of hsp objects
#     :return: None
#     """
#
#     lo_hsp_valid_dir = []
#     for hsp in lo_hsp:
#         #considers the case where the entire strand is found on the contig
#         if not (abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END and hsp.start <= MAX_MM_CONTIG_END):
#             # the the hsp is facing the end of the contig
#             # print('hsp', hsp.name)
#             # print('hsp contig', hsp.contig_name)
#             # print('hsp strand', hsp.strand)
#             # print('identities', hsp.identities)
#             # print('query seq', hsp.query)
#             # print('sbjct seq', hsp.sbjct)
#             # print('gaps', hsp.gaps)
#             # print('alignment', hsp.alignment)
#             # print('first condition... strand = true', abs(hsp.end - len(hsp.sbjct)))
#             # print('second condition... strand = false', abs(hsp.end))
#             # print('start', hsp.start)
#             # print('end', hsp.end)
#             # print('contig end mm', abs(hsp.end - len(hsp.sbjct)))
#             # print('contig end reverse', hsp.end)
#
#             if hsp.strand == True and abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END:
#                 lo_hsp_valid_dir.append(hsp)
#             elif hsp.strand == False and abs(hsp.end) <= MAX_MM_CONTIG_END:
#                 lo_hsp_valid_dir.append(hsp)
#             # else:
#             #     if hsp.strand == True:
#             #         print("hsp was NOT found: start < end and end - len = ", abs(hsp.end - len(hsp.sbjct)))
#             #         print("for debugging, start - len =", abs(hsp.start - len(hsp.sbjct)))
#             #     else:
#             #         print("hsp was NOT found: start > end and start - len = ", abs(hsp.start - len(hsp.sbjct)))
#             #         print("for debugging, end - len =", abs(hsp.end - len(hsp.sbjct)))
#         else:
#             lo_hsp_valid_dir.append(hsp)
#
#     #testing
#     lo_valid_dir_names = []
#     for hsp in lo_hsp_valid_dir:
#         lo_valid_dir_names.append(hsp.name)
#     # print('lo valid dir', lo_valid_dir_names)
#
#     return lo_hsp_valid_dir

    # This is the case where valid_dir is called with a list of tuples! Save for future if needed.
    # # REMEMBER: if the hsp is on the lagging strand, end will be at the beginning of the strand
    # count = 0
    # lo_found = []
    # for tup in lo_hsp_tup:
    #     for hsp in tup:
    #         lo_tup_found = []
    #         if hsp != None:
    #             # considers the case where the entire strand is found
    #             if not (abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END and hsp.start <= MAX_MM_CONTIG_END):
    #                 # the the hsp is facing the end of the contig
    #                 print('start', hsp.start)
    #                 print('end', hsp.end)
    #                 print('contig end mm', abs(hsp.end - len(hsp.sbjct)))
    #                 print('contig end reverse', hsp.end)
    #
    #                 if hsp.strand == True and abs(hsp.end - len(hsp.sbjct)) <= MAX_MM_CONTIG_END:
    #                     lo_tup_found.append(hsp)
    #                 elif hsp.strand == False and abs(hsp.end) <= MAX_MM_CONTIG_END:
    #                     lo_tup_found.append(hsp)
    #                 else:
    #                     if hsp.strand == True:
    #                         print("hsp was NOT found: start < end and end - len = ", abs(hsp.end - len(hsp.sbjct)))
    #                         print("for debugging, start - len =", abs(hsp.start - len(hsp.sbjct)))
    #                     else:
    #                         print("hsp was NOT found: start > end and start - len = ", abs(hsp.start - len(hsp.sbjct)))
    #                         print("for debugging, end - len =", abs(hsp.end - hsp.sbjct))
    #             else:
    #                 lo_tup_found.append(hsp)
    #                 count += 1
    #                 print('the hsp extends along the entire contig! This should happen 40 times')
    #     if lo_tup_found != None:
    #         lo_found.append(tuple(lo_tup_found))

    # print ('count', count)
    # return lo_found

    # for hsp in lo_hsps:
    #     # considers the case where the entire strand is found
    #     if (hsp.end == hsp.query_end or hsp.start == hsp.query_start) and not (hsp.end == hsp.query_end and hsp.start == hsp.query_start):
    #         #Need to consider bp's that are removed here!!!
    #         if hsp.strand == False and abs(hsp.query_end - hsp.end) <= MAX_MM_CONTIG_END:
    #             lo_hsps.remove(hsp)
    #         elif hsp.strand == True and abs(hsp.start - hsp.query_start) <= MAX_MM_CONTIG_END:
    #             lo_hsps.remove(hsp)



# def ehybridization(blast_object:Blastn) -> list:
#     lo_hsp = [hsp for hsp in blast_object.hsp_objects if hsp.length >= CUTOFF_GENE_LENGTH]
#     #Results
#     lo_failures = [hsp for hsp in blast_object.hsp_objects if hsp not in lo_hsp]
#     for hsp in lo_failures:
#         hsp.ehybrid = False
#     lo_glob_hsp_rejects.extend(lo_failures)
#
#     return lo_hsp

# def same_contig_pred(blast_object:Blastn, lo_tup_same_contig:list, dict_f_primers:dict, dict_r_primers:dict):
#     """ Finds query in db through ehybridization followed by epcr.
#         WORKFLOW: ehybridization --> epcr
#
#     :param blast_object: A Blastn object.
#     :param lo_tup_same_contig: A list of tuples with forward and reverse hsp objects
#     :param dict_f_primers: A dict containing forward primers
#     :param dict_r_primers: A dict containing reverse primers
#     :return: A list of lists of forward and reverse HSP objects.
#     """
#
#     lo_queries = []
#
#     #added for testing!!!
#     # lo_blast = []
#     # for hsp in blast_object.hsp_objects:
#     #     lo_blast.append(hsp.name)
#     # print('lo blast', lo_blast)
#
#     lo_hybridized_hsp = ehybridization(blast_object)
#
#     #TODO: added for testing!!!
#     lo_hybridized_hsp_names = []
#     for hsp in lo_hybridized_hsp:
#         lo_hybridized_hsp_names.append(hsp.name)
#     # print('hybridized', lo_hybridized_hsp_names)
#
#     for hsp_tup in lo_tup_same_contig:
#         if any((hsp_tup[0].start == hsp.start or hsp_tup[1].start == hsp.start) for hsp in lo_hybridized_hsp):
#             #Results
#             hsp_tup[0].ehybrid = True
#             hsp_tup[0].ehybrid = True
#
#             if pcr_directly(hsp_tup[0], hsp_tup[1], amplicon_sequences, dict_f_primers, dict_r_primers):
#                 #Results
#                 hsp_tup[0].epcr = True
#                 hsp_tup[1].epcr = True
#
#                 lo_queries.append(list(hsp_tup))
#             else:
#                 #Results
#                 #wasn't found in epcr
#                 hsp_tup[0].epcr = False
#                 hsp_tup[1].epcr = False
#                 lo_glob_hsp_rejects.append(hsp_tup[0])
#                 lo_glob_hsp_rejects.append(hsp_tup[1])
#         elif any ((hsp_tup[0].end == hsp.end or hsp_tup[1].end == hsp.end) for hsp in lo_hybridized_hsp):
#             #Results
#             hsp_tup[0].ehybrid = True
#             hsp_tup[1].ehybrid = True
#
#             if pcr_directly(hsp_tup[0], hsp_tup[1], amplicon_sequences, dict_f_primers, dict_r_primers):
#                 #Results
#                 hsp_tup[0].epcr = True
#                 hsp_tup[1].epcr = True
#
#                 lo_queries.append(list(hsp_tup))
#             else:
#                 #Results
#                 #wasn't found in epcr
#                 hsp_tup[0].epcr = False
#                 hsp_tup[1].epcr = False
#                 lo_glob_hsp_rejects.append(hsp_tup[0])
#                 lo_glob_hsp_rejects.append(hsp_tup[1])
#
#
#    #TODO: added for testing!!! can delete
#     lo_hsp_names = []
#     for quer in lo_queries:
#         for hsp in quer:
#             lo_hsp_names.append(hsp.name)
#     # print('lo_queries', lo_hsp_names)
#     # print('diff', len(lo_hsp_names) - 2 * len(lo_hybridized_hsp_names))
#
#     #TODO: added hybridized results for testing purposes. delete from return value later
#     return (lo_queries, lo_hybridized_hsp_names)


# # TODO: consider case with cj1134 and cj1324!!?
# # TODO: This assumes that it would be highly unlikely that a gene will be found of long enough length on a database more than once!!!
# def diff_contig_pred(blast_object:Blastn, lo_tup_diff_contig:list, dict_f_primers:dict, dict_r_primers:dict):
#     """ Return the hsp objects (could be both forward and reverse or just f or r) that have the correct orientation and long enough length to be considered found.
#
#      Primers located on different contigs.
#         WORKFLOW: entire_gene -> epcr
#         1) both forward and reverse primers found (lo_tup_diff_contig) found
#         2) just forward or just reverse primer found
#         3)
#     """
#     if len(lo_tup_diff_contig) == 0:
#         return []
#
#     lo_hsp_correct_len = find_correct_len_hsps(blast_object, lo_tup_diff_contig)
#
#     #check the hsp's are facing toward the end of the contig and located near the end of the contig.
#     lo_correct_dir = valid_dir(lo_hsp_correct_len)
#     lo_correct_dir_names = [hsp.name for hsp in lo_correct_dir]
#
#     #only considers case when 1 or 2 hsps are of correct length and facing the end of contig.
#     correct_dir_counter = Counter(lo_correct_dir_names)
#     print('keys', correct_dir_counter.values())
#     for key in correct_dir_counter.values() :
#         assert key < 3
#
#     #create tuples out of lo_correct_dir
#     lo_tup_correct_dir = []
#     for hsp in lo_correct_dir:
#         for hsp2 in lo_correct_dir:
#             if hsp.name == hsp2.name:
#                 tup = (hsp, hsp2)
#                 lo_tup_correct_dir.append(tup)
#
#     #TODO: consider case where only one hsp is found of long enough length and facing end of contig!!!
#
#     #testing
#     print_lo_correct_dir_names = []
#     for hsp in lo_correct_dir:
#         print_lo_correct_dir_names.append(hsp.name)
#     print('lo correct dir and length', print_lo_correct_dir_names)
#
#     return lo_tup_correct_dir
#
#     # # #TODO: ensure primer is found too (not just the correct length)!... this will be considered my 'epcr' part of the algorithm!
#     # lo_tup_correct_len = []
#     # for tup in lo_tup_diff_contig:
#     #     if tup[0] == None and tup[1] != None:
#     #         if tup[1].length > CUTOFF_GENE_LENGTH:
#     #             lo_tup_correct_len.append(tup)
#     #     elif tup[0] != None and tup[1] == None:
#     #         if tup[0].length > CUTOFF_GENE_LENGTH:
#     #             lo_tup_correct_len.append(tup)
#     #     elif tup[0] != None and tup[1] != None:
#     #         print("tup contains both hsp's")
#     #         print('tup[0] length', tup[0].length)
#     #         print('tup[1] length', tup[1].length)
#     #         if tup[0].length > CUTOFF_GENE_LENGTH and tup[1].length >= CUTOFF_GENE_LENGTH:
#     #             lo_tup_correct_len.append(tup)
#     #         elif tup[0].length > CUTOFF_GENE_LENGTH:
#     #             one_tup_elm = (tup[0], None)
#     #             lo_tup_correct_len.append(one_tup_elm)
#     #         elif tup[1].length > CUTOFF_GENE_LENGTH:
#     #             # TODO: !!!????????????
#     #             one_tup_elm = (None, tup[1])
#     #             lo_tup_correct_len.append(one_tup_elm)


# def find_correct_len_hsps(blast_object:Blastn, lo_tup_diff_contig:list) -> list:
#     """
#
#     :param blast_object:
#     :param lo_tup_diff_contig:
#     :return:
#     """
#
#     hsp_objects = [hsp for hsp in blast_object.hsp_objects if hsp.length >= CUTOFF_GENE_LENGTH]
#     hsp_objects_iterating_through = []
#     for hsp in hsp_objects:
#         hsp_objects_iterating_through.append(hsp.name)
#         hsp_objects_iterating_through.append(hsp.length)
#
#     for tup in lo_tup_diff_contig:
#         assert len(tup) > 0
#
#     lo_hsp_diff_contigs_0, lo_hsp_diff_contigs_1 = zip(*lo_tup_diff_contig)
#     lo_hsp_diff_contigs_names = []
#     for hsp in lo_hsp_diff_contigs_0:
#         lo_hsp_diff_contigs_names.append(hsp.name)
#
#     lo_hsp_correct_len = []
#     for hsp in hsp_objects:
#         if hsp.name in lo_hsp_diff_contigs_names:
#             lo_hsp_correct_len.append(hsp)
#         elif hsp.name[6:] in lo_hsp_diff_contigs_names:
#             hsp.name = hsp.name[6:]
#             lo_hsp_correct_len.append(hsp)
#
#     return lo_hsp_correct_len

# def one_primer_pred(blast_object:Blastn, lo_hsp_one_primers: list):
#     """
#
#     :param blast_object:
#     :param lo_hsp_one_primers:
#     :return:
#     """
#
#     if len(lo_hsp_one_primers) == 0:
#         return []
#
#     #testing
#     # for hsp in blast_object.hsp_objects:
#     #     if hsp.length >= CUTOFF_GENE_LENGTH:
#     #         print('name and length', hsp.name, hsp.length)
#
#     lo_hsp_len = [hsp for hsp in blast_object.hsp_objects if hsp.length >= CUTOFF_GENE_LENGTH]
#     lo_hsp_one_primers_names = [hsp.name for hsp in lo_hsp_one_primers]
#     print('primer names', lo_hsp_one_primers_names)
#     lo_hsp_correct_len = []
#     for hsp in lo_hsp_len:
#         if hsp.name in lo_hsp_one_primers_names:
#             lo_hsp_correct_len.append(hsp)
#         elif hsp.name[6:] in lo_hsp_one_primers_names:
#             hsp.name = hsp.name[6:]
#             lo_hsp_correct_len.append(hsp)
#
#     #testing
#     # lo_hsp_correct_len_names = []
#     # for hsp in lo_hsp_correct_len:
#     #     lo_hsp_correct_len_names.append(hsp.name)
#     # print('hsp correct len names1', lo_hsp_correct_len_names)
#
#     #check the hsp's are facing toward the end of the contig.
#     lo_correct_dir = valid_dir(lo_hsp_correct_len)
#
#     #if the hsp's are in the middle of the contig, checks if they are of long enough length to be considered found.
#     # lo_incorrect_dir = [hsp for hsp in lo_hsp_correct_len if hsp not in lo_correct_dir]
#     # print(lo_incorrect_dir)
#     # lo_mid_hsps_found = []
#     # for hsp in lo_incorrect_dir:
#     #     print('length', hsp.length)
#     #     print('query length', hsp.query_length)
#     #     if abs(hsp.length - hsp.query_length) <= MAX_MARGIN_AMP:
#     #         if (hsp.identities / hsp.length) >= 0.94 and (hsp.identities / hsp.length) <= 1:
#     #             lo_mid_hsps_found.append(hsp)
#
#     #only considers case when 1 or 2 hsps are of correct length and facing the end of contig.
#     lo_correct_dir_names = [hsp.name for hsp in lo_correct_dir]
#     correct_dir_counter = Counter(lo_correct_dir_names)
#     # print('keys', correct_dir_counter.values())
#     for key in correct_dir_counter.values() :
#         assert key == 1
#
#     # lo_all_hsp_found = list(itertools.chain(lo_correct_dir, lo_mid_hsps_found))
#
#     return lo_correct_dir



# #TODO: Change return value to hash?
# def cgf_prediction(forward_primers:str, reverse_primers:str, database:str, forward_out_file:str, reverse_out_file:str, amplicon_sequences:str, full_out_file:str) -> list:
#     """ Return a list of hsp's that would be predicted as hsp's in vitro...
#
#     :param forward_primers: A Fasta file location containing forward query primers (ids and sequences)
#     :param reverse_primers: A Fasta file location containing reverse query primers (ids and sequences)
#     :param database: A Fasta file path (formatted using makeblastdb) that contains a database (complete or draft)
#     :param forward_out_file: An xml file location to store blastn results from forward primers as queries
#     :param reverse_out_file: An xml file location to store blastn results from reverse primers as queries
#     :param amplicon_sequences: A Fasta file containing amplicon sequences
#     :param full_out_file: An xml file location to store blastn results from amplicon sequences as queries
#     :return: List of HSP's that were found in database.
#     """
#
#     forward_blast = create_blastn_object(forward_primers, database, forward_out_file, True)
#     reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file, True)
#     blast_object = create_blastn_object(amplicon_sequences, database, full_out_file)
#     full_blast_qcov = create_blastn_object(amplicon_sequences, database, full_out_file, True)
#
#     # #testing!!!
#     # forward_blast_names = []
#     # for hsp in forward_blast.hsp_objects:
#     #     forward_blast_names.append(hsp.name)
#     # reverse_blast_names = []
#     # for hsp in reverse_blast.hsp_objects:
#     #     reverse_blast_names.append(hsp.name)
#     # full_blast_names_lengths = []
#     # for hsp in full_blast_qcov.hsp_objects:
#     #     full_blast_names_lengths.append(hsp.name)
#     #     full_blast_names_lengths.append(hsp.length)
#     # blast_names_lengths = []
#     # for hsp in blast_object.hsp_objects:
#     #     blast_names_lengths.append(hsp.name)
#     #     blast_names_lengths.append(hsp.length)
#     # print('forward blast', forward_blast_names)
#     # print('reverse blast', reverse_blast_names)
#     # print('blast (no qcov)', blast_names_lengths)
#     # print('qcov blast', full_blast_names_lengths)
#
#     dict_f_primers = create_primer_dict(forward_primers)
#     dict_r_primers = create_primer_dict(reverse_primers)
#
#     lo_tup_same_queries = [(f_hsp, r_hsp) for f_hsp in forward_blast.hsp_objects for r_hsp in reverse_blast.hsp_objects if f_hsp.name == r_hsp.name]
#     lo_tup_same_contig = [tup for tup in lo_tup_same_queries if tup[0].contig_name == tup[1].contig_name]
#     lo_tup_diff_contig = [tup for tup in lo_tup_same_queries if tup not in lo_tup_same_contig]
#
#     try:
#         lo_f_primers, lo_r_primers = zip(*lo_tup_same_queries)
#     except ValueError:
#         lo_f_primers = []
#         lo_r_primers = []
#     f_hsp_single_primers = (hsp for hsp in forward_blast.hsp_objects if hsp not in lo_f_primers)
#     r_hsp_single_primers = (hsp for hsp in reverse_blast.hsp_objects if hsp not in lo_r_primers)
#     lo_hsp_single_primers = list(itertools.chain(f_hsp_single_primers, r_hsp_single_primers))
#
#     #Results
#     lo_f_same_contig, lo_r_same_contig = zip(*lo_tup_same_contig)
#     for hsp in lo_f_same_contig:
#         hsp.both_primers_found = True
#         hsp.contig = True
#     for hsp in lo_r_same_contig:
#         hsp.both_primers_found = True
#         hsp.contig = True
#     #Results
#     try:
#         lo_f_diff_contig, lo_r_diff_contig = zip(*lo_tup_diff_contig)
#     except ValueError:
#         lo_f_diff_contig = []
#         lo_r_diff_contig = []
#     for hsp in lo_f_diff_contig:
#         hsp.both_primers_found = True
#         hsp.contig = False
#     for hsp in lo_r_same_contig:
#         hsp.both_primers_found = True
#         hsp.contig = False
#     #Results
#     for hsp in lo_hsp_single_primers:
#         hsp.one_primer_found = True
#         hsp.contig = False
#     #Results(false -ve)
#     f_blast_rejects = [hsp for hsp in forward_blast.hsp_objects if hsp not in lo_f_same_contig and hsp not in lo_f_diff_contig and hsp not in f_hsp_single_primers]
#     r_blast_rejects = [hsp for hsp in reverse_blast.hsp_objects if hsp not in lo_r_same_contig and hsp not in lo_r_diff_contig and hsp not in r_hsp_single_primers]
#     print(f_blast_rejects)
#     print(r_blast_rejects)
#     # assert f_blast_rejects == []
#     # assert r_blast_rejects == []
#     lo_glob_hsp_rejects.extend(f_blast_rejects)
#     lo_glob_hsp_rejects.extend(r_blast_rejects)
#
#     #testing
#     lo_hsp_single_primers_names = [hsp.name for hsp in lo_hsp_single_primers]
#
#     #TODO: added for testing purposes. Can delete later.
#     lo_first_hsp_in_tup_diff_contigs = [hsp[0].name for hsp in lo_tup_diff_contig]
#     print('tuples on different contigs', database, lo_first_hsp_in_tup_diff_contigs)
#
#     #predicts when primers are located on the same contig.
#     lo_same_contig_queries_results = same_contig_pred(full_blast_qcov, lo_tup_same_contig, dict_f_primers, dict_r_primers)
#     #TODO returned tuples for testing... remove later
#     lo_same_contig_queries = lo_same_contig_queries_results[0] #list of queries
#     lo_hybridized_results = lo_same_contig_queries_results[1]
#
#     #predicts when primers are located on different contigs.
#     lo_tup_diff_contig_hsp_found = diff_contig_pred(blast_object, lo_tup_diff_contig, dict_f_primers, dict_r_primers)
#
#     #predicts when only one primer is found.
#     lo_one_primer_hsps_found = one_primer_pred(blast_object, lo_hsp_single_primers) #not tuples
#     # lo_one_primer_hsps_tup = []
#     # for hsp in lo_one_primer_hsps_found:
#     #     tup = (hsp, None)
#
#     lo_queries = list(itertools.chain(lo_same_contig_queries, lo_tup_diff_contig_hsp_found, [lo_one_primer_hsps_found]))
#
#     #TODO: added extra return values for testing purposes!!! delete later.
#     return_list = [lo_queries, lo_first_hsp_in_tup_diff_contigs, lo_hybridized_results, lo_hsp_single_primers_names, lo_same_contig_queries, lo_glob_hsp_rejects]    #testing
#
#     # #testing
#     # for query in lo_queries:
#     #     for hsp in query:
#     #         print(hsp.name, hsp.query)
#     #         print(hsp.name, hsp.sbjct)
#
#     return return_list



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