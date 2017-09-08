import cProfile
import copy
import os
import subprocess
from collections import defaultdict

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq

from cgf_pred.Blastn import Blastn
from cgf_pred.HSP import HSP

#BLASTN LIMITATIONS
#db vs primer
MIN_BSR = 0.6         #for primers #TODO!!!
#db vs amp:
E_VALUE_CUTOFF = 0.01   #for ehyb
#other
QCOV_HSP_PERC = 70      #for primers and ehyb on same contig
WORD_SIZE = 7           #word size used for blastn
PERC_ID_CUTOFF = 80     #for both primers and ehyb

#PROGRAM SPECIFIC VALUES
MAX_MARGIN_BTWN_PRIMERS = 50    #max distance that can be +/- the correct distance between primers #TODO: switch to %?
CUTOFF_GENE_LENGTH = 70         #cutoff gene length for ehyb
SNP_THRESHOLD = 4               #5 bp must be an exact match with 3' end of primer
MAX_AMP_PRIMER_ALIGN = 10       #The amount of bp's that can be different btwn the hsp primer and hsp amp when searching for ehyb on same or diff contigs
MAX_PERC_EHYB_PRIMER_ENDS = 0.05
MAX_PERC_END = 0.05             #Max percentage of amp_length that will consider a primer at end of contig.

def blastn_query(query_genes, database, qcov, id):
    """ Outputs stdout from a blastn query using eval, qcov if specified, perc iden, and wordsize in xml format.

    :param query_genes: A path to a Fasta file w/ query genes
    :param database: A path to a fasta file that contains the database that is being searched against
    :param out_file: A xml file that the blastn query
    :restrictions: Database is formatted using makeblastdb
    :return: stdout in xml format
    """
    if qcov == True:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5,
                                             evalue=E_VALUE_CUTOFF, perc_identity=id, qcov_hsp_perc=QCOV_HSP_PERC)

    else:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=database, word_size=WORD_SIZE, outfmt=5,
                                             evalue=E_VALUE_CUTOFF, perc_identity=id)
    stdout, stderr = blastn_cline()
    return stdout


def bs_blast(query_genes, db):
    """ Outputs stdout from blastn in xml format using only perc id and word size

    :param query_genes: A path to a Fasta file w/ query genes
    :param db: A path to a Fasta file that contains the db that is being searched against.
    :return: stdout in xml format
    """
    blastn_cline = NcbiblastnCommandline(query=query_genes, db=db, word_size=WORD_SIZE, outfmt=5, perc_identity=PERC_ID_CUTOFF, qcov_hsp_perc=QCOV_HSP_PERC) #, task='blastn-short')
    stdout, stderr = blastn_cline()
    return stdout

def max_bs_blast(seqn_dir):
    """ Outputs stdout from blastn in xml format using 100% id

    :param seqn_dir: A path to a Fasta file containing sequences to search for bits against itself
    :return: stdout in xml format
    """
    blastn_cline = NcbiblastnCommandline(query=seqn_dir, db=seqn_dir, word_size=WORD_SIZE, outfmt=5, perc_identity=100, qcov_hsp_perc=100) #, task='blastn-short')
    stdout, stderr = blastn_cline()
    return stdout

def create_blastn_object(query_genes:str, database:str, qcov=False,id=PERC_ID_CUTOFF) -> Blastn:
    """ Return a blastn object with initialized blast_records and hsp_records

    :param query_genes: A path to a Fasta File containing query genes.
    :param database: A path to a Fasta file containing db being searched.
    :param qcov: An optional input to include qcov in blastn results.
    :restrictions: Database is formatted using makeblastdb
    :return: Blastn object
    """
    blastn_object = Blastn()
    stdout_xml = blastn_query(query_genes, database, qcov, id)
    blastn_object.create_blast_records(stdout_xml)
    blastn_object.create_hsp_objects(query_genes)
    return blastn_object

def create_blastn_bsr_object(query_genes, db):
    """ Returns a Blastn object with initialized blast_records and hsp_records with cutoff bsr.

    :param query_genes: A path to a fasta file containing all query genes
    :param db: A path to a fasta file containing a single database
    :return: A blast object with initialized blast_records and hsp_records with cutoff bsr
    """
    blastn_object = Blastn()
    stdout_xml = bs_blast(query_genes, db)
    blastn_object.create_blast_records(stdout_xml)
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

def create_dict_from_fasta_seqs(primers):
    """ Return a dictionary that contains primer sequences: key is the primer id and value is the primer sequence.

    :param primers: A fasta file containing primer sequences
    :return: A dictionary
    """
    primer_dict = {}
    for primer in SeqIO.parse(primers, "fasta"):
        primer_dict[primer.id] = primer.seq
    return primer_dict

def is_distance(f_hsp_object, r_hsp_object, amplicon_sequences):
    """ Modifies objects attributes to True if f & r are within MAX_MARGIN_BTWN_PRIMERS to each other.

    :param f_hsp_object: A HSP object
    :param r_hsp_object: A HSP object
    :param amplicon_sequences: A path to a Fasta file that contains amplicon sequences
    :restrictions: f_hsp_object and r_hsp_object must be on the same contig
    :return: None
    """
    assert f_hsp_object.contig_name == r_hsp_object.contig_name
    distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        if f_hsp_object.name in amplicon.id:
            amplicon_length = len(amplicon.seq)
            f_hsp_object.pcr_distance = (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)
            r_hsp_object.pcr_distance = (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)

def epcr(forward_hsp_object:HSP, reverse_hsp_object:HSP, amplicon_sequences, f_primers:dict, r_primers:dict):
    """ Modifies objects attributes to indicate valid strands (f & r hsp facing each other), if SNP located on 3' end
    within SNP_THRESHOLD and correct distance from each other.

    :param forward_hsp_object: A HSP Object
    :param reverse_hsp_object: A HSP Object
    :param amplicon_sequences: A path to a Fasta file that contains amplicon sequences
    :param f_primers: A dict containing f_primer seq'ns
    :param r_primers: A dictionary containing r_primer seq'ns
    :return: None
    """

    valid_strands(forward_hsp_object, reverse_hsp_object)
    is_snp_primer_search(forward_hsp_object, f_primers, r_primers)
    is_snp_primer_search(reverse_hsp_object, f_primers, r_primers)
    is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences)

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


def is_snp_primer_search(hsp_object: HSP, f_primers:dict, r_primers:dict):
    """ Modifies hsp object to indicate any SNP within SNP_THRESHOLD of 3' end

    :param hsp_object: a HSP object
    :param f_primers: A dictionary containing f_primer seq'ns
    :param r_primers: A dictionary containing r_primer seq'ns
    :return:None
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

def valid_dir(hsp: HSP):
    """ Modifies hsp object to determine if facing end of contig and within MAX_PERC_END of the amp len from the end of the contig.

    :param hsp: A HSP object
    :return: None
    """
    if not (abs((hsp.start + hsp.amp_len) - hsp.db_length - 1) <= (MAX_PERC_END * hsp.amp_len) and abs(hsp.start - hsp.amp_len) <= (MAX_PERC_END * hsp.amp_len)):
        if hsp.strand and abs((hsp.start + hsp.amp_len) - hsp.db_length - 1) <= (MAX_PERC_END * hsp.amp_len):
            hsp.location = True
            hsp.end_dist = abs((hsp.start + hsp.amp_len) - hsp.db_length - 1)
        elif not hsp.strand and abs(hsp.start - hsp.amp_len) <= (MAX_PERC_END * hsp.amp_len):
            hsp.location = True
            hsp.end_dist = abs(hsp.end)
        elif hsp.strand:
            hsp.location = False
            hsp.end_dist = abs((hsp.start + hsp.amp_len) - hsp.db_length - 1)
        elif not hsp.strand:
            hsp.location = False
            hsp.end_dist = abs(hsp.start - hsp.amp_len)

    else: #seq'n found over entire amp
        hsp.location = True
        hsp.end_dist = abs((hsp.start + hsp.amp_len) - hsp.db_length - 1)

def ehyb(blast_object: Blastn):
    """ Determines if hsp objects in blast_object are of length >= CUTOFF_GENE_LENGTH

    :param blast_object: a Blastn object containing hsp objects
    :return: A list of HSP's with ehybrid attribute initialized
    """
    lo_ehybrid_hsp = [hsp for hsp in blast_object.hsp_objects if hsp.length >= CUTOFF_GENE_LENGTH]
    lo_failures = [hsp for hsp in blast_object.hsp_objects if hsp.length < CUTOFF_GENE_LENGTH]

    for hsp in lo_ehybrid_hsp:
        hsp.ehybrid = True
    for hsp in lo_failures:
        hsp.ehybrid = False
    lo_ehybrid_results = lo_ehybrid_hsp + lo_failures
    return lo_ehybrid_results

def false_neg_pred(attr_bin_tree:list, hsp:HSP, i:int) -> int:
    """

    :param attr_bin_tree:
    :param hsp:
    :param i:
    :return:
    """
    attr_i = attr_bin_tree[i]
    if attr_i == None:
        assert i != None
        return i
    # elif attr_i == 'ehybrid':
    #     if getattr(hsp, attr_i) == None:
    #         print('amp vs db for', hsp.name, 'not found')
    #change tree to not contain attribute with amp_found!!!
    elif getattr(hsp, attr_i) == None:
        return i
    elif getattr(hsp, attr_i):
        return false_neg_pred(attr_bin_tree, hsp, (2 * i) + 1)
    elif not getattr(hsp, attr_i):
        return false_neg_pred(attr_bin_tree, hsp, (2 * i) + 2)

def max_bs(file_dir:str) -> dict:
    """ Determines the max bs of each genome in file_dir by running it against itself.

    :param file_dir: Path to fasta file containing reference genomes (query) needed for BSR
    :return: A dictionary containing the bits of each genome in file_dir
    """

    dict_bits = {}
    files = (file for file in os.listdir(file_dir) if file.endswith('.fasta'))
    files_paths = []
    for file in files:
        files_paths.append(os.path.abspath(file_dir) + '/' + file)
    for file_path in files_paths:
        stdout_xml = max_bs_blast(file_path)
        blastn_obj = Blastn()
        blastn_obj.create_blast_records(stdout_xml)
        blastn_obj.create_hsp_objects(file_path)
        assert len(blastn_obj.hsp_objects) == 1
        for hsp in blastn_obj.hsp_objects:
            print(hsp.name)
            if '11168_' in hsp.name:
                dict_bits[hsp.name[6:]] = hsp.bits
            else:
                dict_bits[hsp.name] = hsp.bits
    return dict_bits

def bsr(blast_object:Blastn, max_bits_dict:dict):
    """ Removes hsp object from blast_object if BSR is not >= MIN_BSR

    :param blast_object: A Blastn Object
    :param max_bits_dict: A dictionary containing the max_bs of query seq'ns
    :return: None
    """

    for hsp in blast_object.hsp_objects:
        if '11168_' in hsp.name:
            hsp.bsr = hsp.bits / max_bits_dict[hsp.name[6:]]
        else:
            hsp.bsr = hsp.bits / max_bits_dict[hsp.name]
        if hsp.bsr < MIN_BSR:
            # assert hsp in blast_object.hsp_objects

            blast_object.remove_hsp_object_all(hsp)
            assert hsp not in blast_object.hsp_objects

# def sga_ehybrid(lo_hsp_ehybrid):
#     """ Determines the genes in lo_hsp_ehybrid that are specified to only use ehyb
#
#     :param lo_hsp_ehybrid: A list of hsp that were found using ehybrid only
#     :return: dict containing hsp genes from lo_hsp_ehybrid that are specified to use ehyb only.
#     """
#
#     result_dict = defaultdict(list)
#     # ehyb_only = [] #['cj0421c', 'cj0246c', 'cj1294', 'cj1324', 'cj1427c', 'cj1721', 'cj1439', 'cj1552', 'cj1551'] #rm 0755, 'cj1334', 'cj1329'
#     for hsp in lo_hsp_ehybrid:
#         # for gene in ehyb_only:
#         #     if gene in hsp.name:
#                 result_dict[hsp.name].append(hsp)
#     return result_dict
#
# #NOT used
# def sga_epcr(lo_hsp_epcr):
#     """ Determines the genes in lo_hsp_epcr that are specified to only use epcr
#
#     :param lo_hsp_epcr: A list of genes that were found using epcr only
#     :return: dict containing hsp genes from lo_hsp_epcr that are specified to use epcr only.
#     """
#
#     result_dict = defaultdict(list)
#     for hsp in lo_hsp_epcr:
#         # this is where I will determine the genes that only require epcr
#         #if "gene" in hsp.name:
#         result_dict[hsp.name].append(hsp)
#     return result_dict
#
# def sga(hsp_pass_ehyb:list, result_dict:dict):
#     """ Adds gene hsps to result dictionary that are found using the specific gene approach.
#
#     :param hsp_pass_ehyb: A list of gene hsp that were found using ehybrid only
#     :param result_dict: A dict containing results so far (to add results using ehyb only to)
#     :return: result dict with added results with only ehyb
#     """
#
#     ehyb_dict = sga_ehybrid(hsp_pass_ehyb)
#     for name, lo_hsp in ehyb_dict.items():
#         if name not in result_dict:
#             #found using ehybrid
#             if "11168_" in name:
#                 result_dict[name[6:]].append(lo_hsp)
#             else:
#                 result_dict[name].append(lo_hsp)
#     return result_dict

def contig_copy(f_hsp_old, r_hsp_old, max_f_bits_dict, max_r_bits_dict, contig) -> list:
    """ Renders a copy of f_hsp_old and r_hsp_old and assigns attributes for when both primers are found.

    :param f_hsp_old: A HSP object that will be copied and is the partner of r_hsp_old
    :param r_hsp_old: A HSP object that will be copied and is the partner of f_hsp_old
    :param max_f_bits_dict: dict containing max bits for each of the forward primers
    :param max_r_bits_dict: dict containing max bits for each of the reverse primers
    :param contig: Boolean, True if f_hsp_old and r_hsp_old on same contig and False if on different contigs
    :return: None
    """

    f_hsp = copy.copy(f_hsp_old)
    r_hsp = copy.copy(r_hsp_old)
    f_hsp.partner = r_hsp
    r_hsp.partner = f_hsp
    f_hsp.both_primers_found = True
    r_hsp.both_primers_found = True
    f_hsp.contig = contig
    r_hsp.contig = contig
    f_hsp.bsr = f_hsp.bits / max_f_bits_dict[f_hsp.name]
    r_hsp.bsr = r_hsp.bits / max_r_bits_dict[r_hsp.name]
    return [f_hsp, r_hsp]

def ehyb_both_prim_found(blast, f_hsp, r_hsp):
    """ eHybridization when both primers are found.

    :param blast: Blast object using amp as query
    :param f_hsp: A HSP object
    :param r_hsp: A HSP object
    :return: None
    """

    lo_hsp_ehybrid_qcov = ehyb(blast) # assigns ehybrid attributes to each hsp from amp vs db
    ehybrid_qcov_pass = [hsp for hsp in lo_hsp_ehybrid_qcov if hsp.ehybrid == True]
    ehybrid_qcov_fail = [hsp for hsp in lo_hsp_ehybrid_qcov if hsp.ehybrid == False]

    for hsp in ehybrid_qcov_pass:
        # if f_hsp.name in hsp.name and r_hsp.name == hsp.name:
        if abs(f_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) or abs(r_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(f_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) or abs(r_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                and f_hsp.contig_name == hsp.contig_name:
            f_hsp.ehybrid, r_hsp.ehybrid = True, True
            f_hsp.amp_len, r_hsp.amp_len = hsp.length, hsp.length
            f_hsp.amp_sbjct, r_hsp.amp_sbjct = hsp.sbjct, hsp.sbjct
            f_hsp.amp_query, r_hsp.amp_query = hsp.query, hsp.query
    for hsp in ehybrid_qcov_fail:
        # if f_hsp.name in hsp.name and r_hsp.name in hsp.name:
        if abs(f_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) or abs(r_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(f_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) or abs(r_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                and r_hsp.contig_name == hsp.contig_name:
            f_hsp.ehybrid, r_hsp.ehybrid = False, False
            f_hsp.amp_len, r_hsp.amp_len = hsp.length, hsp.length
            f_hsp.amp_sbjct, r_hsp.amp_sbjct = hsp.sbjct, hsp.sbjct
            f_hsp.amp_query, r_hsp.amp_query = hsp.query, hsp.query

def same_contig_pred(lo_tup_same_contig, full_blast_qcov, dict_f_primers, dict_r_primers, max_f_bits_dict, max_r_bits_dict, amp_seq, debug):
    """ Prediction of ecgf when primers are located on the same contig.

    :param lo_tup_same_contig: List of tuples of hsp primers found on the same contig
    :param full_blast_qcov: Blast object using amp as query and qcov included in restriction
    :param dict_f_primers: Dictionary containing forward primers
    :param dict_r_primers: dictionary containing reverse primers
    :param max_f_bits_dict: dict containing max bits for each of the forward primers
    :param max_r_bits_dict: dict containing max bits for each of the reverse primers
    :return: list of results with a dict of ecgf results and a dict of all hsp's found if using debug option.
    """

    if debug == True:
        all_hsp = defaultdict(list)
    result_dict = defaultdict(list)

    for tup in lo_tup_same_contig:
        f_hsp_old = tup[0]
        r_hsp_old = tup[1]

        copied_o = contig_copy(f_hsp_old, r_hsp_old, max_f_bits_dict, max_r_bits_dict, True)
        f_hsp, r_hsp = copied_o[0], copied_o[1]

        epcr(f_hsp, r_hsp, amp_seq, dict_f_primers, dict_r_primers) #assigns valid and snp attributes and pcr_distance and epcr

        ehyb_both_prim_found(full_blast_qcov, f_hsp, r_hsp)

        if f_hsp.epcr and (f_hsp.ehybrid or r_hsp.ehybrid):
            result_dict[f_hsp.name].append(f_hsp)
            result_dict[r_hsp.name].append(r_hsp)

        if debug == True:
            if f_hsp.bsr >= MIN_BSR:
                all_hsp[f_hsp.name].append(f_hsp)
            if r_hsp.bsr >= MIN_BSR:
                all_hsp[r_hsp.name].append(r_hsp)

    if debug == True:
        return [result_dict, all_hsp]
    else:
        return [result_dict]

def diff_contig_pred(lo_tup_diff_contig, max_f_bits_dict, max_r_bits_dict, ehybrid_hsp_pass, ehybrid_hsp_fail, debug):
    """ Prediction of cgf when primers are located on different contigs.

    :param lo_tup_diff_contig: List of tuples containing f and r primers from the same query but on diff contigs
    :param max_f_bits_dict: dict containing max bits for each of the forward primers
    :param max_r_bits_dict: dict containing max bits for each of the reverse primers
    :param ehybrid_hsp_pass: A list of gene hsp that were found using ehybrid only
    :param ehybrid_hsp_fail: A list of gene hsp that were not found using ehybrid only b/c not long enough length
    :param debug: optional debug bool input
    :return: list of results with a dict of ecgf results and a dict of all hsp's found if using debug option.
    """
    if debug == True:
        all_hsp = defaultdict(list)
    result_dict = defaultdict(list)

    for tup in lo_tup_diff_contig:
        f_hsp_old = tup[0]
        r_hsp_old = tup[1]
        copied_o = contig_copy(f_hsp_old, r_hsp_old, max_f_bits_dict, max_r_bits_dict, False)
        f_hsp, r_hsp = copied_o[0], copied_o[1]

        for hsp in ehybrid_hsp_pass:
        # consider depending on strand !!! ie. which direction the primer is facing!!!
            if f_hsp.strand and not r_hsp.strand:
                if f_hsp.start == hsp.start or f_hsp.end == hsp.end and f_hsp.contig_name == hsp.contig_name:
                    f_hsp.ehybrid = True
                    f_hsp.amp_len = hsp.length
                    valid_dir(f_hsp)
                    f_hsp.amp_sbjct = hsp.sbjct
                    f_hsp.amp_query = hsp.query
                if r_hsp.start == hsp.end or r_hsp.end == hsp.start and r_hsp.contig_name == hsp.contig_name and not r_hsp.strand:
                    r_hsp.ehybrid = True
                    r_hsp.amp_len = hsp.length
                    valid_dir(r_hsp)
                    r_hsp.amp_sbjct = hsp.sbjct
                    r_hsp.amp_query = hsp.query
            else:
                if r_hsp.start == hsp.start or r_hsp.end == hsp.end and r_hsp.contig_name == hsp.contig_name and r_hsp.strand:
                    r_hsp.ehybrid = True
                    r_hsp.amp_len = hsp.length
                    valid_dir(r_hsp)
                    r_hsp.amp_sbjct = hsp.sbjct
                    r_hsp.amp_query = hsp.query
                if f_hsp.start == hsp.end or f_hsp.end == hsp.start and f_hsp.contig_name == hsp.contig_name and not f_hsp.strand:
                    f_hsp.ehybrid = True
                    f_hsp.amp_len = hsp.length
                    valid_dir(r_hsp)
                    f_hsp.amp_sbjct = hsp.sbjct
                    f_hsp.amp_query = hsp.query


            # assert len(result_dict[f_hsp.name]) == 0
        if (f_hsp.ehybrid and f_hsp.location) and (r_hsp.ehybrid and r_hsp.location):
            result_dict[f_hsp.name].append(f_hsp)
            result_dict[r_hsp.name].append(r_hsp)
            # print('added ',f_hsp.name, ' b/c on different contigs')
            # if (f_hsp.ehybrid and r_hsp.ehybrid) and (f_hsp.location and r_hsp.location):
            #     result_dict[f_hsp.name].append((f_hsp, r_hsp))
            # elif f_hsp.ehybrid and f_hsp.location:
            #     result_dict[f_hsp.name].append((f_hsp, None))
            # elif r_hsp.ehybrid and r_hsp.location:
            #     result_dict[r_hsp.name].append((None, r_hsp))
                # assert len(result_dict[f_hsp.name]) < 3
        for hsp in ehybrid_hsp_fail:
            if abs(f_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) or abs(f_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                    and f_hsp.contig_name == hsp.contig_name:
                f_hsp.ehybrid = False
                f_hsp.amp_len = hsp.length
                f_hsp.amp_query = hsp.query
                f_hsp.amp_sbjct = hsp.sbjct
            if abs(r_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) or abs(r_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                    and r_hsp.contig_name == hsp.contig_name:
                r_hsp.ehybrid = False
                r_hsp.amp_len = hsp.length
                r_hsp.amp_sbjct = hsp.sbjct
                r_hsp.amp_query = hsp.query
        if debug == True:
            if f_hsp.bsr >= MIN_BSR:
                all_hsp[f_hsp.name].append(f_hsp)
            if r_hsp.bsr >= MIN_BSR:
                all_hsp[r_hsp.name].append(r_hsp)

    if debug == True:
        return [result_dict, all_hsp]
    else:
        return [result_dict]

def single_primer_found(lo_hsp_single_primers, ehybrid_hsp_pass, ehybrid_hsp_fail, debug):
    """ Prediction of cgf when only one primer is found.

    :param lo_hsp_single_primers: list of hsp's where only one primer is found.
    :param ehybrid_hsp_pass: A list of gene hsp that were found using ehybrid only
    :param ehybrid_hsp_fail: A list of gene hsp that were not found using ehybrid only b/c not long enough length
    :param debug: optional debug boolean input
    :return: list of results with a dict of ecgf results and a dict of all hsp's found if using debug option.
    """
    result_dict = defaultdict(list)

    if debug == True:
        all_hsp  = defaultdict(list)
    for single_hsp in lo_hsp_single_primers:
        single_hsp.both_primers_found = False
        single_hsp.contig = False

        for blast_hsp in ehybrid_hsp_pass:
            if abs(single_hsp.start - blast_hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) or abs(single_hsp.end - blast_hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) and single_hsp.contig_name == blast_hsp.contig_name:
                # f_hsp.amp_found = True
                single_hsp.ehybrid = True
                single_hsp.amp_len = blast_hsp.length
                single_hsp.amp_sbjct = blast_hsp.sbjct
                single_hsp.amp_query = blast_hsp.query
                valid_dir(single_hsp)
                if single_hsp.location == True:
                    result_dict[single_hsp.name].append(single_hsp)
                    # print(single_hsp.name, 'added b/c of single primer found')
                    assert len(result_dict[single_hsp.name]) < 2
        for blast_hsp in ehybrid_hsp_fail:
            if abs(single_hsp.start - blast_hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) or abs(single_hsp.end - blast_hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) and single_hsp.contig_name == blast_hsp.contig_name:
                # f_hsp.amp_found = True
                single_hsp.ehybrid = False
                single_hsp.amp_len = blast_hsp.length
                single_hsp.amp_sbjct = blast_hsp.sbjct
                single_hsp.amp_query = blast_hsp.query

        if debug == True:
            if single_hsp.bsr >= MIN_BSR:
                all_hsp[single_hsp.name].append(single_hsp)
    if debug == True:
        return [result_dict, all_hsp]
    else:
        return [result_dict]

#TODO: for multiprocessing !!!!!
# output = mp.Queue()

def ecgf(forward_primers:str, reverse_primers:str, database:str, amp_sequences:str, max_f_bits_dict:dict, max_r_bits_dict:dict, max_amp_bits_dict:dict, debug=False) -> list:
    """ Predicts in vitro cgf (eCGF)

    :param forward_primers: A string containing the path to all the 40 primer sequences in a fasta file
    :param reverse_primers: A string containing the path to all the 40 reverse primer sequences in a fasta file
    :param database: A string containing the path to a directory that contains strains of interest in fasta files (formatted using makeblastdb)
    :param amp_sequences: A string containing the path to all the 40 amplicon sequences in a fasta file
    :param max_f_bits_dict: dict containing max bits for each of the forward primers from the forward_primers file.
    :param max_r_bits_dict: dict containing max bits for each of the reverse primers from the reverse_primers file.
    :param max_amp_bits_dict: dict containing max bits for each of the amplicon sequences in the amp_sequences file.
    :param debug: An debug option that returns all hits that failed in the program.
    :restrictions: fasta files in database must be formatted using makeblastdb.
    :return: debug=False: A list containing a dictionary of gene found (key = hsp.name, val = hsp)
             debug=True: A list containing a dictionary of genes found and all the hits that failed in the program.
    """

    # forward_blast = create_blastn_object(forward_primers, database, forward_out_file, True)
    # reverse_blast = create_blastn_object(reverse_primers, database, reverse_out_file, True)
    #TODO: multiprocessing when creating blastn objects?
    forward_blast_bsr = create_blastn_bsr_object(forward_primers, database)
    reverse_blast_bsr = create_blastn_bsr_object(reverse_primers, database)
    blast_object = create_blastn_object(amp_sequences, database)
    full_blast_qcov = create_blastn_object(amp_sequences, database, True)
    ehyb_blast = create_blastn_object(amp_sequences, database, False, 70)

    bsr(forward_blast_bsr, max_f_bits_dict)
    bsr(reverse_blast_bsr, max_r_bits_dict)
    dict_f_primers = create_dict_from_fasta_seqs(forward_primers)
    dict_r_primers = create_dict_from_fasta_seqs(reverse_primers)
    dict_amp = create_dict_from_fasta_seqs(amp_sequences)

    lo_tup_same_queries = [(f_hsp,r_hsp) for f_hsp in forward_blast_bsr.hsp_objects for r_hsp in reverse_blast_bsr.hsp_objects if f_hsp.name == r_hsp.name]
    lo_tup_same_contig = [tup for tup in lo_tup_same_queries if tup[0].contig_name == tup[1].contig_name]
    # lo_tup_diff_contig = [tup for tup in lo_tup_same_queries if tup not in lo_tup_same_contig]

    #same contig prediction
    same_contig = same_contig_pred(lo_tup_same_contig, full_blast_qcov, dict_f_primers, dict_r_primers, max_f_bits_dict, max_r_bits_dict, amp_sequences, debug)
    results_dict_same_contig = same_contig[0]
    all_hsp_same_contig = {}
    if debug == True:
        all_hsp_same_contig = same_contig[1]

    # print('same contig results', results_dict_same_contig)

    #doesn't look for genes on different contigs if they were already found on same contig
    lo_tup_diff_contig = [tup for tup in lo_tup_same_queries if tup not in lo_tup_same_contig and tup[0].name not in results_dict_same_contig]
    # lo_tup_diff_contig_not_found = [tup for tup in lo_tup_diff_contig if tup[0].name not in same_contig_result_keys]

    #ehybrid results (qcov not required)
    lo_hsp_ehybrid = ehyb(blast_object)  # assigns ehybrid attributes to each hsp from amp vs db
    ehybrid_hsp_pass = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == True]
    ehybrid_hsp_fail = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == False]

    #Different contig prediction
    diff_contig = diff_contig_pred(lo_tup_diff_contig, max_r_bits_dict, max_r_bits_dict, ehybrid_hsp_pass, ehybrid_hsp_fail, debug)
    results_diff_contig = diff_contig[0]
    all_hsp_diff_contig = {}
    if debug == True:
        all_hsp_diff_contig = diff_contig[1]

    # f_hsp_single_primers = [hsp for hsp in forward_blast_bsr.hsp_objects if hsp not in lo_f_primers]
    diff_contig_results_keys = [key for key, val in results_diff_contig.items()]
    same_contig_result_keys = [key for key,value in results_dict_same_contig.items()]

    f_hsp_single_primers = [hsp for hsp in forward_blast_bsr.hsp_objects if hsp.name not in same_contig_result_keys and hsp.name not in diff_contig_results_keys]
    r_hsp_single_primers = [hsp for hsp in reverse_blast_bsr.hsp_objects if hsp.name not in same_contig_result_keys and hsp.name not in diff_contig_results_keys]

    #testing
    f_hsp_test = [hsp.name for hsp in f_hsp_single_primers]
    r_hsp_test = [hsp.name for hsp in r_hsp_single_primers]

    #assigns bsr !!!
    lo_hsp_single_primers = []
    for f_hsp in f_hsp_single_primers:
        f_hsp.bsr = f_hsp.bits / max_f_bits_dict[f_hsp.name]
        if f_hsp.bsr >= MIN_BSR:
            lo_hsp_single_primers.append(f_hsp)
    for r_hsp in r_hsp_single_primers:
        r_hsp.bsr = r_hsp.bits / max_r_bits_dict[r_hsp.name]
        if r_hsp.bsr >= MIN_BSR:
            lo_hsp_single_primers.append(r_hsp)
        # assert r_hsp.bsr > MIN_BSR
    # lo_hsp_single_primers = list(itertools.chain(f_hsp_single_primers, r_hsp_single_primers))
    for hsp in lo_hsp_single_primers:
        assert hsp.bsr >= MIN_BSR

    #CHANGED TO LOOK FOR ALL HSP'S THAT WERE NOT FOUND IN same contig or diff contig results!!!
    #One primer found prediction
    one_primer = single_primer_found(lo_hsp_single_primers, ehybrid_hsp_pass, ehybrid_hsp_fail, debug)
    results_one_primer = one_primer[0]
    all_hsp_one_primer = defaultdict(list)
    if debug == True:
        all_hsp_one_primer = one_primer[1]

    #combine all results
    result_dict = defaultdict(list)
    all_hsp = defaultdict(list)
    if debug == True:
        for k,v in all_hsp_same_contig.items():
            for hsp in v:
                all_hsp[k].append(hsp)
        for k,v in all_hsp_diff_contig.items():
            for hsp in v:
                all_hsp[k].append(hsp)
        for k,v in all_hsp_one_primer.items():
            for hsp in v:
                all_hsp[k].append(hsp)
    for k,v in results_dict_same_contig.items():
        result_dict[k].append(v)
    for k,v in results_diff_contig.items():
        result_dict[k].append(v)
    for k,v in results_one_primer.items():
        result_dict[k].append(v)

    #CASE 4:
    lo_hsp_ehybrid = ehyb(ehyb_blast)  # assigns ehybrid attributes to each hsp from amp vs db
    ehybrid_pass = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == True]
    ehyb_pos = [hsp for hsp in ehybrid_pass if hsp.name[6:] not in result_dict.keys()]
    fourth_case_check(ehyb_pos, result_dict, dict_f_primers, dict_r_primers, dict_amp)
    # result_dict = sga(ehybrid_qcov_pass, result_dict) #add genes that are found using sga if not already found.

    if debug == True:
        results_list = [result_dict, all_hsp]
    else:
        #TODO: remove ehyb_pos from here... I added it for testing.
        results_list = [result_dict, ehyb_pos]

    # print('epcr only results', epcr_only_results)
    #TODO !!! comment out for multiprocessing
    return results_list
    # output.put(result_dict)

def fourth_case_check(ehyb_pos, result_dict, f_primer_dict, r_primer_dict, amp_dict):
    file = open("/home/sfisher/Sequences/11168_test_files/fourth_case_check/8_per_genome_trial", "a")
    ehyb_pos_names = [hsp.name[6:] for hsp in ehyb_pos]
    file.write('\n \n Genes not found using eCGF but found using eHYB ' + str(ehyb_pos_names))
    file.write('\n Number of genes found in eHYB only ' + str(len(ehyb_pos)))
    file.write('\n Number of genes found using eCGF (Probability of being TRULY +ve is ~high) '+ str(len(result_dict.keys())))
    file.write('\n Number of genes not found at all (Probability of being TRULY -ve is high) '+ str(abs(40 - len(result_dict.keys()) - len(ehyb_pos))))
    file.write('\n Below is more info on the genes found using eHYB only (was -ve in eCGF but +ve in eHYB): ')
    for hsp in ehyb_pos:
        hsp_name = hsp.name[6:]
        perc_id = hsp.identities / hsp.length
        file.write('\n \n' + hsp_name)
        file.write('\n % id ' + str(perc_id))
        file.write('\n qcov (% of query seqn that overlaps sbjct seqn) ' + str(len(hsp.query) / (len(amp_dict[hsp.name]) + 1)))
        file.write('\n Entire Query Seq      ' + str(len(amp_dict[hsp.name]) + 1) + " 1 " + str(len(amp_dict[hsp.name])) + " " + str(amp_dict[hsp.name]))
        file.write('\n Matching Query Seq    ' + str(len(hsp.query)) + " " + str(hsp.query_start) + " " + str(hsp.query_end) + " "*(hsp.query_start + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.query)
        file.write('\n Matching Sbjct Seq    ' + str(len(hsp.sbjct)) + " " + str(hsp.query_start) + " " + str(hsp.query_end) + " "*(hsp.query_start + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.sbjct)
        file.write('\n Forward primer         ' + str(len(f_primer_dict[hsp_name])) + "       " + str(f_primer_dict[hsp_name]))
        file.write('\n Reverse primer         ' + str(len(r_primer_dict[hsp_name])) + " "*(len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]) + 7) + str(amp_dict[hsp.name][len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]): ] ))

        # if (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
        file.write('\n Forward Primer qcov: ' +  str((len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name])) + " (CUTOFF: " + str((QCOV_HSP_PERC / 100)) + " )")
        # if (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
        file.write('\n Reverse Primer qcov: ' + str((len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name])) + " (CUTOFF:  " + str((QCOV_HSP_PERC / 100)) + " )")
        # if perc_id * 100 < PERC_ID_CUTOFF:
        file.write('\n Perc ID for ehyb in eCGF: '+ str(perc_id) + " (CUTOFF: " + str((PERC_ID_CUTOFF / 100)) + " )")
    file.close()

    # ehyb_pos_names = [hsp.name[6:] for hsp in ehyb_pos]
    # print('\nGenes not found using eCGF but found using eHYB', ehyb_pos_names)
    # print('Number of genes found in eHYB only', len(ehyb_pos))
    # print('Number of genes found using eCGF (Probability of being TRULY +ve is ~high)', len(result_dict.keys()))
    # print('Number of genes not found at all (Probability of being TRULY -ve is high)', abs(40 - len(result_dict.keys()) - len(ehyb_pos)))
    # print('Below is more info on the genes found using eHYB only (was -ve in eCGF but +ve in eHYB): ')
    # for hsp in ehyb_pos:
    #     hsp_name = hsp.name[6:]
    #     perc_id = hsp.identities / hsp.length
    #     print('\n', hsp_name)
    #     print('% id', perc_id)
    #     print('qcov (% of query seqn that overlaps sbjct seqn)', str(len(hsp.query) / (len(amp_dict[hsp.name]) + 1)))
    #     print('Entire Query Seq      ', len(amp_dict[hsp.name]) + 1, "1 ", len(amp_dict[hsp.name]), amp_dict[hsp.name])
    #     print('Matching Query Seq    ', len(hsp.query) , hsp.query_start, hsp.query_end, " "*(hsp.query_start-1 + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))), hsp.query)
    #     print('Matching Sbjct Seq    ', len(hsp.sbjct) , hsp.query_start, hsp.query_end, " "*(hsp.query_start-1 + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))), hsp.sbjct)
    #     print('Forward primer        ', len(f_primer_dict[hsp_name]) , "       ", f_primer_dict[hsp_name])
    #     print('Reverse primer        ', len(r_primer_dict[hsp_name]), " "*(len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]) + 7), amp_dict[hsp.name][len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]): ] )
    #
    #     if (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
    #         print('Forward Primer qcov: ', (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name]), " < ", (QCOV_HSP_PERC / 100))
    #     if (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
    #         print('Reverse Primer qcov: ', (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]), " < ", (QCOV_HSP_PERC / 100))
    #     if perc_id * 100 < PERC_ID_CUTOFF:
    #         print('Perc ID for ehyb in eCGF:', perc_id, " < ", (PERC_ID_CUTOFF / 100))

def per_gene_fourth_case_check(all_ehyb_pos, forward_primers, reverse_primers, amplicon_sequences):

    f_primer_dict = create_dict_from_fasta_seqs(forward_primers)
    r_primer_dict = create_dict_from_fasta_seqs(reverse_primers)
    amp_dict = create_dict_from_fasta_seqs(amplicon_sequences)
    gene_ehyb_pos_dict = defaultdict(list)

    for ehyb_pos in all_ehyb_pos.values():
        for hsp in ehyb_pos:
            gene_ehyb_pos_dict[hsp.name].append(hsp)

    file = open("/home/sfisher/Sequences/11168_test_files/fourth_case_check/8_per_gene_trial", "a")
    for gene_name, lo_hsp in gene_ehyb_pos_dict.items():
        file.write("\n \n  " + gene_name + " not found using eCGF but found using eHYB " + str(len(lo_hsp)) + " times.")
        for hsp in lo_hsp:
            hsp_name = hsp.name[6:]
            perc_id = hsp.identities / hsp.length
            file.write('\n \n' + hsp_name)
            file.write('\n % id ' + str(perc_id))
            file.write('\n qcov (% of query seqn that overlaps sbjct seqn) ' + str(len(hsp.query) / (len(amp_dict[hsp.name]) + 1)))
            file.write('\n Entire Query Seq      ' + str(len(amp_dict[hsp.name]) + 1) + " 1 " + str(len(amp_dict[hsp.name])) + " " + str(amp_dict[hsp.name]))
            file.write('\n Matching Query Seq    ' + str(len(hsp.query)) + " " + str(hsp.query_start) + " " + str(hsp.query_end) + " "*(hsp.query_start + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.query)
            file.write('\n Matching Sbjct Seq    ' + str(len(hsp.sbjct)) + " " + str(hsp.query_start) + " " + str(hsp.query_end) + " "*(hsp.query_start + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.sbjct)
            file.write('\n Forward primer         ' + str(len(f_primer_dict[hsp_name])) + "       " + str(f_primer_dict[hsp_name]))
            file.write('\n Reverse primer         ' + str(len(r_primer_dict[hsp_name])) + " "*(len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]) + 7) + str(amp_dict[hsp.name][len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]): ] ))

            # if (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
            file.write('\n Forward Primer qcov: ' + str((len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name])) + " (CUTOFF: " + str((QCOV_HSP_PERC / 100)) + " )")
            # if (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
            file.write('\n Reverse Primer qcov: ' + str((len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name])) + " (CUTOFF:  " + str((QCOV_HSP_PERC / 100)) + " )")
            # if perc_id * 100 < PERC_ID_CUTOFF:
            file.write('\n Perc ID for ehyb in eCGF: ' + str(perc_id) + " (CUTOFF: " + str((PERC_ID_CUTOFF / 100)) + " )")

    file.close()



#TODO: changed db_dir to file_path when calling on each genome.
def main(db_fasta, f_primers_fasta, r_primers_fasta, amp_fasta):
    """

    :param db_directory: The location of the directory with fasta database files contained (already formatted using makeblastdb)
    :param f_primers_fasta: The fasta file location of the forward primers
    :param r_primers_fasta: The fasta file location of the reverse primers
    :param amplicon_sequences: The fasta file location of the amplicon sequences
    :return: A dictionary
    """
    #makeblastdb all inputs if not already done
    blastdb_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb'.format(db_fasta)
    blastfp_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb_fp'.format(f_primers_fasta)
    blastrp_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb_rp'.format(r_primers_fasta)
    blastamp_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb_amp'.format(amp_fasta)
    DB_process = subprocess.Popen(blastdb_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    DB_process.wait()
    DB_process = subprocess.Popen(blastfp_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    DB_process.wait()
    DB_process = subprocess.Popen(blastrp_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    DB_process.wait()
    DB_process = subprocess.Popen(blastamp_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    DB_process.wait()

    #bsr
    f_bs_primer_dir = "/home/sfisher/Sequences/BSR/f_primers/"
    r_bs_primer_dir = "/home/sfisher/Sequences/BSR/r_primers/"
    amp_bs_dir = "/home/sfisher/Sequences/BSR/amp_seq/"
    max_f_bits_dict = max_bs(f_bs_primer_dir)
    max_r_bits_dict = max_bs(r_bs_primer_dir)
    max_amp_bits_dict = max_bs(amp_bs_dir)

    files = (file for file in os.listdir(db_fasta) if file.endswith('.fasta'))
    files_paths = []
    for file in files:
        files_paths.append(os.path.abspath(db_fasta) + '/' + file)
    cgf_predictions_dict = {}

    #TODO: added dict for testing
    all_ehyb_pos = {}

    #TODO: comment out for multiprocessing!
    for file_path in files_paths:
        file_name = file_path.partition(db_fasta + "/")[2]
        result = ecgf(f_primers_fasta, r_primers_fasta, file_path, amp_fasta, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)
        cgf_predictions_dict[file_name] = result[0]
        #TODO: delete below... added for testing:
        all_ehyb_pos[file_name] = result[1]

    per_gene_fourth_case_check(all_ehyb_pos, f_primers_fasta, r_primers_fasta, amp_fasta)
    # multiprocessing
    # processes = [mp.Process(target=ecgf, args=(f_primers_fasta, r_primers_fasta, file_path, amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)) for file_path in files_paths]
    # for p in processes:
    #     p.start()
    # for p in processes:
    #     p.join()
    # results = [output.get() for p in processes]
    # print('results!!!', results)

    #multiprocessing using pool
    # pool = mp.Pool(processes=4)
    # results = [pool.apply(ecgf, args=(f_primers_fasta, r_primers_fasta, file_path, amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)) for file_path in files_paths]
    # print('results', results)

    print(cgf_predictions_dict)
    return cgf_predictions_dict


# TODO: Commented out so I could use as package in github. (running main from __main__.py
# if __name__ == "__main__":
#     # test_pcr_prediction = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_pcr_prediction"
#     # test_bp_removed = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_bp_removed"
#     # test_gene_annotation = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_gene_annotation_error"
#     # test_contig_trunc = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_contig_trunc"
#     # test_valid_dir = "/home/sfisher/Sequences/amplicon_sequences/individual_amp_seq/test_valid_dir"
#     # test_11168 = "/home/sfisher/Sequences/11168_complete_genome"
#     forward_primers = "/home/sfisher/Sequences/cgf_forward_primers.fasta" #fasta file with primer id's and primer sequences
#     reverse_primers = "/home/sfisher/Sequences/cgf_reverse_primers.fasta" #fasta file with primer id's and primer sequences
#     amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
#     all_genomes = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
#     debug_cases = "/home/sfisher/Sequences/11168_test_files/debug_genes"
#     multipro_tests = "/home/sfisher/Sequences/11168_test_files/multipro_tests"
#     # test_memory = "/home/sfisher/Sequences/11168_test_files/memory_trial"
#     # test_memory2 = "/home/sfisher/Sequences/11168_test_files/memory_trial2"
#     # test_cprofile_file = "/home/sfisher/Sequences/11168_test_files/memory_trial/06_2855.fasta"
#
#     #TODO: check for error for if makeblastdb hasn't been run and return an error message to the user indicating the error.
#
#     main(all_genomes, forward_primers, reverse_primers, amplicon_sequences)
#     # cProfile.run('cgf_prediction_trial(forward_primers, reverse_primers, test_cprofile_file, amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)')
#     # cProfile.run('main(all_genomes, forward_primers, reverse_primers, amplicon_sequences); print')
#     # main(test_11168_cases, forward_primers, reverse_primers, amplicon_sequences)



#TODO: With database input, makeblastdb or run through command line?