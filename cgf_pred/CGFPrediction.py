import cProfile
from memory_profiler import profile
import copy
import os
import subprocess
import gc
from collections import defaultdict

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq

from cgf_pred.Blastn import Blastn
from cgf_pred.HSP import HSP
from pathlib import Path
from math import log
import csv

#BLASTN LIMITATIONS
#db vs primer
MIN_BSR = 0.60         # for primers #TODO!!!
#db vs amp:
E_VALUE_CUTOFF = 0.01   # for ehyb
#other
QCOV_HSP_PERC = 70      # for primers and ehyb on same contig
WORD_SIZE = 7           # word size used for blastn
PERC_ID_CUTOFF = 80     # for both primers and ehyb

#PROGRAM SPECIFIC VALUES
MAX_MARGIN_BTWN_PRIMERS = 50    #max distance from correct distance between primers
CUTOFF_GENE_LENGTH = 70         #cutoff gene length for ehyb
SNP_THRESHOLD = 4               #SNP_THRESHOLD bp must be an exact match with 3' end of primer
MIN_SNP_HAM_DIST = 2
MAX_AMP_PRIMER_ALIGN = 10       #The amount of bp's that can be different btwn the hsp primer and hsp amp when searching for ehyb on same or diff contigs
MAX_PERC_EHYB_PRIMER_ENDS = 0.05
MAX_PERC_END = 0.10             #Max percentage of amp_length that will consider a primer at end of contig.
SINGLE_PRIMER_ID = 0.90

#bit score calculation
LAMBDA = 1.28
K = 0.460

GENE_LIST = ['11168_cj0008', '11168_cj0033', '11168_cj0035', '11168_cj0057', '11168_cj0177', '11168_cj0181',
             '11168_cj0264c', '11168_cj0297c', '11168_cj0298c', '11168_cj0307', '11168_cj0421c', '11168_cj0483',
             '11168_cj0486', '11168_cj0566', '11168_cj0569', '11168_cj0570', '11168_cj0625', '11168_cj0728',
             '11168_cj0733', '11168_cj0736', '11168_cj0755', '11168_cj0860', '11168_cj0967', '11168_cj1134',
             '11168_cj1136', '11168_cj1141', '11168_cj1294', '11168_cj1324', '11168_cj1329', '11168_cj1334',
             '11168_cj1427c', '11168_cj1431c', '11168_cj1439', '11168_cj1550c', '11168_cj1551', '11168_cj1552',
             '11168_cj1585', '11168_cj1679', '11168_cj1721', '11168_cj1727c']

def blastn_query1(query_genes, db, qcov=False, evalue=False, id=PERC_ID_CUTOFF):
    """ Runs blastn with parameters specified and puts results into xml format

    :param query_genes: A fasta file containing all of the queries
    :param db: A fasta file containing the database to query search
    :param qcov: Optional qcov for Blastn
    :param evalue: Optional evalue for Blastn
    :param id: Optional percent identity cutoff for Blastn
    :return: stdout xml format
    """

    #TODO: ONLY RUN MAKEBLASTDB IF THE DB IS NOT ALREADY FORMATTED!!!
    blastdb_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb'.format(db)
    DB_process = subprocess.run(blastdb_cmd,
                                  shell=True,
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                check=True)
    # DB_process.wait()

    if qcov and evalue:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=db, word_size=WORD_SIZE, outfmt=5,
                                             perc_identity=id, qcov_hsp_perc=QCOV_HSP_PERC, evalue=E_VALUE_CUTOFF)
    elif qcov:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=db, word_size=WORD_SIZE, outfmt=5,
                                             perc_identity=id, qcov_hsp_perc=QCOV_HSP_PERC)
    elif evalue:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=db, word_size=WORD_SIZE, outfmt=5,
                                             perc_identity=id, evalue=E_VALUE_CUTOFF)
    else:
        blastn_cline = NcbiblastnCommandline(query=query_genes, db=db, word_size=WORD_SIZE, outfmt=5,
                                             perc_identity=id)
    stdout, stderr = blastn_cline()
    return stdout

def blastn_query_exceptions(query_gene, db, qcov):
    blastdb_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb'.format(db)
    DB_process = subprocess.run(blastdb_cmd,
                                  shell=True,
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                check=True)

    blastn_cline = NcbiblastnCommandline(query = query_gene, db = db, word_size = WORD_SIZE, outfmt = 5,
                                         perc_identity = PERC_ID_CUTOFF, qcov_hsp_perc = qcov)
    stdout, stderr = blastn_cline()
    return stdout

def create_blastn_object_exceptions(query_gene:str, db:str, qcov):
    blastn_object = Blastn()
    stdout_xml = blastn_query_exceptions(query_gene, db, qcov)
    blastn_object.create_blast_records(stdout_xml)
    blastn_object.create_hsp_objects(query_gene)
    return blastn_object


def create_blastn_object(query_genes:str, db:str, qcov=False,id=PERC_ID_CUTOFF):
    """ Return a blastn object with initialized blast_records and hsp_records

    :param query_genes: A path to a Fasta File containing query genes.
    :param database: A path to a Fasta file containing db being searched.
    :param qcov: An optional input to include qcov in blastn results.
    :restrictions: Database is formatted using makeblastdb
    :return: Blastn object
    """
    blastn_object = Blastn()
    stdout_xml = blastn_query1(query_genes, db, qcov=qcov, id=id) #TODO: CHANGED FROM blastn_query
    blastn_object.create_blast_records(stdout_xml)
    blastn_object.create_hsp_objects(query_genes)
    return blastn_object.hsp_objects

def create_blastn_bsr_object(query_genes, db):
    """ Returns a Blastn object with initialized blast_records and hsp_records with cutoff bsr.

    :param query_genes: A path to a fasta file containing all query genes
    :param db: A path to a fasta file containing a single database
    :return: A blast object with initialized blast_records and hsp_records with cutoff bsr
    """
    blastn_object = Blastn()
    stdout_xml = blastn_query1(query_genes, db, qcov=True) #TODO: CHANGED FROM bs_blast
    blastn_object.create_blast_records(stdout_xml)
    blastn_object.create_hsp_objects(query_genes)
    return blastn_object

def valid_strands(first_hsp_object: HSP, second_hsp_object: HSP) -> None :
    """ Assigns valid attributes to first and second hsp object.

    :param first_hsp_object: A HSP object to compare with second_hsp_object
    :param second_hsp_object: A HSP object to compare with first_hsp_object
    :return: None
    """

    if first_hsp_object.name == second_hsp_object.name:
        if (first_hsp_object.strand or second_hsp_object.strand) \
                and not (first_hsp_object.strand and second_hsp_object.strand):
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
        if '11168_' not in primer.id:
            primer_dict['11168_' + primer.id] = primer.seq
        else:
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
    # assert f_hsp_object.contig_name == r_hsp_object.contig_name
    distance = abs(f_hsp_object.start - r_hsp_object.start) + 1
    for amplicon in SeqIO.parse(amplicon_sequences, "fasta"):
        if f_hsp_object.name in amplicon.id:
            amplicon_length = len(amplicon.seq)
            f_hsp_object.pcr_distance = (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)
            r_hsp_object.pcr_distance = (abs(amplicon_length - distance) <= MAX_MARGIN_BTWN_PRIMERS)

def extend_sbjct(hsp_object, database, primer_dict):
    start = hsp_object.start
    end = hsp_object.end
    query_start = hsp_object.query_start
    query_end = hsp_object.query_end
    len_missing =  len(primer_dict[hsp_object.name]) - abs(end - start)
    begin_missing = query_start - 1
    end_missing = abs(len(primer_dict[hsp_object.name]) - query_end)
    if len_missing > 0:
        with open(database, 'r') as fasta:
            for contig in SeqIO.parse(fasta, 'fasta'):
                if contig.name == hsp_object.contig_name:
                    if start > end:
                        seq_found = contig[end - end_missing - 1: start + begin_missing]
                        hsp_object.sbjct = seq_found.reverse_complement().seq
                    else:
                        seq_found = contig[start - begin_missing - 1: end + end_missing]
                        hsp_object.sbjct = seq_found.seq

def epcr(forward_hsp_object:HSP, reverse_hsp_object:HSP, amplicon_sequences, f_primers:dict, r_primers:dict, database):
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
    extend_sbjct(forward_hsp_object, database, f_primers)
    extend_sbjct(reverse_hsp_object, database, r_primers)
    is_snp_primer_search(forward_hsp_object, f_primers, r_primers)
    is_snp_primer_search(reverse_hsp_object, f_primers, r_primers)
    is_distance(forward_hsp_object, reverse_hsp_object, amplicon_sequences)

    assert forward_hsp_object.valid == reverse_hsp_object.valid
    if forward_hsp_object.valid == True:
        if (forward_hsp_object.snp == False and reverse_hsp_object.snp == False and
                (forward_hsp_object.snp_match != reverse_hsp_object.snp_match or
                (forward_hsp_object.snp_match == None and reverse_hsp_object.snp_match == None))):
            forward_hsp_object.epcr = reverse_hsp_object.pcr_distance
            reverse_hsp_object.epcr = reverse_hsp_object.pcr_distance
        else:
            forward_hsp_object.epcr, reverse_hsp_object.epcr = False, False
    else:
        forward_hsp_object.epcr, reverse_hsp_object.epcr = False, False

def hamming_dist(seq1, seq2):
    """ Calculate the hamming distance between two sequences
    return: distance between the two sequences
    """
    #TODO: check out this case!!!
    # assert len(seq1) == len(seq2)

    dist = sum(x != y for x, y in zip(seq1, seq2))
    return(dist)

def f_primer_on_left_snp_check(f_prim_end, r_prim_r_comp, db_3_end_left, db_3_end_right):
    # Determines if a SNP is present in the case where the forward primer is oriented on the left

    return(hamming_dist(f_prim_end, db_3_end_left) <= MIN_SNP_HAM_DIST or
           hamming_dist(r_prim_r_comp, db_3_end_right) <= MIN_SNP_HAM_DIST)

def r_primer_on_left_snp_check(r_primer_end, f_prim_r_comp, db_3_end_left, db_3_end_right):
    #Determines if a SNP is present in the case where the reverse primer is oriented on the left

    return(hamming_dist(r_primer_end, db_3_end_left) <= MIN_SNP_HAM_DIST or
           hamming_dist(f_prim_r_comp, db_3_end_right) <= MIN_SNP_HAM_DIST)

def get_orientation(f_primer, r_primer, hsp_object):
    """ Determines the orientation of a sbjct on a contig, according to the extension of the 3'end of the sbjct.

    :param f_primer:
    :param r_primer:
    :param hsp_object:
    :return:
    """
    f_prim_left = Seq(str(f_primer[len(f_primer) - int(len(f_primer) / 2) : ]), generic_dna)
    f_prim_right = f_prim_left.reverse_complement()
    r_prim_left = Seq(str(r_primer[len(r_primer) - int(len(f_primer) / 2) : ]), generic_dna)
    r_prim_right = r_prim_left.reverse_complement()

    db_3_end_left = hsp_object.sbjct[len(hsp_object.sbjct) - int(len(f_primer) / 2) : ]
    db_3_end_right = hsp_object.sbjct[ : int(len(f_primer) / 2)]

    ham_f_prim_on_left = min(hamming_dist(f_prim_left, db_3_end_left),
                             hamming_dist(r_prim_right, db_3_end_right))
    ham_r_prim_on_left = min(hamming_dist(r_prim_left, db_3_end_left),
                             hamming_dist(f_prim_right, db_3_end_right))

    if ham_f_prim_on_left < ham_r_prim_on_left and \
                    abs(ham_f_prim_on_left - ham_r_prim_on_left) > 1:
        return True
    elif ham_f_prim_on_left > ham_r_prim_on_left and \
                    abs(ham_f_prim_on_left - ham_r_prim_on_left) > 1:
        return False
    else:
        return None

def is_snp_primer_search(hsp_object: HSP, f_primers:dict, r_primers:dict):
    """ Modifies hsp object to indicate any SNP within SNP_THRESHOLD of 3' end. Exceptions are made

    :param hsp_object: a HSP object
    :param f_primers: A dictionary containing f_primer seq'ns
    :param r_primers: A dictionary containing r_primer seq'ns
    :return:None
    """

    if hsp_object.name in f_primers and hsp_object.name in r_primers:
        f_primer = f_primers[hsp_object.name]
        r_primer = r_primers[hsp_object.name]
        f_prim_left = Seq(str(f_primer[len(f_primer) - SNP_THRESHOLD : ]), generic_dna)
        f_prim_right = f_prim_left.reverse_complement()
        r_prim_left = Seq(str(r_primer[len(r_primer) - SNP_THRESHOLD : ]), generic_dna)
        r_prim_right = r_prim_left.reverse_complement()

        db_3_end_left = hsp_object.sbjct[len(hsp_object.sbjct) - SNP_THRESHOLD : ]
        db_3_end_right = hsp_object.sbjct[ : SNP_THRESHOLD]

        f_primer_on_left_ham_check = f_primer_on_left_snp_check(f_prim_left, r_prim_right, db_3_end_left, db_3_end_right)
        r_primer_on_left_ham_check = r_primer_on_left_snp_check(r_prim_left, f_prim_right, db_3_end_left, db_3_end_right)
        orientation = get_orientation(f_primer, r_primer, hsp_object)

        hsp_object.snp_match = orientation
        if orientation == True:
            if f_primer_on_left_ham_check:
                hsp_object.snp = False
                if 'cj0307' in hsp_object.name and (db_3_end_right == "TTCC" or db_3_end_left == "GGAA"):
                    hsp_object.snp = True
                elif 'cj0177' in hsp_object.name and (db_3_end_right == "TTGA" or db_3_end_left == "TCAA"):
                    hsp_object.snp = True
            else:
                hsp_object.snp = True
        elif orientation == False:
            if r_primer_on_left_ham_check:
                hsp_object.snp = False
                if 'cj0307' in hsp_object.name and (db_3_end_right == "GGAA" or db_3_end_left == "TTCC"):
                    hsp_object.snp = True
                elif 'cj0177' in hsp_object.name and (db_3_end_right == "TCAA" or db_3_end_left == "TTGA"):
                    hsp_object.snp = True
            else:
                hsp_object.snp = True
        else:
            if (r_primer_on_left_ham_check or f_primer_on_left_ham_check):
                hsp_object.snp = False
            else:
                hsp_object.snp = True


def valid_dir(hsp: HSP):
    """ Modifies hsp object to determine if facing end of contig and
        within MAX_PERC_END of the amp len from the end of the contig.

    :param hsp: A HSP object
    :return: None
    """
    dist_end = abs((hsp.db_length + 1) - hsp.start - hsp.amp_len)
    dist_start = abs(hsp.start - hsp.amp_len)

    #TODO: could there be a case where the sequence is found over the entire amp!?
    if dist_end <= (MAX_PERC_END * hsp.amp_len):
        hsp.location = True
        hsp.end_dist = dist_end
    elif dist_start <= (MAX_PERC_END * hsp.amp_len):
        hsp.location = True
        hsp.end_dist = dist_start
    else:
        hsp.location = False

def ehyb(blast_object, cutoff_gene_len=CUTOFF_GENE_LENGTH):
    """ Determines if hsp objects in blast_object are of length >= CUTOFF_GENE_LENGTH

    :param blast_object: a Blastn object containing hsp objects
    :return: A list of HSP's with ehybrid attribute initialized
    """
    lo_ehybrid_hsp = [hsp for hsp in blast_object if hsp.length >= cutoff_gene_len]
    lo_failures = [hsp for hsp in blast_object if hsp.length < cutoff_gene_len]

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
    elif getattr(hsp, attr_i) == None:
        return i
    elif getattr(hsp, attr_i):
        return false_neg_pred(attr_bin_tree, hsp, (2 * i) + 1)
    elif not getattr(hsp, attr_i):
        return false_neg_pred(attr_bin_tree, hsp, (2 * i) + 2)

def max_bits(seqn_dir:str) -> dict:
    max_bits_dict = {}
    for record in SeqIO.parse(seqn_dir, "fasta"):
        if "11168_" not in record.id:
            max_bits_dict['11168_' + record.id] = (LAMBDA * len(record.seq) - log(K)) / log(2)
        else:
            max_bits_dict[record.id] = (LAMBDA * len(record.seq) - log(K)) / log(2)
    return max_bits_dict


def cj0181_missing_seq(hsp_object, primer_dict, database, chimeric_seq) -> bool:
    """ Determines if the missing cj0181 sequence is the chimeric primer sequence for its exception

    :param hsp_object:
    :param primer_dict:
    :param database:
    :param chimeric_seq:
    :return:
    """
    start = hsp_object.start
    end = hsp_object.end
    query_start = hsp_object.query_start
    query_end = hsp_object.query_end
    len_missing = len(primer_dict[hsp_object.name]) - abs(end - start) - 1
    begin_missing = query_start - 1
    end_missing = abs(len(primer_dict[hsp_object.name]) - query_end)
    if len_missing == 7:
        with open(database, 'r') as fasta:
            for contig in SeqIO.parse(fasta, 'fasta'):
                if contig.name == hsp_object.contig_name:
                    if start > end:
                        seq_found = contig[end - end_missing - 1: start + begin_missing]
                        hsp_object.sbjct = str(seq_found.reverse_complement().seq)
                    else:
                        seq_found = contig[start - begin_missing - 1: end + end_missing]
                        hsp_object.sbjct = str(seq_found.seq)
    missing_seq_found = hsp_object.sbjct[ : 7]
    r_comp_chimeric_seq = str(Seq(chimeric_seq).reverse_complement())
    return missing_seq_found == chimeric_seq or missing_seq_found == r_comp_chimeric_seq


def bsr(blast_object:Blastn, max_bits_dict:dict):
    """ Removes hsp object from blast_object if BSR is not >= MIN_BSR

    :param blast_object: A Blastn Object
    :param max_bits_dict: A dictionary containing the max_bs of query seq'ns
    :return: None
    """

    for hsp in blast_object.hsp_objects:
        hsp.bsr = hsp.bits / max_bits_dict[hsp.name]

        if hsp.bsr < MIN_BSR:
            blast_object.remove_hsp_object_all(hsp)

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
        if abs(f_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(r_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(f_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(r_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                and f_hsp.contig_name == hsp.contig_name:
            f_hsp.ehybrid, r_hsp.ehybrid = True, True
            f_hsp.amp_len, r_hsp.amp_len = hsp.length, hsp.length
            f_hsp.amp_sbjct, r_hsp.amp_sbjct = hsp.sbjct, hsp.sbjct
            f_hsp.amp_query, r_hsp.amp_query = hsp.query, hsp.query
    for hsp in ehybrid_qcov_fail:
        # if f_hsp.name in hsp.name and r_hsp.name in hsp.name:
        if abs(f_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(r_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(f_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                or abs(r_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
                and r_hsp.contig_name == hsp.contig_name:
            f_hsp.ehybrid, r_hsp.ehybrid = False, False
            f_hsp.amp_len, r_hsp.amp_len = hsp.length, hsp.length
            f_hsp.amp_sbjct, r_hsp.amp_sbjct = hsp.sbjct, hsp.sbjct
            f_hsp.amp_query, r_hsp.amp_query = hsp.query, hsp.query

def same_contig_pred(lo_tup_same_contig, full_blast_qcov, dict_f_primers, dict_r_primers,
                     max_f_bits_dict, max_r_bits_dict, amp_seq, database):
    """ Prediction of ecgf when primers are located on the same contig.

    :param lo_tup_same_contig: List of tuples of hsp primers found on the same contig
    :param full_blast_qcov: Blast object using amp as query and qcov included in restriction
    :param dict_f_primers: Dictionary containing forward primers
    :param dict_r_primers: dictionary containing reverse primers
    :param max_f_bits_dict: dict containing max bits for each of the forward primers
    :param max_r_bits_dict: dict containing max bits for each of the reverse primers
    :return: list of results with a dict of ecgf results.
    """
    result_dict = defaultdict(list)

    for tup in lo_tup_same_contig:
        f_hsp_old = tup[0]
        r_hsp_old = tup[1]

        copied_o = contig_copy(f_hsp_old, r_hsp_old, max_f_bits_dict, max_r_bits_dict, True)
        f_hsp, r_hsp = copied_o[0], copied_o[1]

        epcr(f_hsp, r_hsp, amp_seq, dict_f_primers, dict_r_primers, database) #assigns valid and snp attributes and pcr_distance and epcr

        ehyb_both_prim_found(full_blast_qcov, f_hsp, r_hsp)

        if f_hsp.epcr and (f_hsp.ehybrid or r_hsp.ehybrid):
            result_dict[f_hsp.name].append(f_hsp)
            result_dict[r_hsp.name].append(r_hsp)

    return result_dict

def diff_contig_pred(lo_tup_diff_contig, max_f_bits_dict, max_r_bits_dict, ehybrid_hsp_pass):
    """ Prediction of cgf when primers are located on different contigs.

    :param lo_tup_diff_contig: List of tuples containing f and r primers from the same query but on diff contigs
    :param max_f_bits_dict: dict containing max bits for each of the forward primers
    :param max_r_bits_dict: dict containing max bits for each of the reverse primers
    :param ehybrid_hsp_pass: A list of gene hsp that were found using ehybrid only
    :return: list of results with a dict of ecgf results.
    """

    result_dict = defaultdict(list)

    for tup in lo_tup_diff_contig:
        f_hsp_old = tup[0]
        r_hsp_old = tup[1]
        copied_o = contig_copy(f_hsp_old, r_hsp_old, max_f_bits_dict, max_r_bits_dict, False)
        f_hsp, r_hsp = copied_o[0], copied_o[1]

        #TODO: make helper functions
        for hsp in ehybrid_hsp_pass:
            if f_hsp.start == hsp.start or f_hsp.end == hsp.end \
                    and f_hsp.contig_name == hsp.contig_name:
                f_hsp.ehybrid = True
                f_hsp.amp_len = hsp.length
                valid_dir(f_hsp)
                f_hsp.amp_sbjct = hsp.sbjct
                f_hsp.amp_query = hsp.query
            if r_hsp.start == hsp.end or r_hsp.end == hsp.start \
                    and r_hsp.contig_name == hsp.contig_name:
                r_hsp.ehybrid = True
                r_hsp.amp_len = hsp.length
                valid_dir(r_hsp)
                r_hsp.amp_sbjct = hsp.sbjct
                r_hsp.amp_query = hsp.query
            if r_hsp.start == hsp.start or r_hsp.end == hsp.end \
                    and r_hsp.contig_name == hsp.contig_name :
                r_hsp.ehybrid = True
                r_hsp.amp_len = hsp.length
                valid_dir(r_hsp)
                r_hsp.amp_sbjct = hsp.sbjct
                r_hsp.amp_query = hsp.query
            if f_hsp.start == hsp.end or f_hsp.end == hsp.start \
                    and f_hsp.contig_name == hsp.contig_name :
                f_hsp.ehybrid = True
                f_hsp.amp_len = hsp.length
                valid_dir(f_hsp)
                f_hsp.amp_sbjct = hsp.sbjct
                f_hsp.amp_query = hsp.query

        if (f_hsp.ehybrid and f_hsp.location) and (r_hsp.ehybrid and r_hsp.location):
            result_dict[f_hsp.name].append(f_hsp)
            result_dict[r_hsp.name].append(r_hsp)

        # for hsp in ehybrid_hsp_fail:
        #     if abs(f_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
        #             or abs(f_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
        #             and f_hsp.contig_name == hsp.contig_name:
        #         f_hsp.ehybrid = False
        #         f_hsp.amp_len = hsp.length
        #         f_hsp.amp_query = hsp.query
        #         f_hsp.amp_sbjct = hsp.sbjct
        #     if abs(r_hsp.start - hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
        #             or abs(r_hsp.end - hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * hsp.length) \
        #             and r_hsp.contig_name == hsp.contig_name:
        #         r_hsp.ehybrid = False
        #         r_hsp.amp_len = hsp.length
        #         r_hsp.amp_sbjct = hsp.sbjct
        #         r_hsp.amp_query = hsp.query

    return result_dict

def single_primer_found(lo_hsp_single_primers, ehybrid_hsp_pass):
    """ Prediction of cgf when only one primer is found.

    :param lo_hsp_single_primers: list of hsp's where only one primer is found.
    :param ehybrid_hsp_pass: A list of gene hsp that were found using ehybrid only
    :param ehybrid_hsp_fail: A list of gene hsp that were not found using ehybrid only b/c not long enough length
    :return: list of results with a dict of ecgf results .
    """
    result_dict = defaultdict(list)

    for single_hsp in lo_hsp_single_primers:
        single_hsp.both_primers_found = False
        single_hsp.contig = False

        for blast_hsp in ehybrid_hsp_pass:
            if abs(single_hsp.start - blast_hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) or \
                                    abs(single_hsp.end - blast_hsp.end) <= \
                                    (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) and \
                                    single_hsp.contig_name == blast_hsp.contig_name:
                single_hsp.ehybrid = True
                single_hsp.amp_len = blast_hsp.length
                single_hsp.amp_sbjct = blast_hsp.sbjct
                single_hsp.amp_query = blast_hsp.query
                valid_dir(single_hsp)
                if single_hsp.location == True:
                    result_dict[single_hsp.name].append(single_hsp)
                    print(single_hsp.name, 'added b/c of single primer found')

        # for blast_hsp in ehybrid_hsp_fail:
        #     if abs(single_hsp.start - blast_hsp.start) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) or \
        #                             abs(single_hsp.end - blast_hsp.end) <= (MAX_PERC_EHYB_PRIMER_ENDS * blast_hsp.length) \
        #                     and single_hsp.contig_name == blast_hsp.contig_name:
        #         single_hsp.ehybrid = False
        #         single_hsp.amp_len = blast_hsp.length
        #         single_hsp.amp_sbjct = blast_hsp.sbjct
        #         single_hsp.amp_query = blast_hsp.query

    return result_dict

#TODO: delete file_dict from input (only needed for fourth case check)
#TODO: make exceptions into helper functions
# @profile
def ecgf(forward_primers:str, reverse_primers:str, database:str, amp_sequences:str, cj0181_f_primer, cj0181_r_primer) -> list:
    """ Predicts in vitro cgf (eCGF)

    :param forward_primers: A string containing the path to all the 40 primer sequences in a fasta file
    :param reverse_primers: A string containing the path to all the 40 reverse primer sequences in a fasta file
    :param database: A string containing the path to a directory that contains strains of interest in fasta files
    :param amp_sequences: A string containing the path to all the 40 amplicon sequences in a fasta file
    :return: A list containing a dictionary of gene found (key = hsp.name, val = hsp)
    """


    #bsr
    max_f_bits_dict = max_bits(forward_primers)
    max_r_bits_dict = max_bits(reverse_primers)

    #blastn
    forward_blast_bsr = create_blastn_bsr_object(forward_primers, database)
    reverse_blast_bsr = create_blastn_bsr_object(reverse_primers, database)
    blast_object_hsps = create_blastn_object(amp_sequences, database)
    full_blast_qcov_hsps = create_blastn_object(amp_sequences, database, True)
    #for 4th case exception check
    ehyb_blast_hsps = create_blastn_object(amp_sequences, database, False, 70)

    #variable setup
    dict_f_primers = create_dict_from_fasta_seqs(forward_primers)
    dict_r_primers = create_dict_from_fasta_seqs(reverse_primers)
    # dict_amp = create_dict_from_fasta_seqs(amp_sequences)
    bsr(forward_blast_bsr, max_f_bits_dict)
    bsr(reverse_blast_bsr, max_r_bits_dict)

    #predictions
    results_list = four_branch_prediction(forward_blast_bsr, reverse_blast_bsr, full_blast_qcov_hsps,
                                              dict_f_primers, dict_r_primers, max_f_bits_dict, max_r_bits_dict,
                                              amp_sequences, database, blast_object_hsps)


    result_dict = results_list[0]
    single_primer_results_dict = results_list[1]


    #CASE 4: exceptions
    lo_hsp_ehybrid = ehyb(ehyb_blast_hsps)  # assigns ehybrid attributes to each hsp from amp vs db
    ehybrid_pass = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == True]
    ehyb_pos = [hsp for hsp in ehybrid_pass if hsp.name not in result_dict.keys()
                if hsp.name not in single_primer_results_dict.keys()]
    # TODO: below is so I could see if I could find F+ causes for Steven's analysis
    # ehyb_pos = [hsp for hsp in ehybrid_pass]
    ehyb_names = [hsp.name for hsp in ehyb_pos]

    #cj0181 case exception - AGGATTA 7bp 5'end
    exception_result_dict = {}
    if '11168_cj0181' in ehyb_names:
        f_blast_exc = create_blastn_object_exceptions(cj0181_f_primer, database, 60)
        r_blast_exc = create_blastn_object_exceptions(cj0181_r_primer, database, 60)
        lo_f_primers = []
        lo_r_primers = []

        for hsp in f_blast_exc.hsp_objects:
            if cj0181_missing_seq(hsp, dict_f_primers, database, "AGGATTA"):
                lo_f_primers.append(hsp)
        for hsp in r_blast_exc.hsp_objects:
            lo_r_primers.append(hsp)

        f_blast_exc.hsp_objects = lo_f_primers
        r_blast_exc.hsp_objects = lo_r_primers

        if len(lo_f_primers) > 0 and len(lo_r_primers) > 0 :
            exception_result_dict = four_branch_prediction(f_blast_exc, r_blast_exc, ehyb_pos, dict_f_primers, dict_r_primers, max_f_bits_dict,
                                       max_r_bits_dict, amp_sequences, database, ehyb_pos)[0]

    #TODO: delete print statements
    # for hsp in ehyb_pos:
    #     print(hsp.name, hsp.sbjct)
    print('genome: ', database)
    print('found using ecgf: ', result_dict.keys())
    print('found using ehyb: ', ehyb_names)
    print('found using exception: ', exception_result_dict.keys())
    print(os.path.basename(database))
    print('len genes found including exceptions', len(result_dict) + len(exception_result_dict))

    # if os.path.basename(database) in file_dict.keys():
        # ehyb_pos_with_exceptions = [hsp for hsp in ehyb_pos if hsp.name not in exception_result_dict.keys()]
        # fourth_case_check(ehyb_pos_with_exceptions, result_dict, dict_f_primers, dict_r_primers, dict_amp,
        #                   os.path.basename(database), file_dict)


    return [results_list[0], ehyb_pos, exception_result_dict, single_primer_results_dict]

def assign_bsr(f_hsp_single_primers, r_hsp_single_primers, max_f_bits_dict, max_r_bits_dict):

    lo_hsp_single_primers = []
    for f_hsp in f_hsp_single_primers:
        f_hsp.bsr = f_hsp.bits / max_f_bits_dict[f_hsp.name]
        if f_hsp.bsr >= MIN_BSR:
            lo_hsp_single_primers.append(f_hsp)
    for r_hsp in r_hsp_single_primers:
        r_hsp.bsr = r_hsp.bits / max_r_bits_dict[r_hsp.name]
        if r_hsp.bsr >= MIN_BSR:
            lo_hsp_single_primers.append(r_hsp)
    return lo_hsp_single_primers

#TODO: split into helper functions!!!
def four_branch_prediction(forward_blast, reverse_blast, full_blast_qcov_hsps, dict_f_primers, dict_r_primers, max_f_bits_dict,
                           max_r_bits_dict, amp_sequences, database, blast_object_hsps):

    lo_tup_same_queries = [(f_hsp,r_hsp) for f_hsp in forward_blast.hsp_objects
                           for r_hsp in reverse_blast.hsp_objects
                           if f_hsp.name == r_hsp.name]
    lo_tup_same_contig = [tup for tup in lo_tup_same_queries if tup[0].contig_name == tup[1].contig_name]

    results_dict_same_contig = same_contig_pred(lo_tup_same_contig, full_blast_qcov_hsps, dict_f_primers,
                                                dict_r_primers, max_f_bits_dict, max_r_bits_dict, amp_sequences, database)

    #doesn't look for genes on different contigs if they were already found on same contig
    lo_tup_diff_contig = [tup for tup in lo_tup_same_queries
                          if tup not in lo_tup_same_contig and tup[0].name not in results_dict_same_contig]

    lo_hsp_ehybrid = ehyb(blast_object_hsps, 150)
    ehybrid_hsp_pass = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == True]
    # ehybrid_hsp_fail = [hsp for hsp in lo_hsp_ehybrid if hsp.ehybrid == False]

    #Different contig prediction
    # diff_contig = diff_contig_pred(lo_tup_diff_contig, max_r_bits_dict, max_r_bits_dict, ehybrid_hsp_pass)
    # results_diff_contig = diff_contig[0]
    results_diff_contig = diff_contig_pred(lo_tup_diff_contig, max_r_bits_dict, max_r_bits_dict, ehybrid_hsp_pass)
    diff_contig_results_keys = [key for key in results_diff_contig.keys()]
    same_contig_result_keys = [key for key in results_dict_same_contig.keys()]

    f_hsp_single_primers = [hsp for hsp in forward_blast.hsp_objects
                            if hsp.name not in same_contig_result_keys and hsp.name not in diff_contig_results_keys]
    r_hsp_single_primers = [hsp for hsp in reverse_blast.hsp_objects
                            if hsp.name not in same_contig_result_keys and hsp.name not in diff_contig_results_keys]

    ehybrid_hsp_pass_single = [hsp for hsp in ehybrid_hsp_pass if (hsp.identities / hsp.length) >= SINGLE_PRIMER_ID]
    # ehybrid_hsp_fail_single = [hsp for hsp in ehybrid_hsp_fail if (hsp.identities / hsp.length) >= SINGLE_PRIMER_ID]

    lo_hsp_single_primers = assign_bsr(f_hsp_single_primers, r_hsp_single_primers, max_f_bits_dict, max_r_bits_dict)
    results_one_primer = single_primer_found(lo_hsp_single_primers, ehybrid_hsp_pass_single)

    #combine all results
    result_dict = defaultdict(list)
    single_primer_results_dict = defaultdict(list)
    for k,v in results_dict_same_contig.items():
        result_dict[k].append(v)
    for k,v in results_diff_contig.items():
        result_dict[k].append(v)
    for k,v in results_one_primer.items():
        single_primer_results_dict[k].append(v)

    results_list = [result_dict, single_primer_results_dict]

    gc.collect()
    return results_list

# #TODO: remove from program or make seperate case for debugging in the future!
# def fourth_case_check(ehyb_pos, result_dict, f_primer_dict, r_primer_dict, amp_dict, genome_name, lab_result_file_dict):
#     file = open("/home/sfisher/eCGF/all_genomes/Results/delete", "a")
#     ehyb_pos_names = [hsp.name for hsp in ehyb_pos]
#     file.write('\n \n Genome name: ' + str(genome_name))
#     file.write('\n Genes not found using eCGF but found using eHYB ' + str(ehyb_pos_names))
#     file.write('\n Number of genes found in eHYB only ' + str(len(ehyb_pos)))
#     file.write('\n Number of genes found using eCGF (Probability of being TRULY +ve is ~high) ' +
#                str(len(result_dict.keys())))
#     file.write('\n Number of genes not found at all (Probability of being TRULY -ve is high) ' +
#                str(abs(40 - len(result_dict.keys()) - len(ehyb_pos))))
#     print(lab_result_file_dict)
#     if '.fasta' not in genome_name:
#         lab_results = lab_result_file_dict[genome_name + '.fasta']
#     else:
#         lab_results = lab_result_file_dict[genome_name]
#     file.write('\n Number of genes found in lab ' + str(len(lab_results)))
#     file.write('\n Below is more info on the genes found using eHYB only (was -ve in eCGF but +ve in eHYB): ')
#     print(lab_results)
#     for hsp in ehyb_pos:
#         hsp_name = hsp.name
#         perc_id = hsp.identities / hsp.length
#         lab_result = "Found" if hsp.name in lab_results else "Not Found"
#         file.write('\n \n' + hsp_name)
#         file.write('\n Lab Result: ' + lab_result)
#         file.write('\n contig: ' + hsp.contig_name)
#         file.write('\n % id ' + str(perc_id))
#         file.write('\n BSR: ' + str(hsp.bsr))
#         file.write('\n End Dist: ' + str(hsp.end_dist))
#         file.write(
#             '\n qcov (% of query seqn that overlaps sbjct seqn) ' + str(len(hsp.query) / (len(amp_dict[hsp.name]) + 1)))
#         file.write('\n Entire Query Seq      ' + str(len(amp_dict[hsp.name]) + 1) + " 1 " + str(
#             len(amp_dict[hsp.name])) + " " + str(amp_dict[hsp.name]))
#         file.write('\n Matching Query Seq    ' + str(len(hsp.query)) + " " + str(hsp.query_start) + " " + str(
#             hsp.query_end) + " " * (hsp.query_start + (
#         7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.query)
#         file.write('\n Matching Sbjct Seq    ' + str(len(hsp.sbjct)) + " " + str(hsp.query_start) + " " + str(
#             hsp.query_end) + " " * (hsp.query_start + (
#         7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.sbjct)
#         file.write(
#             '\n Forward primer         ' + str(len(f_primer_dict[hsp_name])) + "       " + str(f_primer_dict[hsp_name]))
#         file.write('\n Reverse primer         ' + str(len(r_primer_dict[hsp_name])) + " " * (
#             len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]) + 7) + str(
#                 amp_dict[hsp.name][len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]):]))
#
#         # if (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
#         file.write('\n Forward Primer qcov: ' + str(
#             (len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name])) + " (CUTOFF: " + str(
#             (QCOV_HSP_PERC / 100)) + " )")
#         # if (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
#         file.write('\n Reverse Primer qcov: ' + str(
#             (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(
#                 r_primer_dict[hsp_name])) + " (CUTOFF:  " + str((QCOV_HSP_PERC / 100)) + " )")
#         # if perc_id * 100 < PERC_ID_CUTOFF:
#         file.write('\n Perc ID for ehyb in eCGF: ' + str(perc_id) + " (CUTOFF: " + str((PERC_ID_CUTOFF / 100)) + " )")
#     file.close()
#
# #TODO: remove from program or make seperate for debugging in the future!!
# def per_gene_fourth_case_check(all_ehyb_pos, forward_primers, reverse_primers, amplicon_sequences, lab_results_dict):
#
#     f_primer_dict = create_dict_from_fasta_seqs(forward_primers)
#     r_primer_dict = create_dict_from_fasta_seqs(reverse_primers)
#     amp_dict = create_dict_from_fasta_seqs(amplicon_sequences)
#     gene_ehyb_pos_dict = defaultdict(list)
#
#     for genome, ehyb_pos in all_ehyb_pos.items():
#         for hsp in ehyb_pos:
#             print('genome', genome)
#             gene_ehyb_pos_dict[hsp.name].append((hsp,genome))
#
#     file = open("/home/sfisher/eCGF/all_genomes/Results/delete", "a")
#     for gene_name, lo_hsp_tup in gene_ehyb_pos_dict.items():
#         file.write("\n \n  " + gene_name + " not found using eCGF but found using eHYB " + str(len(lo_hsp_tup)) + " times.")
#         for hsp_tup in lo_hsp_tup:
#             hsp, genome = hsp_tup
#             try:
#                 found = lab_results_dict[genome + '.fasta']
#             except:
#                 found = lab_results_dict[genome]
#             hsp_name = hsp.name
#             perc_id = hsp.identities / hsp.length
#             lab_result = "Found" if hsp.name in found else "Not Found"
#             file.write('\n \n' + hsp_name)
#             file.write('\n Lab Result: '  + lab_result)
#             file.write('\n Contig: ' + hsp.contig_name)
#             file.write('\n % id ' + str(perc_id))
#             file.write('\n qcov (% of query seqn that overlaps sbjct seqn) ' + str(len(hsp.query) / (len(amp_dict[hsp.name]) + 1)))
#             file.write('\n Entire Query Seq      ' + str(len(amp_dict[hsp.name]) + 1) + " 1 " + str(len(amp_dict[hsp.name])) + " " + str(amp_dict[hsp.name]))
#             file.write('\n Matching Query Seq    ' + str(len(hsp.query)) + " " + str(hsp.query_start) + " " + str(hsp.query_end) + " "*(hsp.query_start + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.query)
#             file.write('\n Matching Sbjct Seq    ' + str(len(hsp.sbjct)) + " " + str(hsp.query_start) + " " + str(hsp.query_end) + " "*(hsp.query_start + (7 - len(str(len(hsp.query))) - len(str(hsp.query_start)) - len(str(hsp.query_end)))) + hsp.sbjct)
#             file.write('\n Forward primer         ' + str(len(f_primer_dict[hsp_name])) + "       " + str(f_primer_dict[hsp_name]))
#             file.write('\n Reverse primer         ' + str(len(r_primer_dict[hsp_name])) + " "*(len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]) + 7) + str(amp_dict[hsp.name][len(amp_dict[hsp.name]) - len(r_primer_dict[hsp_name]): ] ))
#             file.write('\n Forward Primer qcov: ' + str((len(f_primer_dict[hsp_name]) - hsp.query_start + 1) / len(f_primer_dict[hsp_name])) + " (CUTOFF: " + str((QCOV_HSP_PERC / 100)) + " )")
#             # if (len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name]) < (QCOV_HSP_PERC / 100):
#             file.write('\n Reverse Primer qcov: ' + str((len(r_primer_dict[hsp_name]) - (len(amp_dict[hsp.name]) - hsp.query_end)) / len(r_primer_dict[hsp_name])) + " (CUTOFF:  " + str((QCOV_HSP_PERC / 100)) + " )")
#             # if perc_id * 100 < PERC_ID_CUTOFF:
#             file.write('\n Perc ID for ehyb in eCGF: ' + str(perc_id) + " (CUTOFF: " + str((PERC_ID_CUTOFF / 100)) + " )")
#
#     file.close()

def find_closest_fingerprint(bin_results:list, cgftypes_file:str) -> list:
    """ Return the closest fingerprint match of the binary results based off of the fingerprints in cgftypes_file

    """

    with open(cgftypes_file, "r") as fingerprint_track:
        csvReader = csv.reader(fingerprint_track)
        header = next(csvReader)
        cgf_type_index = header.index("cgf.type")
        strain_freq_index = header.index("num.strains")

        dict_gene_indicies = {gene : header.index(gene) for gene in GENE_LIST}
        found_bin_dict = {gene: bin_results[index] for index, gene in enumerate(GENE_LIST)}

        largest_lo_tup = []
        largest_strain_freq = 0
        for row in csvReader:
            row_bin_dict = {gene: int(row[index]) for gene, index in dict_gene_indicies.items()}
            dict_differences = {val : found_bin_dict[val] - row_bin_dict[val] for val in found_bin_dict
                                if val in row_bin_dict and found_bin_dict[val] - row_bin_dict[val] != 0}

            if largest_lo_tup == []:
                largest_lo_tup = [(row[cgf_type_index], dict_differences)]
                largest_strain_freq = row[strain_freq_index]

            elif len(dict_differences) < len(largest_lo_tup[0][1]):
                largest_lo_tup = [(row[cgf_type_index], dict_differences)]

            elif len(dict_differences) == len(largest_lo_tup[0][1]):
                strain_freq = row[strain_freq_index]
                if strain_freq > largest_strain_freq:
                    largest_strain_freq = strain_freq
                    largest_lo_tup = [(row[cgf_type_index], dict_differences)] + largest_lo_tup
                else:
                    largest_lo_tup.append((row[cgf_type_index], dict_differences))

    print(largest_lo_tup)
    return largest_lo_tup

def make_binary(lo_result):
    if 2 not in lo_result:
        return list(lo_result)
    else:
        for result in lo_result:
            if result == 2:
                count = lo_result.index(result)
                tmp = list(lo_result)
                tmp[count] = 0
                lo_result_0 = make_binary(tmp)
                tmp[count] = 1
                lo_result_1 = make_binary(tmp)
                return [lo_result_0, lo_result_1]

def flatten(items, seqtypes=(list, tuple)):
    for i, x in enumerate(items):
        while i < len(items) and isinstance(items[i], seqtypes):
            items[i:i+1] = items[i]
    return items

def chunk(iterable, length):
    assert len(iterable)  % length == 0
    unit = []
    for i in iterable:
        unit.append(i)
        if len(unit) == length:
            yield unit
            unit = []

def closest_matches(genome, lo_result):
    if 2 in lo_result:
        make_binary.lo_solutions = []
        chained_lo_bin_results = list(chunk(flatten(make_binary(lo_result)), 40))
    else:
        chained_lo_bin_results = [lo_result]

    closest_match = []
    same_sim_2_case = []
    lo_result_to_display = []
    #TODO: include fprints in fastas or ask user to input them?... maybe have option to include??
    for lo_bin_result in chained_lo_bin_results:
        match = find_closest_fingerprint(lo_bin_result,
                                         "/home/sfisher/eCGF/all_genomes/cgf_types_fprints_26149.csv")
        if closest_match == []:
            closest_match = match
            lo_result_to_display = lo_bin_result
        elif len(closest_match[0][1]) == len(match[0][1]):
            for x in match:
                same_sim_2_case.append(x)
        elif len(closest_match[0][1]) > len(match[0][1]):
            closest_match = match
            lo_result_to_display = lo_bin_result

    return [closest_match, lo_result_to_display, same_sim_2_case]

def na_empty_lists(list):
    if len(list) < 1:
        return ['NA']
    else:
        return list

def write_result_to_file(cgf_predictions_dict, exception_dict_result, single_primer_results_dict,
                         error_rates_file, out_results):

    error_rate_dict = {}
    with open(error_rates_file, 'r') as error_rates :
        csvReader = csv.reader(error_rates)
        header = next(csvReader)
        gene_index = header.index("gene")
        error_rate_index = header.index("error_rate")

        for row in csvReader:
            gene = row[gene_index]
            error_rate_val = row[error_rate_index]
            error_rate_dict['11168_' + gene] = error_rate_val

    results = {genome: [1 if gene in cgf_predictions_dict[genome]
                            else 2 if gene in exception_dict_result[genome]
                            else 2 if gene in single_primer_results_dict[genome]
                            else 0 for gene in GENE_LIST] for genome in cgf_predictions_dict.keys()}
    print('results', results)

    with open(out_results, 'a') as file:
        csv_writer = csv.writer(file, dialect='excel')

        csv_writer.writerow(['Genes:'] + GENE_LIST +
                            ['Nearest Match',
                             'Number of Similarities (/40)',
                             'Genes with Disagreement / Error Rate (avg: 0.00425)',
                             'Genes with 2:',
                             'Other cgf.type with same # of matches from 2 Case: ',
                             'Other cgf.type with same # of matches: '])


        for genome, lo_result in results.items():
            highest_sim_matches = closest_matches(genome, lo_result)
            closest_match = highest_sim_matches[0]
            fingerprint = highest_sim_matches[1]
            other_closest_matches = na_empty_lists([tup[0] for tup in closest_match[2:]])
            same_sim_2_case = na_empty_lists([tup[0] for tup in closest_matches(genome, lo_result)[2]])

            cgftype_match = closest_match[0][0]
            dict_mismatches = closest_match[0][1]
            num_sim = str(40 - len(dict_mismatches))
            print('dict mm', dict_mismatches.keys())
            print('error rate', error_rate_dict)
            contrasting_genes = na_empty_lists([(key, error_rate_dict[key]) for key in dict_mismatches.keys()])
            genes_with_2 = na_empty_lists([GENE_LIST[ind] for ind, x in enumerate(lo_result) if x == 2])

            row = [genome] + \
                  fingerprint + \
                  [cgftype_match, num_sim] + \
                  [contrasting_genes] + \
                  [genes_with_2] + \
                  [same_sim_2_case] + \
                  [(other_closest_matches)]
            csv_writer.writerow(row)


#TODO: changed db_dir to file_path when calling on each genome.
#TODO: put CASE 4: exceptions into another function and call that separately
def main(db_fasta, f_primers_fasta, r_primers_fasta, amp_fasta, cj0181_f_primer, cj0181_r_primer,
         out_results, error_rates_file):
    """

    :param db_directory: The location of the directory with fasta database files contained (already formatted using makeblastdb)
    :param f_primers_fasta: The fasta file location of the forward primers
    :param r_primers_fasta: The fasta file location of the reverse primers
    :param amplicon_sequences: The fasta file location of the amplicon sequences
    :return: A dictionary
    """
    #TODO: delete up to file.close() (for validation...fourth case check!)

    # lab_binary_results = "/home/sfisher/eCGF/all_genomes/Results/nov_1_labcgf_1046genomes.txt"
    # with open(lab_binary_results, "r") as file:
    #     next(file)
    #     lines = file.readlines()
    #     file_gene_dict = {}
    #     for line in lines:
    #         genes_expected = []
    #         count = 0
    #         for word in line.split(","):
    #             if word == '1' or word == '1\n':
    #                 genes_expected.append(GENE_LIST[count])
    #                 count += 1
    #             elif word == '0' or word == '0\n':
    #                 count += 1
    #         genome_name = line.split(',', 1)[0]
    #         if '.fasta' not in genome_name:
    #             file_gene_dict[genome_name + '.fasta'] = genes_expected
    #         else:
    #             file_gene_dict[genome_name] = genes_expected

    files = (file for file in os.listdir(db_fasta) if file.endswith('.fasta'))
    files_paths = []
    for file in files:
        files_paths.append(os.path.abspath(db_fasta) + '/' + file)
    cgf_predictions_dict = {}
    exception_dict_result = {}
    single_primer_results_dict = {}

    #TODO: added dict for testing
    # all_ehyb_pos_dict = {}

    for file_path in files_paths:
        file_name = Path(file_path).stem
        print('file_name', file_name)
        result = ecgf(f_primers_fasta, r_primers_fasta, file_path, amp_fasta, cj0181_f_primer, cj0181_r_primer)
        cgf_predictions_dict[file_name] = result[0]
        exception_dict_result[file_name] = result[2]
        single_primer_results_dict[file_name] = result[3]

        # #TODO: delete below... added for testing:
        # all_ehyb_pos_dict[file_name] = result[1]

        gc.collect()

    # per_gene_fourth_case_check(all_ehyb_pos_dict, f_primers_fasta, r_primers_fasta, amp_fasta, file_gene_dict)

    write_result_to_file(cgf_predictions_dict, exception_dict_result, single_primer_results_dict,
                         error_rates_file, out_results)

    print(cgf_predictions_dict)
    return cgf_predictions_dict


# TODO: Commented out so I could use as package in github. (running main from __main__.py)
if __name__ == "__main__":

    small_testset = "/home/sfisher/Sequences/11168_test_files/246_gnomes_2nd_tests"
    # debug_cases = "/home/sfisher/Sequences/11168_test_files/debug_genes/other"
    # stevens_genomes = "/home/sfisher/eCGF/genomes_for_Steven"
    # lab_binary_results = "/home/sfisher/Sequences/11168_test_files/cgf40_v2.txt"
    puppy_genomes = "/home/sfisher/eCGF/puppy_genomes/assemblies"
    all_genomes = "/home/sfisher/eCGF/all_genomes/genomes_with_lab_data"
    debug_genomes = "/home/sfisher/eCGF/all_genomes/debug_smaller"
    # out_results = "/home/sfisher/eCGF/all_genomes/testing_results/new_single_cell_Test.csv"
    out_results = "/home/sfisher/eCGF/all_genomes/Results/final_results_new.csv"
    # out_results = "/home/sfisher/eCGF/puppy_genomes/eCGF_results.csv"

    #TODO: add lab binary results for reference??
    lab_binary_results = "/home/sfisher/eCGF/all_genomes/Results/nov_1_labcgf_1046genomes.txt"
    error_rates_file = "/home/sfisher/eCGF/all_genomes/error_rates.txt"
    forward_primers = "/home/sfisher/eCGF/primer_amp_fastas/cgf_forward_primers.fasta"
    reverse_primers = "/home/sfisher/eCGF/primer_amp_fastas/cgf_reverse_primers.fasta"
    amplicon_sequences = "/home/sfisher/eCGF/primer_amp_fastas/amplicon_sequences.fasta"
    f_primer_file = "/home/sfisher/Sequences/BSR/f_primers/cj0181.fasta"
    r_primer_file = "/home/sfisher/Sequences/BSR/r_primers/cj0181.fasta"

    main(all_genomes, forward_primers, reverse_primers, amplicon_sequences, f_primer_file, r_primer_file, out_results, error_rates_file)


    # cProfile.run('cgf_prediction_trial(forward_primers, reverse_primers, test_cprofile_file, amplicon_sequences, max_f_bits_dict, max_r_bits_dict, max_amp_bits_dict)')
    # cProfile.run('main(all_genomes, forward_primers, reverse_primers, amplicon_sequences); print')
    # main(test_11168_cases, forward_primers, reverse_primers, amplicon_sequences)