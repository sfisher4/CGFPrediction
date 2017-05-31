from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

HSP_THRESHOLD = 1.0 #99% identity

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
# amplicons contains sequences to compare the length of the hsp with.
def create_hsp_records(blast_records, amplicons):
    dict_hsp = {}

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            #very naive... don't continue to parse once you found the primer
            for amplicon in SeqIO.parse(amplicons, "fasta"):
                if alignment.title.find(amplicon.name) != -1:
                    lo_hsp = [hsp for hsp in alignment.hsps if HSP_THRESHOLD == (hsp.identities / len(amplicon.seq))]
                    if (len(lo_hsp) != 0): #lo_hsp is not empty
                        dict_hsp[alignment] = lo_hsp
                    break
    return dict_hsp

#returns a list of keys in the dictionary input
def dict_to_list(dict):
    list = []
    for key in dict.keys():
        list.append(key)
    return list

def search_gene_directly(query_genes, database, out_file):
    blastn_query(query_genes, database, out_file)
    blast_records = create_blast_record(out_file)
    dict_hsp = create_hsp_records(blast_records, database)
    return dict_to_list(dict_hsp)

def main():
    # assume database is complete or primers are on the same contigs
    database = "/home/sfisher/Documents/example_genomes/complete/IA3902.fasta"  # contains id and a complete genome sequence
    entire_seq_out_file = "/home/sfisher/Sequences/blast_record/entire_seq_blast.xml"
    amplicon_sequences = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
    amplicon_sequence = "/home/sfisher/Sequences/amplicon_sequences/cj0008.fasta" #TODO: inputting contigs!!!


    # using amplicon_sequences as database for testing purposes
    print(len(search_gene_directly(amplicon_sequences, amplicon_sequences, entire_seq_out_file))) #TODO: returning 40? Why?
    print(search_gene_directly(amplicon_sequences, amplicon_sequence, entire_seq_out_file))

if __name__ == "__main__": main()
