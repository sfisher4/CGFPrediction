

# TODO: change so you do not create a new file each time!!!


def create_multiple_fasta(query_file, type:str):
    file = open(query_file)
    line_num = 0
    lines = file.readlines()
    for line in lines:
        if ">" in line:
            line_name = line[1:-1]
            print(line_name)
            file = open("/home/sfisher/Sequences/BSR/" + type + '/'+ line_name + ".fasta", 'w')
            file.write(line)
            next_line = lines[line_num + 1]
            file.write(next_line)
            line_num += 2
            file.close()
    file.close()

f_query_file = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
r_query_file = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
amp_query_file = "/home/sfisher/Sequences/amplicon_sequences/amplicon_sequences.fasta"
create_multiple_fasta(f_query_file, 'f_primers')
create_multiple_fasta(r_query_file, 'r_primers')
create_multiple_fasta(amp_query_file, 'amp_seq')