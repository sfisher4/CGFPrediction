import string

#Creates two fasta files, one containing the forward sequence and one containing the reverse sequence from a text file (specific format given in text file)
def create_primer_files(forward_output, reverse_output, text_file):
    with open(text_file) as textfile:

        forward_primer_file = open(forward_output, "w")
        reverse_primer_file = open(reverse_output, "w")

        word_list = []
        for word in textfile.read().split():
            word_list.append(str(word))
        #print(word_list)
        i=0
        for word in word_list:
            #gene
            #if  i % 6 == 0:
                #print(word)
            #amplicon size (bp)
            #if (i - 1) % 6 == 0:
                #print(word)
            #forward primer
            if (i - 2) % 6 == 0:
                forward_primer_file.write(word_list[i - 2] + "\n"+ word + "\n")
            #reverse primer
            if (i - 3) % 6 == 0:
                reverse_primer_file.write(word_list[i - 3] + "\n" + word + "\n")
            #strand
            #if (i - 4) % 6 ==0:
                #print(word)
                """
            #amplicon 
            if (i - 5) % 6 == 0:
                file.write(word_list[i-5] + " " + word_list[i-4] + " " + word_list[i-1] + "\n" + word + "\n")
                """
            i += 1
        forward_primer_file.close()
        reverse_primer_file.close()

#Creates a fasta file from a text file (specific format given in text file) containing all elements of the text file
def create_primers_amplicon_fasta(output, text_file):
    with open(text_file) as textfile:

        file = open(output, "w")

        word_list = []
        for word in textfile.read().split():
            word_list.append(str(word))
        # print(word_list)
        i = 0
        for word in word_list:
            # gene
                # if  i % 6 == 0:
                    # print(word)
            # amplicon size (bp)
                # if (i - 1) % 6 == 0:
                    # print(word)
            # forward primer
            if (i - 2) % 6 == 0:
                file.write(word_list[i - 2] + " forward strand" "\n" + word + "\n")
            # reverse primer
            if (i - 3) % 6 == 0:
                file.write(word_list[i - 3] + " reverse strand" "\n" + word + "\n")
                # strand
                # if (i - 4) % 6 ==0:
                # print(word)

            #amplicon
            if (i - 5) % 6 == 0:
                file.write(word_list[i-5] + " " + word_list[i-4] + " " + word_list[i-1] + "\n" + word + "\n")

            i += 1
        file.close()

def text_file_split():
    file = open("/home/sfisher/Sequences/11168_test_files/cgf40_results_modified.txt")
    lines = file.readlines()

    file_list = []
    for line in lines:
        new_line = line.split()
        file_list.append(new_line)
    file.close()
    print(file_list)

    dict_result = {}
    for result_list in file_list:
        dict_result[result_list[0]] = result_list[1:]





forward_output = "/home/sfisher/Sequences/cgf_forward_primers.fasta"
reverse_output = "/home/sfisher/Sequences/cgf_reverse_primers.fasta"
text_file = "/home/sfisher/Sequences/cgf_primers_amplicons_modified.txt"
create_primer_files(forward_output, reverse_output, text_file)

complete_output = "/home/sfisher/Sequences/cgf_primers_amplicons_write.fasta"
create_primers_amplicon_fasta(complete_output, text_file)
text_file_split()