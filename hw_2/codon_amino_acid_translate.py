import sys
def translate(table_file, amino_output, codon_input):
    translate_dict = {} #Key: codon, Value: amino acid, for translating codons to amino acids
    with open(table_file, "r") as tf:
        tf.readline(); #consume firstline
        cur_line = tf.readline()
        while(cur_line != ""):
            line_arr = cur_line.strip().split('\t')
            codon = line_arr[0]
            amino = line_arr[2]
            translate_dict[codon] = amino
            cur_line = tf.readline()
    #translation dictionary complete

    #begin creating a new sequence of amino acids, translated from dna sequence using dictionary created above.
    amino_sqn = ""
    with open(codon_input, "r") as input_sequence:
        header = input_sequence.readline() #consume and save the header
        codon = input_sequence.read(3) #open fna file and read 3 bits (3 characters)
        found_start_codon = False
        gc_count = 0
        total_count = 0

        while(True):
            #print("Current codon: " + str(codon))
            for base in codon: #checks the codon for GC content in regards to Q10
                if(base == 'G' or base == 'C'):
                    gc_count+=1
                total_count+=1
                
            if len(codon) != 3: #have reached the end of the file, ignore last 1-2 codons and exit 
                break
            elif not(found_start_codon):
                if(codon == "ATG"): #Reached the beginning of the protein encoding
                    amino_sqn += translate_dict[codon]
                    found_start_codon = True
            else:
                amino_sqn += translate_dict[codon]
                if(codon == "TGA"): #Reached the end of the protein encoding
                    break 
            codon = input_sequence.read(3)

    gc_content = (gc_count/total_count) * 100
    print("GC content for " + header)
    print("-> " + str(gc_content) +"%\n")
    with open(amino_output, "w") as out_sequence:
        out_sequence.write(header)
        out_sequence.write(amino_sqn)
        
        


if __name__ == "__main__":
    if(len(sys.argv) != 4):
        print("Expecting 3 files, (1) a txt file for conversion help (2) a text file for amino acid sequence output (3) codon text file that needs translating")
    else:
        translate(sys.argv[1], sys.argv[2], sys.argv[3])