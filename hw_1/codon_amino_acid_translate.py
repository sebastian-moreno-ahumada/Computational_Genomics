import sys
import csv
def translate(table_file, csv_output, csv_input):
    translate_dict = {} #Key: codon, Value: amino acid, for translating codons to amino acids
    with open(table_file, "r") as tf:
        tf.readline(); #consume firstline
        cur_line = tf.readline()
        while(cur_line != ""):
            line_arr = cur_line.strip().split('\t')
            codon = line_arr[0]
            amino = line_arr[1]
            translate_dict[codon] = amino
            cur_line = tf.readline()
    #translation dictionary complete

    #begin creating an amino acid dictionary w occurrences based on translations.
    amino_dict = {} #Key: amino acid, Value: # of Occurrences
    with open(csv_input, "r") as csv_in:
        cur_in_line = csv_in.readline()
        while(cur_in_line != ""):
            cur_in_line_arr = cur_in_line.strip().split(',')
            current_codon = cur_in_line_arr[0]
            current_codon_amnt = cur_in_line_arr[1]
            if translate_dict[current_codon] in amino_dict: #the amino acid already has a codon(s) and its(their) occurnces
                amino_dict[translate_dict[current_codon]] += int(current_codon_amnt)
            else:
                amino_dict[translate_dict[current_codon]] = int(current_codon_amnt)
            cur_in_line = csv_in.readline()
    #Start writing amino acid totals to CSV for output.
    c = open(csv_output, "w")
    csv_write = csv.writer(c)
    for k, v in amino_dict.items():
        csv_write.writerow([k, v])
    c.close()
    return

if __name__ == "__main__":
    if(len(sys.argv) != 4):
        print("Expecting 2 files, (1) a txt file for conversion help (2) a csv file for amino acid output (3) codon file that needs translating")
    else:
        translate(sys.argv[1], sys.argv[2], sys.argv[3])