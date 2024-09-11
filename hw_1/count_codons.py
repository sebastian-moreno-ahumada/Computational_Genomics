#This program takes in a fna file and a CSV
#Counts the number of respective codon patterns from fna file and writes to CSV

import sys
import csv
import re
def count_codons(fna_file, csv_file):
    codon_pattern_count = {} #Dictionary for storing the number of appearences for a codon pattern
    pattern = r"^[ACGT]{3}$" #regex for detecting non-codon patterns
    try:
        fna = open(fna_file, "r")
        codon = fna.read(3) #open fna file and read 3 bits (3 characters)
        while(True):
            if len(codon) != 3: #have reached the end of the file, ignore last 1-2 codons and exit 
                break
            elif not(re.match(pattern, codon)): #if we encounter buffer info, not codons, we read the line. applies to seperate genomes misc info.
                codon = fna.read(3)
                continue
            elif codon in codon_pattern_count:  #if we encounter a codon pattern previously encountered, increment the count
                codon_pattern_count[codon] += 1
            else:
                codon_pattern_count[codon] = 1 #if we encounter a codon pattern not previously enountered, introduce it to dictionary
            codon = fna.read(3) #Move on to the next codon/base pattern.
        fna.close()
    except FileNotFoundError:
        print("fna file not found")
        return
    c = open(csv_file, "w", newline="") #Open CSV file for writing
    csv_write = csv.writer(c) # Use CSV lib for easy CSV writing
    for k, v in codon_pattern_count.items():
        csv_write.writerow([k, v]) # write dictionary information into CSV file with CSV Writer
    c.close()
    return

if __name__ == "__main__": #Main that takes in terminal arguments and calls count_codons based on arguments 
    if(len(sys.argv) != 3):
        print("Expecting 2 arguments: (1) FNA file containing codons/bases (2) CSV file for writing\n")
    else:
        count_codons(sys.argv[1], sys.argv[2])