import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#extract sequences from fna files and store/return in a list
def read_fna_sequences(path):
    #this list will hold all sequences being extracted from the text file  
    sequences = []
    with open(path, 'r') as file:
        for line in file:
            # make sure we not starting with >, not a sequences
            if not ('>' in line):
                # Strip the line clean for \n
                line = line.strip()
                sequences.append(line)
    return sequences

#Calculate conservation rates at each position of sequences
def calculate_conservation(sequences):
    #assume that all sequences are of the same length
    conservation_rates = []
    conservation_bases = []
    # Go through all positions of the sequences  
    for i in range(len(sequences[0])):
        # dynamically build a dictionary that keeps track of the count for each base or '-' encountered
        base_counts = {}
        for cur_seq in sequences:
            base = cur_seq[i]
            # check if the current base or - character is in the dictionary yet
            if base in base_counts:
                #already exists, increment count
                base_counts[base] += 1
            else:
                #doesnt exists
                base_counts[base] = 1

        # Find the maximum count out of the bases (A T C G)
        max_count = -1
        maxbase = None
        for base in base_counts.keys():
            if base_counts[base] > max_count and (base == 'A' or base == 'T' or base == 'C' or base == 'G'):
                max_count = base_counts[base]
                maxbase = base
        # create the conservation rate based on the greatest occuring base
        conservation_rate = max_count/len(sequences) #note length of sequences is synonymous to the total number of bases and '-'
        conservation_rates.append(conservation_rate)
        conservation_bases.append(maxbase)
    
    # uncomment line below for returning non smoothed data and bases for writing out to file for Q1   
    #return conservation_rates, conservation_bases

    # create a data frame w/ pandas, referencing the conservation_rates as 'rate' column
    df = pd.DataFrame(conservation_rates, columns=['rate'])
    # rolling average window of 90 was choosen in order to base smoothing over the another 90 genes positions out of about 1500
    # approximately 6% dependency on other values for current value smoothing, can be increased for an even smoother graph
    # create another column in the dataframe for smoothed values, smoothing is done using pandas 
    df['window_rates'] = df['rate'].rolling(window=90).mean()
    conservation_rates_window = df['window_rates'].tolist()
    
    #Question 3 plotting
    df.plot(y=['window_rates'], kind='line', figsize=(16,12))
    plt.title("Conservation Plot")
    plt.xlabel("Position in 16S RNA gene")
    plt.ylabel("Pct Conserved")
    
    #plt.show()
    return conservation_rates_window, df

#Find variable regions and also threshold for considering a variable region
def variable_regions(data, step_size):
    # Decided to use the average of all smoothed conservation percentages as the threshold point
    #    A differing statistic could be applied here easily but I just thought that the average reflected well 
    #    to the number of variability regions I could depict by eye. Moreover it gives 8 variable regions, only differs by one for the number 
    #    of widely accepted amount of 9.  
    variability_threshold = np.nanmean(data)
    #print("Variability threshold:", variability_threshold)

    variable_regions = []   # holds a tuple of variable regions
    in_region = False       # Depicts if we are currently or just previuly in a variable region (below the threshold)
    start = None            # start coordinate 
    end = None              # end coordinate

    #Go through the list of conserved rates at a step to avoid small 'squigles' in the data and detecting regions not actually variable
    for i in range(0, len(data), step_size):
        # Disregard outside nan values from the list due to the sliding window average
        if np.isnan(data[i]):
            continue
        
        # If the current data point is below our threshold and we are currently not in considered variable region
            # save the start position within the data and notify that we are now currently in a variable region
        if data[i] < variability_threshold and not in_region:
            in_region = True
            start = i
        # If the current data point is greater than or equal to threshold and we are currently in the considered variable region  \
            # save the end position within the data and notify that are now in a non-variable region, above the the threshold
            # we also now have a complete 
        elif data[i] >= variability_threshold and in_region:
            end = i
            variable_regions.append((start, end))
            in_region = False
        #other conditional cases can be disconsidered since they have no effect since they either still in the region or still not in the region

    # Makes sure that if we finish going through the data IN a variable region, we record that end position ~ len(data).
    # save the current start and final (len(data)) position
    if in_region:
        variable_regions.append((start, len(data)))

    return variability_threshold, variable_regions

#Write conservations to an output to a file
def conservations_to_file(conservation_rates, conservation_bases, outfile):
    base_index = 0
    with open(outfile, 'w') as file:
        for rate in conservation_rates:
            #          Conserved base                   Conservation rate
            file.write(f'{conservation_bases[base_index]}: \t{rate}\n')
            base_index+=1

#Write variable regions to an output to a file, tab spaced
def regions_to_file(variable_regions, filename):
    with open(filename, 'w') as file:
        #unpakc tuple of start end regions and write to outfile
        for start, end in variable_regions:
            out_string = str(start) + "\t" + str(end) + "\n"
            file.write(out_string)

if __name__ == "__main__":
    INPUT_FILE = 'Homework4-seqs-with-primers.fna'
    OUTPUT_FILE = 'solution-problem-test.txt'
    OUTPUT_FILE3 = 'solution-problem-test3.txt'
    
    sequences = read_fna_sequences(INPUT_FILE)
    
    #uncomment the 2 line below for Q1 work
    #conservation, bases = calculate_conservation(sequences)
    #conservations_to_file(conservation, bases, OUTPUT_FILE)
    
    conservation, df = calculate_conservation(sequences)
                        #List of conservation values, step size     
    threshold, regions = variable_regions(conservation, 5)
    #print(regions)
    regions_to_file(regions, OUTPUT_FILE3)

    #Question 4 plotting
    df.plot(y=['window_rates'], kind='line', figsize=(16,12))
    plt.title("Conservation Plot")
    plt.xlabel("Position in 16S RNA gene")
    plt.ylabel("Pct Conserved")
    
    #add the threshold as a dashed line 
    first = True
    for start, end in regions:
        if first:
            plt.hlines(y=threshold, xmin=start, xmax=end, colors='r', linewidth=2, label="Variable region threshold")
            first = False
        else:
            plt.hlines(y=threshold, xmin=start, xmax=end, colors='r', linewidth=2)
    
    plt.legend()
    plt.show()
