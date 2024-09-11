import sys
import csv
import os

########
#TODO
## Go through the ct file and store the pairings for each respective index
## Each position/index within the list will store the corresponding index it is paired to
########
def get_pairings(input_ct_file):
    pairings = []
    with open(input_ct_file, 'r') as file:
        first_line = True   #For ignoring the first line
        for line in file:
            if first_line: #ignore the first line of the file
                first_line = False
                continue
            values = line.split()
            pairings.append(int(values[4])-1) #store the index in which the the corresponding index into the list has paired with
    return pairings 
        

########
#TODO
## Go through the ct file and extract the rna sequence that was folded
## Extract consecutive sequences of no pairing as a start and end index of size gap_sizes. list of tuples  
########
def calculate_gaps(input_ct_file, gap_sizes):
    rna_sequence = "" #This is for mainting the sequence of the ct file
    gap_indexes = []  #this is for holding the index positions into the rna sequence in which there are a serious of no pairings

    # Open the ct file and read line by line
        # Column 1: position within the sequence
        # Column 2: base
        # Column 3: index within the sequence 
        # Column 5: what that position is pairing to (looking for 0's) 
    with open(input_ct_file, 'r') as file:
        first_line = True   #For ignoring the first line
        
        gap_start_index = None  #keep track of starting index
        current_gap_size = 0    #keeps track of gap sizes
        gaps = []

        #iterate through all the lines in the ct file
        for line in file:
            
            if first_line: #ignore the first line of the file
                first_line = False
                continue
            values = line.split()
            
            #build out the sequence 
            rna_sequence += values[1]

            # encounter a non-paired base gap 
            if values[4] == '0':
                # currently not building a range of gap
                if gap_start_index == None:
                    gap_start_index = int(values[2]) #set starting positon of building a gap
                current_gap_size += 1
            
            # encounter a paired base
            else:
                # end building range of gap
                if gap_start_index != None:
                    # found a consecutive set of gaps of at least the size we are looking for 
                    if current_gap_size >= gap_sizes:
                        #           Start(inclusive)  End(exclusive)  
                        gaps.append((gap_start_index, int(values[2])))
                    # reset gap tracking values
                    current_gap_size = 0
                    gap_start_index = None
    #print(rna_sequence)
    #print(gaps)
    return rna_sequence, gaps

########
#TODO
## Given rna sequence and no paired gaps within the rna sequnce, find pseudoknots
## Return a list of pair ranges of indexes that form possible psuedoknots 
########
def find_knots(rna_sequence, gap_ranges, minSizeKnot):
    knots = [] # list of pairings
    for i in range(len(gap_ranges)-1):
        for j in range(i+1, len(gap_ranges)):
            # compare 2 sequences to find a possible pseudonot within the rna sequence
            longerRange, shorterRange = longer_sequence(gap_ranges[i], gap_ranges[j])
            # calculate number of frameshifts are needed to fully compare the 2 gap sequences
            frameShifts = (longerRange[1]-longerRange[0]) - (shorterRange[1] - shorterRange[0]) + 1
            for frameshift in range(frameShifts):
                knots += (knots_in_sequences(rna_sequence, (longerRange[0] + frameshift, longerRange[1] - (frameShifts - frameshift - 1)), shorterRange, minSizeKnot))
    return knots

# Helper function for find_knots for a given frameshift
#   takes in 2 ranges of unpaired sequences in rna and checks to see if a knot would form between them
#assumes that the ranges are of the same size!
def knots_in_sequences(sequence, range1, range2, minSizeKnot):
    knots = []
    range1Iterator = range1[0]
    range2Iterator = range2[1]-1
    
    buildingKnot = False
    knotSize = 0
    range1Start = None
    range2Start = None
    range1End = None
    range2End = None
    #loop through range1 ASCENDINGLY and loop through range2 DESCENDINGLY to allow for inverse comparing
    while range1Iterator < range1[1] and range2Iterator >= range2[0]:
        isPair = correct_base_pairing(sequence[range1Iterator], sequence[range2Iterator])

        #found a pair and no knot is buing built
        if(isPair and (not buildingKnot)):
            range1Start = range1Iterator
            range2End = range2Iterator + 1
            knotSize += 1
            buildingKnot = True
            
        #found a pair and currently building knot
        elif(isPair and buildingKnot):
            knotSize += 1
        
        #did not find a pair and no knot is being built (nothing to do)
        elif((not isPair) and (not buildingKnot)):
            pass

        #did not find a pair and currently building knot
        elif((not isPair) and buildingKnot):
            range1End = range1Iterator
            range2Start = range2Iterator + 1
            buildingKnot = False
            if(knotSize >= minSizeKnot):
                if range1Start < range2Start:
                    knots.append(((range1Start, range1End), (range2Start, range2End), knotSize))
                else:
                    knots.append(((range2Start, range2End), (range1Start, range1End), knotSize))
            knotSize = 0
        
        #move to next bases
        range1Iterator+=1
        range2Iterator-=1

    # finish knot that was being built
    if(buildingKnot):
        range1End = range1Iterator
        range2Start = range2Iterator + 1
        buildingKnot = False
        if(knotSize >= minSizeKnot):
            if range1Start < range2Start:
                knots.append(((range1Start, range1End), (range2Start, range2End), knotSize))
            else:
                knots.append(((range2Start, range2End), (range1Start, range1End), knotSize))
    return knots

def longer_sequence(range1, range2):
    if range1[1] - range1[0] >= range2[1] - range2[0]:
        return range1, range2
    else:
        return range2, range1

def correct_base_pairing(base1, base2):
    return ((base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A')) or ((base1 == 'G' and base2 == 'C') or (base1 == 'C' and base2 == 'G'))

def print_gaps(sequence, gaps):
    for gap in gaps:
        print(str(gap[0]) + "-" + str(gap[1]) + ": " + sequence[gap[0]:gap[1]])

def print_knots(sequence, knots):
    for knot in knots:
        print(knot)
        print("Sequence 1: " + str(sequence[knot[0][0]: knot[0][1]]))
        print("Sequence 2: " + str(sequence[knot[1][0]: knot[1][1]]))
    print()

def write_knots_to_csv(knots, output_file, filepath):
    # Check if the file exists and if it's empty
    file_exists = os.path.isfile(output_file)
    file_empty = True if not file_exists or os.stat(output_file).st_size == 0 else False

    # Open the output file in append mode
    with open(output_file, 'a', newline='') as file:
        writer = csv.writer(file)

        # Write the header only if the file is empty
        if file_empty:
            writer.writerow(['Range1 Start', 'Range1 End', 'Range2 Start', 'Range2 End', 'Knot Size', 'ct file/species'])
        writer.writerow(['', '', '', '', '', filepath])
        # Write the data
        for knot in knots:
            range1, range2, knot_size = knot
            writer.writerow([range1[0], range1[1], range2[0], range2[1], knot_size])


# Produce a visual representation of possible pseudo-knots and also what pairing/bonding is formed before the found regions of pseudo knot
def visualize_knots(sequence, knots, pairings, offset):
    # *note: offset determines how many bases are needed to show for extensive matching besides the pseudo knot range
    for knot in knots:
        #unpack ranges
        range1, range2, knot_size = knot
        range1Start = range1[0]
        range1End = range1[1]
        range2Start = range2[0]
        range2End = range2[1]

        adjustedRange1Start = max(range1Start - offset, 0)
        adjustedRange2Start = max(range2Start - offset, 0)
        
        if(adjustedRange1Start == 0):
            range1offset = 0 - range1Start - offset
        else:
            range1offset = offset
        
        if(adjustedRange2Start == 0):
            range2offset = 0 - range2Start - offset
        else:
            range2offset = offset

        range1CounterpartPairing = "" # row 1
        range1CounterpartPairingSymbols = "" #row 2
        for i in range(adjustedRange1Start, range1Start):
            if(pairings[i] == -1): #case in which there is no pairing for the counterpart
                range1CounterpartPairing += "-"
                range1CounterpartPairingSymbols += " "
            else:
                range1CounterpartPairing += sequence[pairings[i]]
                range1CounterpartPairingSymbols += "|"
        

        range1String = sequence[adjustedRange1Start:range1End] # row 3
        
        range1range2Symbols = "" # row 4
        for i in range(range1offset):
            range1range2Symbols += " "
        for i in range(knot_size):
            range1range2Symbols += "|" # can assume pairings since we found the knots based on pairings
        
        range2String = sequence[adjustedRange2Start:range2End] # row 5
        for i in range(range2offset):
            range2String += " "
        range2String = range2String[::-1] # reverse row 5 for lining up with row 3

        range2CounterpartPairingSymbols = "" #row 6
        range2CounterpartPairing = "" # row 7
        
        for i in range(adjustedRange2Start, range2Start):
            if(pairings[i] == -1): #case in which there is no pairing for the counterpart
                range2CounterpartPairing += "-"
                range2CounterpartPairingSymbols += " "
            else:
                range2CounterpartPairing += sequence[pairings[i]]
                range2CounterpartPairingSymbols += "|"
        for i in range(range1offset+knot_size):
            range2CounterpartPairing += " "
            range2CounterpartPairingSymbols += " "
        
        range2CounterpartPairing = range2CounterpartPairing[::-1]
        range2CounterpartPairingSymbols = range2CounterpartPairingSymbols[::-1]

        print("Knot range(s): " + str(range1Start) + "-" + str(range1End) + " with " + str(range2Start) + "-" + str(range2End))
        # print built rows for visualization
        print(range1CounterpartPairing)
        print(range1CounterpartPairingSymbols)
        print(range1String)
        print(range1range2Symbols)
        print(range2String)
        print(range2CounterpartPairingSymbols)
        print(range2CounterpartPairing)
        print()
            
# main for running calculate_gaps, find_knots, and visualize_knots for finding possible pseudoknots in an rna sequence
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Expected inputs (1) CT file for extracting pairing gaps and producing pseudo knots (2) minimum gap size we are looking for")
    else:
        sequence, gaps = calculate_gaps(sys.argv[1], int(sys.argv[2]))
        #print_gaps(sequence, gaps)
        knots = find_knots(sequence, gaps, 5) # a tuple holding a pair of tuples that correspond to the range of the indexes that produce a potential knot.
        print_knots(sequence, knots)
        pairings = get_pairings(sys.argv[1])
        visualize_knots(sequence, knots, pairings, 20)
        #write_knots_to_csv(knots, "knots.csv", sys.argv[1])

        
##############
# PAST ITERATIONS 
##############
#original working visualizer, doesnt take into account index out of range for, indeces < 0
def visualize_knots_2(sequence, knots, pairings, offset):
    # *note: offset determines how many bases are needed to show for extensive matching besides the pseudo knot range
    for knot in knots:
        #unpack ranges
        range1, range2, knot_size = knot
        range1Start = range1[0]
        range1End = range1[1]

        range2Start = range2[0]
        range2End = range2[1]

        range1CounterpartPairing = "" # row 1
        range1CounterpartPairingSymbols = "" #row 2
        for i in range(range1Start-offset, range1Start):
            if(pairings[i] == -1): #case in which there is no pairing for the counterpart
                range1CounterpartPairing += "-"
                range1CounterpartPairingSymbols += " "
            else:
                range1CounterpartPairing += sequence[pairings[i]]
                range1CounterpartPairingSymbols += "|"
        

        range1String = sequence[range1Start-offset:range1End] # row 3
        
        range1range2Symbols = "" # row 4
        for i in range(offset):
            range1range2Symbols += " "
        for i in range(knot_size):
            range1range2Symbols += "|" # can assume pairings since we found the knots based on pairings
        
        range2String = sequence[range2Start-offset:range2End] # row 5
        for i in range(offset):
            range2String += " "
        range2String = range2String[::-1] # reverse row 5 for lining up with row 3

        range2CounterpartPairingSymbols = "" #row 6
        range2CounterpartPairing = "" # row 7
        
        for i in range(range2Start-offset, range2Start):
            if(pairings[i] == -1): #case in which there is no pairing for the counterpart
                range2CounterpartPairing += "-"
                range2CounterpartPairingSymbols += " "
            else:
                range2CounterpartPairing += sequence[pairings[i]]
                range2CounterpartPairingSymbols += "|"
        for i in range(offset+knot_size):
            range2CounterpartPairing += " "
            range2CounterpartPairingSymbols += " "
        
        range2CounterpartPairing = range2CounterpartPairing[::-1]
        range2CounterpartPairingSymbols = range2CounterpartPairingSymbols[::-1]

        print("Knot range(s): " + str(range1Start) + "-" + str(range1End) + " with " + str(range2Start) + "-" + str(range2End))
        print(range1CounterpartPairing)
        print(range1CounterpartPairingSymbols)
        print(range1String)
        print(range1range2Symbols)
        print(range2String)
        print(range2CounterpartPairingSymbols)
        print(range2CounterpartPairing)
        print()



################################################################################################################
### previous solution for finding knots, not correct
################################################################################################################

#TODO
## Given rna sequence and no paired gaps within the rna sequnce, find pseudoknots
## Return a list of pair ranges of indexes that form possible psuedoknots 
########
#EXTENSION, if the pseudoknot is not at least k bases long, dont recognize it? 
def finding_knots(rna_sequence, gap_ranges):
    knots = [] # list of pairings

    for i in range(len(gap_ranges)-1):
        for j in range(i, len(gap_ranges)):
            # compare 2 sequences to find a possible pseudonot within the rna sequence
            longerRange, shorterRange = longer_sequence(gap_ranges[i], gap_ranges[j])
            # calculate number of frameshifts are needed to fully compare the 2 gap sequences
            frameShifts = (longerRange[1]-longerRange[0]) - (shorterRange[1] - shorterRange[0]) + 1
            
            #Keep track of the longest ranges
            maxNumberOfPairings = -1
            maxShorterRangePairings = None
            maxLongerRangePairings = None

            for frameshift in range(frameShifts):
                #iterate through the other range indexes in reverse to get alignment replicating a pseudoknot
                longRangeStartIndex = frameshift + longerRange[0]
                longRangeIterator = frameshift + longerRange[0]
                shortRangeEndIndex = shorterRange[1]
                numberOfPairings = 0  
                for shortRangeIterator in range(shorterRange[1]-1, shorterRange[0]-1, -1):
                    # if match, increment tally for number of matches
                    if correct_base_pairing(rna_sequence[shortRangeIterator], rna_sequence[longRangeIterator]):
                       numberOfPairings += 1
                       longRangeIterator+=1
                    # end comparing sequences for that frameshift
                    else:
                        break
                #check to see if we found a new longest pairing in between 2 sequence ranges
                if numberOfPairings > maxNumberOfPairings and numberOfPairings > 0:
                    maxNumberOfPairings = numberOfPairings
                    maxLongerRangePairings = (longRangeStartIndex, longRangeStartIndex+numberOfPairings)
                    maxShorterRangePairings = (shortRangeEndIndex-numberOfPairings, shortRangeEndIndex)
            
            # add possible found pseudo knot to the list 
            if(maxLongerRangePairings != None and maxShorterRangePairings != None and maxNumberOfPairings != -1):
                knots.append((maxLongerRangePairings, maxShorterRangePairings, maxNumberOfPairings))

    for knot in knots:
        print(knot)
    return knots