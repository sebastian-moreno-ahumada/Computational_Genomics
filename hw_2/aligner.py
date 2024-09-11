import argparse
#NOTE
###
# I discussed programming solutions with Galaan Omar (Omar0243)
#
# I also used AI (ChatGPT LLM) for...
# - high level pseudocode/process explanation for NW algorithm
# - The following specific chunk of code/logic
# - // Traceback to find the optimal alignment
Traceback_Logic = """
alignment_sequence1 = ""
alignment_sequence2 = ""
visualization = ""
i = length(sequence1)
j = length(sequence2)

while i > 0 or j > 0:
    current_score = score_matrix[i][j]
    if i > 0 and score_matrix[i-1][j] + gap_penalty == current_score:
        // Move upwards (gap in sequence1)
        alignment_sequence1 = sequence1[i-1] + alignment_sequence1
        alignment_sequence2 = "_" + alignment_sequence2
        visualization = " " + visualization
        i = i - 1
    else if j > 0 and score_matrix[i][j-1] + gap_penalty == current_score:
        // Move to the left (gap in sequence2)
        alignment_sequence1 = "_" + alignment_sequence1
        alignment_sequence2 = sequence2[j-1] + alignment_sequence2
        visualization = " " + visualization
        j = j - 1
    else:
        // Diagonal move (match/mismatch)
        alignment_sequence1 = sequence1[i-1] + alignment_sequence1
        alignment_sequence2 = sequence2[j-1] + alignment_sequence2
        
        // Check for match or mismatch
        if sequence1[i-1] == sequence2[j-1]:
            visualization = "|" + visualization
        else:
            visualization = "X" + visualization

        i = i - 1
        j = j - 1
"""
#Ultimately the traceback algorithm was helpful, but required significant modification; can be seen during the traceback portion of my code

#returns None on error 
def needleman_wunsch(sequence1_file, sequence2_file, outputFile, noStartEndPenalty, gapPenalty, mismatchPenalty, matchScore):
    sequence1 = ""
    sequence2 = ""
    sequence1_header = ""
    sequence2_header = ""
    try:
        with open(sequence1_file, 'r') as sequence:
            sequence1_header = sequence.readline()#consume the first line for the title
            sequence1 = sequence.readline().strip()
    except FileNotFoundError:
        print("Incorrect file path for sequence 1\n")
        return None
    try:
        with open(sequence2_file, 'r') as sequence:
            sequence2_header = sequence.readline()#consume the first line for the title
            sequence2 = sequence.readline().strip()
    except FileNotFoundError:
        print("Incorrect file path for sequence 2\n")
        return None
    #sequence on side/vertically
    n = len(sequence1)
    #print("Sequnce1 length: " + str(n))
    #sequnce on top/horizontally
    m = len(sequence2)
    #print("Sequnce2 length: " + str(m))
    #initialize 2-d matrix for scoring w/ 0
    matrix_score = []
    for i in range(n+1):
        row = []
        for j in range(m+1):
            row.append(0)
        matrix_score.append(row)
    #initialize first row and column to account for gap penalty 
    if not(noStartEndPenalty):
        for i in range(len(matrix_score)):
            matrix_score[i][0] = i*gapPenalty
        for j in range(len(matrix_score[0])):
            matrix_score[0][j] = j*gapPenalty
    #fill in score_matrix with optimal scoring for all the cells based on previous cells
    for i in range(1, n+1):
        for j in range(1, m+1):
            #calculate the running score for a match/mismatch
            if(sequence1[i-1] == sequence2[j-1]): # We have a match
                match_mismatch_score = matrix_score[i-1][j-1] + matchScore
            else:                                               # We have a mismatch
                match_mismatch_score = matrix_score[i-1][j-1] + mismatchPenalty
            
            #calculate the running score for a gap from either sequence1 or sequence2
            #take into consideration of no start or end gap penalties
            if noStartEndPenalty and i==n:#case in which seq1 has been fully consumed and want to add gaps at the end of it (no penalty)
                seq1Gap = matrix_score[i][j-1]
            else:
                seq1Gap = matrix_score[i][j-1] + gapPenalty
            
            if noStartEndPenalty and j==m:#case in which seq2 has been fully consumed and want to add gaps at the end of it (no penalty)
                seq2Gap = matrix_score[i-1][j]
            else:
                seq2Gap = matrix_score[i-1][j] + gapPenalty

            matrix_score[i][j] = max(match_mismatch_score ,seq1Gap, seq2Gap)
    # double check score_matrix computation here
    #for row in matrix_score:
    #    print(row)

    #perform traceback to find the optimal score 
    #start at very bottom right
    i = n
    j = m
    optimalScore = matrix_score[i][j]
    alignment_sequence1 = ""
    alignment_sequence2 = ""
    visualization = ""
    mismatch_count = 0
    while(i > 0) or (j > 0):
        if noStartEndPenalty: #no start end penalty alignment is being incorporated
            #Go up, source is 1 up from current position: consume from sequence1 and put gaps in sequence2
            if i > 0 and (((j==m or j==0) and matrix_score[i-1][j] == matrix_score[i][j]) or (matrix_score[i-1][j] + gapPenalty == matrix_score[i][j])):
                alignment_sequence1 = sequence1[i-1] + alignment_sequence1
                alignment_sequence2 = "_" + alignment_sequence2
                visualization = " " + visualization
                i=i-1
            #Go left, source is 1 left from current position: consume from sequence2 and put gaps in sequence1
            elif j > 0 and (((i==n or i==0) and matrix_score[i][j-1] == matrix_score[i][j]) or (matrix_score[i][j-1] + gapPenalty == matrix_score[i][j])):
                alignment_sequence1 = "_" + alignment_sequence1
                alignment_sequence2 = sequence2[j-1] + alignment_sequence2
                visualization = " " + visualization
                j=j-1
            #Go diagonal, source is up-left diagonally: consume from sequence1 and sequence2. update visualizer respectively
            
            else: #dangerous default value?
                alignment_sequence1 = sequence1[i-1] + alignment_sequence1
                alignment_sequence2 = sequence2[j-1] + alignment_sequence2
                #visualizer update based on match or mismatch
                if sequence1[i-1] != sequence2[j-1]:
                    visualization = "x" + visualization
                    mismatch_count+=1
                else:
                    visualization = "|" + visualization
                i=i-1
                j=j-1

        else: #start end penalty alignment is being incorporated
            #Go up, source is 1 up from current position: consume from sequence1 and put gaps in sequence2
            if i > 0 and (matrix_score[i-1][j] + gapPenalty == matrix_score[i][j]):
                alignment_sequence1 = sequence1[i-1] + alignment_sequence1
                alignment_sequence2 = "_" + alignment_sequence2
                visualization = " " + visualization
                i=i-1
            #Go left, source is 1 left from current position: consume from sequence2 and put gaps in sequence1
            elif j > 0 and (matrix_score[i][j-1] + gapPenalty == matrix_score[i][j]):
                alignment_sequence1 = "_" + alignment_sequence1
                alignment_sequence2 = sequence2[j-1] + alignment_sequence2
                visualization = " " + visualization
                j=j-1

            else:
            #Go diagonal, source is up-left diagonally: consume from sequence1 and sequence2. update visualizer respectively
            #dangerous default value?
                alignment_sequence1 = sequence1[i-1] + alignment_sequence1
                alignment_sequence2 = sequence2[j-1] + alignment_sequence2
                #visualizer update based on match or mismatch
                if sequence1[i-1] != sequence2[j-1]:
                    visualization = "x" + visualization
                    mismatch_count+=1 
                else:
                    visualization = "|" + visualization
                i=i-1
                j=j-1

    #begin writting respective information to the output file
    #Add a flag for outputting mismatch count?
    with open(outputFile, "w") as outfile:
        outfile.write(str(optimalScore) + "\n")
        outfile.write(sequence1_header)
        outfile.write(alignment_sequence1 + "\n")
        outfile.write(visualization + "\n")
        outfile.write(alignment_sequence2 + "\n")
        outfile.write(sequence2_header) 
        # outfile.write("Mismatch count: " + str(mismatch_count) + "\n") - was for counting the number of mismatches for Q6
        

if __name__ == "__main__":
    #main method to setup function call to 
    parser = argparse.ArgumentParser(description='Needleman-Wunsch Alignment')
    # set up argument parser
    parser.add_argument('-q', '--query', required=True, type=str, help='FNA file containing query sequence')
    parser.add_argument('-r', '--reference', required=True, type=str, help='FNA file containing reference sequence')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output file for optimal alignment information and visualization')
    parser.add_argument('-g', '--gap_penalty', required=True, type=int, help='Gap penalty')
    parser.add_argument('-p', '--mismatch_penalty', required=True, type=int, help='Mismatch penalty')
    parser.add_argument('-m', '--match_score', required=True, type=int, help='Match score')
    parser.add_argument('--ignore_outer_gaps', action='store_true', help='Do not include start or end gap penalties')
    
    args = parser.parse_args()
    needleman_wunsch(args.query, args.reference, args.output, args.ignore_outer_gaps, args.gap_penalty, args.mismatch_penalty, args.match_score)

    #previous command line parser
    og_CL_parser = """
    if(len(sys.argv) < 4):
        #NOTE, include extra possible parameters for user customization of penalties and alignment scoring?
        ("Expecting arguments: (1)FNA file containing sequence I (2)FNA file containing sequence II (3)Output file for optimal alignment information and visualization (4 optional) -NP *for no start-end gap penalties* \n")
    else:
        noStartEndPenalty = False
        for arg in sys.argv:
            if arg == "-NP":
                print("Not including start or end penalties")
                noStartEndPenalty = True
        needleman_wunsch(sys.argv[1], sys.argv[2], sys.argv[3], noStartEndPenalty, -2, -1, 1)
"""