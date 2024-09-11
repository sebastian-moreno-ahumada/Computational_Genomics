import sys
import os
import random
from skbio.tree import TreeNode as tn

#import pdb; pdb.set_trace()

###################################################################################
# Code for bootstrapping
###################################################################################

#Q6 steps *for my sake of understanding*
# make a new folder
# create a base ID to Sequence dictionary, similar to the one in distances
# for loop n times, base sequence files and tree files on current iteration
    # Make a new dictionary for the new base ID to Sequence dictionary

#####
#TODO
# take in a root of a tree and produce and get all the children tips NAME
def partitioning(beginningTreeNode):
    tips = []
    for n in beginningTreeNode.preorder():
        if n.children == []:
            tips.append(n.name)
    return tips

####
# TODO
# create bootstrap fna files, bootstrap trees, and calculate support percentages and output to bootstrap.txt
def bootstrap(input_file, numberOfBootstraps, originalTreeRoot):
    #Begin Q6
    FOLDER_NAME_BOOTSTRAP = "Bootstrap_files"
    #make a list with all the sequences, similar to the one used for question 2 (*should probably use a method for coupling...)
    baseIdSeqDict = {}
    try:
        with open(input_file, 'r') as f:
            #loop through the fna file until reach an empty line
            lengthOfSequences = None  #represent the global length of the sequences 
            while True:
                lineID = f.readline() #read the first line for the ID
                if not lineID: # If the line is empty, we've reached the end of the file
                    break
                lineSeq = f.readline() #read the second line for the Sequence
                #clean and strip the the lines of any whitespace
                lineID = lineID.strip()
                lineID = lineID[1:]         # Remove the leading '>'
                lineSeq = lineSeq.strip()
                baseIdSeqDict[lineID] = lineSeq
                lengthOfSequences = len(lineSeq)
    except FileNotFoundError:
        print("The specified file was not found, unable to create output file(s).")
        return
    
    # Create a new directory if it doesn't exist
    if not os.path.exists(FOLDER_NAME_BOOTSTRAP):
        os.makedirs(FOLDER_NAME_BOOTSTRAP)

    # Maintain a list of all the roots of bootstrap trees
    bootstrapTrees = []

    #Create n bootstrap distance matrices and trees
    for i in range(numberOfBootstraps):
        
        #create a differing dictionary for bootstrapped sequences, fill with existing ID's and set sequences to empty
        otherIdSeqDict = {}
        for ID in baseIdSeqDict.keys():
            otherIdSeqDict[ID] = ""

        for j in range(lengthOfSequences): # Need m characters from randomColumnIndex to build new sequence
            # Generate a random integer between 0 and greatest possible index into sequences
            randomColumnIndex = random.randint(0, lengthOfSequences-1)
            # iterate through the list of keys/ID's and append their respective random column character
            for currentID in baseIdSeqDict.keys():
                # Updated sequence          # Existing sequence         # Base respective to sequence and random column chosen
                otherIdSeqDict[currentID] = otherIdSeqDict[currentID] + baseIdSeqDict[currentID][randomColumnIndex]

        #Completed making new bootstrap sequence, output to file with fna structure and then create a tree
        FOLDER_NAME_BOOTSTRAP = "Bootstrap_files"
        outputFileSequences = os.path.join(FOLDER_NAME_BOOTSTRAP, f"sequences-bootstrap_{i+1}.txt")
        outputFileTree = os.path.join(FOLDER_NAME_BOOTSTRAP, f"edges-bootstrap_{i+1}.txt")

        # create a bootstrap fna file
        idSequenceString = ""
        for currentID in otherIdSeqDict.keys(): #go through newly built sequences for each ID and output to file in fna format
            idSequenceString += ">" + str(currentID) + "\n"
            idSequenceString += otherIdSeqDict[currentID] + "\n"

        # Open the an output file for writing the newly computed bootstrap sequences
        with open(outputFileSequences, 'w') as file:
            file.write(idSequenceString)
        
        # Build a tree, give it the boolean flags for not writting original distance matrix and updated output file name
        # also save the built bootstrap tree to the list of bootstrap trees
        bootstrapTrees.append(geneticDistance(outputFileSequences, outputFileTree, False))
    
    #Begin Q7
    # Outer loop: Go through all internal nodes in the main tree, starting at the root...
    supportPercentages = []
    for node in originalTreeRoot.preorder():
        # Do not perform any partitioning if we are at a tip in the main tree...
        if node.children == []:
            # Continue on to the next node
            continue
        
        mainPartitionSet = []
        nodeSupportTally = 0
        # Boolean to help handle root or other internal nodes
        twoChildren = True 
        if len(node.children) == 3:
            twoChildren = False
        #calculate the current partition set of the node we are on in the main tree
        mainPartitionSet.append((partitioning(node.children[0])))
        mainPartitionSet.append((partitioning(node.children[1])))
        if not twoChildren:
            mainPartitionSet.append(partitioning(node.children[2]))
        
        # sort the partitions for easy comparing later
        for partitionIndex in range(len(mainPartitionSet)):
            mainPartitionSet[partitionIndex].sort()
        
        # Loop: For each tree in the list of bootstrap trees...
        for tree in bootstrapTrees:
            # Inner loop, Go through all internal nodes of the currrently selected bootstrap tree
            for otherNode in tree.preorder():

                if otherNode.children == []:
                    # Continue on to the next node
                    continue

                # create the partition for the current bootstrap node
                otherParitionSet = []
                # at the root for the main tree
                if not twoChildren:
                    # case in which we are at the root for both trees
                    if len(otherNode.children) == 3:
                        otherParitionSet.append(partitioning(otherNode.children[0]))
                        otherParitionSet.append(partitioning(otherNode.children[1]))
                        otherParitionSet.append(partitioning(otherNode.children[2]))
                    else:
                        #No need to continue, have already passed the root case for the given bootstrap
                        break
                
                # Not at the root for the main tree 
                else:
                    #at the root for the bootstrap tree, 
                    if len(otherNode.children) == 3:
                        # Move on to the next node within the bootstrap tree
                        continue
                    # case in which we are at internal nodes for both trees
                    else:
                        otherParitionSet.append(partitioning(otherNode.children[0]))
                        otherParitionSet.append(partitioning(otherNode.children[1]))
                
                # sort the partitions for easy comparing later
                for partitionIndex in range(len(otherParitionSet)):
                    otherParitionSet[partitionIndex].sort()
                
                # compare the partitions 
                # case for which we are at the root for both trees
                # *both mainParitionSet and otherParitionSet should have the same size at this point*
                
                #print("Main Partition: ")
                #print(mainPartitionSet)
                #print("Other Partition: ")
                #print(otherParitionSet)

                correctTally = 0
                for partition in mainPartitionSet:
                    for otherPartition in otherParitionSet:
                        if partition == otherPartition:
                            correctTally += 1
                            break
                
                if not twoChildren and correctTally == 3: # Found equivelant partitions for an internal node in the main tree
                    nodeSupportTally += 1
                    break # move on to the next bootstrap tree for partition matching
                elif twoChildren and correctTally == 2:
                    nodeSupportTally += 1
                    break # move on to the next bootstrap tree for partition matching
                    
        
        # Once we've gone through all the trees, append a tuple to a list containing current internal node name and its support percentage
        supportPercentages.append((node.name, (nodeSupportTally/float(numberOfBootstraps))))
    
    # Calculated all support percentages, write out to file 
    SUPPORT_FILE_NAME = "bootstrap.txt"
    with open(SUPPORT_FILE_NAME, 'w') as fileHandle:
        for currSupport in supportPercentages:
            strBuild = str(currSupport[0]) + "\t" + str(currSupport[1]) + "\n"
            fileHandle.write(strBuild)


###################################################################################
# Code for distances
###################################################################################
#TODO
# Take in a file to output tree edges based off of a pre-built tree given a root
# Note, the the given names for the internal nodes will not matter, only the names to the internal tips
def writeTreeToFile(outfile, root, internalEdgeNumbering, IdToOriginalIdx):
    edgesString = "" # String for building output
    edgesFile = open(outfile, "w")
    for n in root.preorder():
        if n == root: #root doesnt have a parent, dont need to print but reassign its ID to #tips + 1
            n.name = str(internalEdgeNumbering)
            internalEdgeNumbering += 1
        else: #not at a root
            if(n.children != []): # Internal node that is not the root and is not tip, need to give it a new ID based on pre-order
                n.name = str(internalEdgeNumbering)
                internalEdgeNumbering+=1
                #               Ancestor Node             Descendant Node       #Lenght in between
                edgesString += str(n.parent.name) + "\t" + str(n.name) + "\t" + str(n.length) + "\n"
            else: #We are at a tip
                n.name = IdToOriginalIdx[n.name] #Translate original ID to the sequences order within fna
                edgesString += str(n.parent.name) + "\t" + str(n.name) + "\t" + str(n.length) + "\n"

    edgesFile.write(edgesString)

#TODO 
# Create an initial list of tip-nodes with correct ID names
def originalTips(tipsIDs):
    tips = []
    for ID in tipsIDs:
        tips.append(tn(ID))
    return tips

#TODO:
#prints a matrix for debugging
def mprint(matrix):
    for row in matrix:
        print(row)

###########
#TODO: 
#BUILDS A disimilarity/distance matrix, also output the file into a matrix
################################################
def makeD(input_file, outputToAFile=True):
    OUTPUT_FILE="genetic-distances.txt"
    #make a dictionary with all the sequences
    IdSeqDict = {}
    sequenceOrdering = {}
    orderingCounter = 1
    try:
        with open(input_file, 'r') as f:
            #loop through the fna file until reach an empty line
            while True:
                lineID = f.readline() #read the first line for the ID
                if not lineID: # If the line is empty, we've reached the end of the file
                    break
                lineSeq = f.readline() #read the second line for the Sequence
                #clean and strip the the lines of any whitespace
                lineID = lineID.strip()
                lineID = lineID[1:]         # Remove the leading '>'
                lineSeq = lineSeq.strip()
                IdSeqDict[lineID] = lineSeq
                #For writing to tree to file later for changing sequenceID to their respective ordering in the fna
                sequenceOrdering[lineID] = str(orderingCounter)
                orderingCounter += 1
    except FileNotFoundError:
        print("The specified file was not found, unable to create output file.")
        return
    #Now that we have a dictionary of all the string for comparing, iterate through them all and compute genetic distances
    #First, create the header row string  
    headerRowString = ""
    taxonOrderingID = list(IdSeqDict.keys()) #maintains a list of the respective sequences with their respective ID's
    for id_index in range(len(taxonOrderingID)):
        headerRowString += taxonOrderingID[id_index]
        if(id_index == len(taxonOrderingID)-1): # when at last sequence, new line
            headerRowString += "\n"
        else:
            headerRowString += "\t"

    #Secondly, iterate through all the strings in the list, calculate their dissimilarity, and write to output file.
        #Create a 2-d matrix to represent the respective dissimilarity values 
    tableString = "" 
    disimilarityPercentage = None
    disimilarityTable = []                              #build inital 2-d distance matrix
    for i in range(len(taxonOrderingID)):
        tableString += str(taxonOrderingID[i]) + "\t"    # String representation for a certain sequnces distances to others
        disimilarityRow = []
        for j in range(len(taxonOrderingID)):                       #compare strings of sequences
            if(i == j):
               #dont calculate dissimilarity
               disimilarityPercentage = 0
            else:
                #calculate dissimilarity
                totalCounter = 0
                disimilarityCounter = 0
                currentSequence = IdSeqDict[taxonOrderingID[i]]
                otherSequence = IdSeqDict[taxonOrderingID[j]]
                for k in range(len(currentSequence)):
                    if currentSequence[k] != otherSequence[k]:
                        disimilarityCounter+=1
                    totalCounter += 1
                disimilarityPercentage = disimilarityCounter/ float(totalCounter)
            
            disimilarityRow.append(disimilarityPercentage) #add dissimilarity percentages to the current row
            if(j == len(taxonOrderingID)-1):
                tableString += str(disimilarityPercentage) + "\n"
            else:
                tableString += str(disimilarityPercentage) + "\t"
        disimilarityTable.append(disimilarityRow)              #add built row to dissimilarity table/ distance matrix
    
    if outputToAFile: # For wanting to ouptu to a file
        with open(OUTPUT_FILE, "w") as file: #write string representation of the initial distance matrix to the out-file
            file.write(headerRowString)
            file.write(tableString)

    return (disimilarityTable, taxonOrderingID, sequenceOrdering) #return the newly built 2-d distance matrix, respectvie sequence ID headers, and ID to ordering

###########
#TODO: 
#BUILDS A QMATRIX based on distance matrix
###########
def makeQ(d): #d is the distance/dissimilarity 2-d matrix at the time
    Qmatrix = []
    
    #initialize a empty/0 Q-Matrix
    for i in range(len(d)):
        row = []
        for j in range(len(d)):
            row.append(0)
        Qmatrix.append(row)
        
    #build a Qmatrix
    for i in range(len(d)-1): #iterate through the top diagonal half of the the dissimilarity matrix
        for j in range(i+1, len(d[i])): #ab, ac, ad, ae -> bc, bd, be -> cd, ce -> de
            
            iSummation = 0
            for k in d[i]: #iterate through the i-row
                iSummation += k

            jSummation = 0
            for l in d[j]: #iterate through the j-row
                jSummation += l

            # Formula for Q matrix computation and setting the Qmatrix to its respective values 
            Qn_ij = ((len(d) - 2) * d[i][j]) - iSummation - jSummation
            Qmatrix[i][j] = Qn_ij
            Qmatrix[j][i] = Qn_ij
    
    return Qmatrix

###########
#TODO: 
#FIND THE INDEXES (i,j) FOR THE SMALLEST VALUE IN QMATRIX 
###########
def QMin_ij(q): #q is the Q 2-d matrix
    smallestDistance = sys.maxsize
    minIJ = (None, None)
    for i in range(len(q)-1):
        for j in range(i+1, len(q)):
            if(q[i][j]) < smallestDistance:
                smallestDistance = q[i][j]
                minIJ = (i, j)
    return minIJ

###########
#TODO
# Take in an fna file, containing multiple rows of differing rna sequences
# Output can be written to edges.txt for built tree unless specified otherwise
# Output should be written to genetic-distances.txt for the original distance matrix, or can be flagged to not do so
def geneticDistance(input_file, treeOutFile="edges.txt", outputOriginalDistanceMatrix=True):

    # Initial Distance matrix is created based on fna file 
    disimilarityTable, taxonOrderingID, IdToOriginalIdx = makeD(input_file, outputOriginalDistanceMatrix)
    # Contains a list of tips, and future internal nodes, for tree building purposes 
    tips = originalTips(taxonOrderingID)
    
    print("Here is the initial distance matrix")
    mprint(disimilarityTable)
    print("Here is the initial taxonIDs in the distance matrix")
    print(taxonOrderingID)
    print("Here is the initial tip nodes of all the sequences")
    print(tips)
    
    nextAvailableID = int(len(taxonOrderingID)) + 1 #if there are taxon sequences 1,2,3,4 initially, then the new taxon/node ID for merging is 5
    nextAvailableID2 = int(len(taxonOrderingID)) + 1 #secondary copy for later use

    while(len(disimilarityTable) > 3): #neigbor join until there are only 3 nodes left
        print("The size of the table is n*n where n=" + str(len(disimilarityTable)))
        
        #build a Qmatrix
        print("Qmatrix is: ")
        Qmatrix = makeQ(disimilarityTable)
        mprint(Qmatrix)
        smallest_IJ = QMin_ij(Qmatrix) #tuple containing the index for the smallest Q value in the matrix
        smallestI = smallest_IJ[0]
        smallestJ = smallest_IJ[1]

        #join taxons with indexes for the smallest value within Qmatrix; smallestI and smallestJ
        print("Smallest value in the Qmatrix")
        print("smallestI: " + str(smallestI))
        print("smallestJ: " + str(smallestJ))

        #perform branch length estimation for joining nodes to a new ancestor
        smallIsum = 0
        smallJsum = 0
        for k in range(len(disimilarityTable)):
            smallIsum += disimilarityTable[smallestI][k]
            smallJsum += disimilarityTable[smallestJ][k]

        length_iu = (disimilarityTable[smallestI][smallestJ] / 2.0) + (1.0 / ((2.0 * (len(disimilarityTable)-2))) * (smallIsum - smallJsum))
        length_ju = disimilarityTable[smallestI][smallestJ] - length_iu

        #create a new internal node and have its children be the correspnding smallest tips that are being joined
        newInternalNode = tn(nextAvailableID, children=[tips[smallestI], tips[smallestJ]])
        tips[smallestI].parent = newInternalNode
        tips[smallestJ].parent = newInternalNode
        
        # Set the calculated lengths
        tips[smallestI].length = length_iu
        tips[smallestJ].length = length_ju

        # remake a tips list to keep track of the tips that are leaves, remove recently joined leaves smallestI and smallestJ and add new parent leaf
        newTips = []
        for tipIdx in range(len(tips)):
            if(tipIdx == smallestI) or (tipIdx == smallestJ):
                continue
            else:
                newTips.append(tips[tipIdx])
        newTips.insert(0, newInternalNode)
        # update tips to reflect all current leaves
        tips = list(newTips)
        
        # make updates to distance matrix for the next iteration of this
        # make a newTaxonOrderingID to facilitate and keep track of the taxon/sequence ID's within the distance matrix still
        # EX. old: [A, B, C, D, E] -> [U, C, D, E]
        newTaxonOrderingID = []
        for idIndex in range(len(taxonOrderingID)):
            if(idIndex == smallestI) or (idIndex == smallestJ):
                continue
            else:
                newTaxonOrderingID.append(taxonOrderingID[idIndex])
        newTaxonOrderingID.insert(0, nextAvailableID)           #add new node to the front of the newTaxonOrderingID list
        
        print("TAXON ORDERING ID's: before and after")
        print(taxonOrderingID)
        print(newTaxonOrderingID)   

        #Create a new distance/dissimilarity table to hold new and carrying over distances, size (n-1) * (n-1)
        newDisimilarityTable = []
        for i in range(len(newTaxonOrderingID)): #initialize a 0 filled updated dissimilarity matrix
            row = []
            for j in range(len(newTaxonOrderingID)):
                row.append(0)
            newDisimilarityTable.append(row)
        
        #print("ENSURE THAT BOTH NEW TABLE AND NEW TAXON ORDERING ID HAVE SAME SIZES")
        #print("new Taxon ordering ID size: " + str(len(newTaxonOrderingID)))
        #print("new disimilarity table size: " + str(len(newDisimilarityTable)))
        
        #add the new estimated distance for the new node based on the others 
        #first row and column of the new distance/disimilarity matrix
        for i in range(1, len(newTaxonOrderingID)):
            currentDistanceAproxTaxonID = taxonOrderingID.index(newTaxonOrderingID[i])
            newDisimilarityTable[0][i] = 0.5 * (disimilarityTable[smallestI][currentDistanceAproxTaxonID] + disimilarityTable[smallestJ][currentDistanceAproxTaxonID] - disimilarityTable[smallestI][smallestJ])
            newDisimilarityTable[i][0] = newDisimilarityTable[0][i]
        
        print("NEW DISTANCE MATRIX AFTER WITH NEW DISTANCES FOR U-NODE *not complete yet*\n")
        mprint(newDisimilarityTable)

        #Incorporate old dissimilarity matrix values into the new one
        for i in range(1, len(newTaxonOrderingID)-1):
            for j in range(i+1, len(newTaxonOrderingID)):
                
                for k in range(len(taxonOrderingID)):
                    if taxonOrderingID[k] == newTaxonOrderingID[i]:
                        iReferenceIndex = k
                    if taxonOrderingID[k] == newTaxonOrderingID[j]:
                        jReferenceIndex = k
                newDisimilarityTable[i][j] = disimilarityTable[iReferenceIndex][jReferenceIndex] 
                newDisimilarityTable[j][i] = disimilarityTable[iReferenceIndex][jReferenceIndex]
        
        #update variables for next iteration
        print("NEW DISTANCE MATRIX:\n")
        mprint(newDisimilarityTable)
        disimilarityTable = []
        for i in range(len(newDisimilarityTable)):
                disimilarityTable.append(list(newDisimilarityTable[i]))
        taxonOrderingID = list(newTaxonOrderingID)
        nextAvailableID+=1
    
    #OUT OF LOOP - At this point, we have a disimilarity matrix of n*n where n = 3
    print("Last Iteration!!!! here is the final Q matrix")
    lastQ = makeQ(disimilarityTable) 
    mprint(lastQ)
    #ASSUME THE MATRIX IS N*N where n=3, a bit of hardcoding here based on wiki implementation of algorithm   

    #merge i (idx=1) and j(idx=2) to root aswell as last node
        #calculate final distances
    length_wi = (0.5 * disimilarityTable[1][2]) + ((0.5 * (1.0 / (len(disimilarityTable)-2))) * ((disimilarityTable[1][0] + disimilarityTable[1][2]) - (disimilarityTable[2][0] + disimilarityTable[2][1])))
    length_wj = disimilarityTable[1][2] - length_wi
    length_wu = disimilarityTable[1][0] - length_wi

    # Update tree                            Node u    Node i   Node j
    rootNode = tn(nextAvailableID, children=[tips[0], tips[1], tips[2]])
    tips[0].parent = rootNode
    tips[1].parent = rootNode
    tips[2].parent = rootNode
    # Set the calculated lengths
    tips[0].length = length_wu
    tips[1].length = length_wi
    tips[2].length = length_wj
    # Set last internal node to root
    tr = rootNode

    #print final lengths 
    print("Joining parent: " + str(nextAvailableID) + ", With " + str(taxonOrderingID[1]) +  " Their distance is: " + str(length_wi))
    print("Joining parent: " + str(nextAvailableID) + ", With " + str(taxonOrderingID[2]) +  " Their distance is: " + str(length_wj))
    print("Joining parent: " + str(nextAvailableID) + ", With " + str(taxonOrderingID[0]) +  " Their distance is: " + str(length_wu))

    writeTreeToFile(treeOutFile, tr, nextAvailableID2, IdToOriginalIdx)

    #Return the created tree for reusage
    return tr

if __name__ == "__main__": #Main that takes in terminal arguments and calls dissimilarity based on arguments 
    if(len(sys.argv) != 2):
        print("Expecting 1 arguments: (1) FNA file containing sequences for calculating dissimilarity\n")
    else:
        print(sys.argv[1])
        originalTree = geneticDistance(sys.argv[1])
        bootstrap(sys.argv[1], 100, originalTree)