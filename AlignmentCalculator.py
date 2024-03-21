import sys

# prints all the rows of a 2d array (used for debugging)
# def printMatrix(matrix):
#     for row in matrix: print(row)  

# Class for a calculator that generates an alignment matrix for the
# input sequences and generates an alignment for them
class AlignmentCalculator():
    # initializes all variables for the calculator
    def __init__(self, matchScore, mismatchScore, gapPenalty):
        self.mismatchScore = mismatchScore
        self.gapPenalty = gapPenalty
        self.matchScore = matchScore

    # generates and returns an alignment matrix made from the input sequences
    def generateAlignmentMatrix(self, seq1, seq2):
        width = len(seq1)+1
        height = len(seq2)+1

        # creating matrix with initial values
        alignmentMatrix = [[0]*width for _ in range(height)]

        for y in range(height): alignmentMatrix[y][0] = self.gapPenalty * y
        for x in range(width): alignmentMatrix[0][x] = self.gapPenalty * x

        # calculating the score for each cell in the matrix
        for y in range(1,height):
            for x in range(1,width):
                # setting the score value to either the match score or mismatch score
                score = self.matchScore if seq1[x-1] == seq2[y-1] else self.mismatchScore

                topLeftScore = alignmentMatrix[y-1][x-1]
                leftScore = alignmentMatrix[y][x-1]
                topScore = alignmentMatrix[y-1][x]

                # putting the max calculated value into the current cell
                value = max(topScore + self.gapPenalty, leftScore + self.gapPenalty, topLeftScore + score)
                alignmentMatrix[y][x] = value

        return alignmentMatrix
    
    # returns a tuple containing two strings. the alignment for seq1 and for seq2
    def generateAlignments(self, seq1, seq2, matrix):
        alignedSeq1 = []
        alignedSeq2 = []

        y = len(seq2)
        x = len(seq1)

        # iterating backwards through the matrix matrix
        while x > 0 or y > 0:
            currentScore = matrix[y][x]

            # trying to see if the calculation using the left cell resulted in the max value
            if matrix[y][x-1] + self.gapPenalty == currentScore:
                alignedSeq1.append(seq1[x-1]) # seq1 appends the letter
                alignedSeq2.append('-') # meanwhile for seq2, we add a gap
                x -= 1
            # trying to see if the calculation using the top cell resulted in the max value
            elif matrix[y-1][x] + self.gapPenalty == currentScore:
                alignedSeq2.append(seq2[y-1]) # seq2 appends the letter
                alignedSeq1.append('-') # meanwhile for seq1, we add a gap
                y -= 1
            # when there is no gap, just append the next letter in the sequence
            else:
                alignedSeq1.append(seq1[x-1])
                alignedSeq2.append(seq2[y-1])
                x -= 1
                y -= 1
        
        # reversing the list to get the proper alignment and joining all the values into a string
        # it is done this way because it is more efficient than appending constantly to a string.
        # this is because strings are immutable and every time you "append" to a string, you're
        # actually copying all the characters over to a new string
        return (''.join(alignedSeq1[::-1]), ''.join(alignedSeq2[::-1]))
        
    # returns the score of the input matrix
    @classmethod
    def getScore(cls, matrix): 
        return matrix[-1][-1]

# returns a list of all the sequence pairs in tuples
def getSequencePairs():
    if len(sys.argv) < 2: return ()
    path = sys.argv[1]

    with open(path) as csv:
        # using line breaks as separators in order to get a list of all the rows in the csv.
        # the list is spliced from index 1 onwards in order to remove the csv titles/headers.
        # rstrip is used to prevent errors caused by trailing newline characters
        pairs = csv.read().rstrip().split('\n')[1:]
        
        # creating a tuple for each row which contains both
        # sequences by using the comma as the separator.
        return (pair.split(',') for pair in pairs) 

def main():
    MATCH_SCORE=1
    MISMATCH_SCORE=-1
    GAP_PENALTY=-2

    # getting list of sequence pairs from csv file
    sequencePairs = getSequencePairs()
    
    # computing and printing the alignment values for all the sequences
    for seq1, seq2 in sequencePairs:
        # initializing calculator with necessary parameters
        calculator = AlignmentCalculator(MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY)

        alignmentMatrix = calculator.generateAlignmentMatrix(seq1, seq2)
        alignedSeq1, alignedSeq2 = calculator.generateAlignments(seq1, seq2, alignmentMatrix)

        print(alignedSeq1, alignedSeq2, AlignmentCalculator.getScore(alignmentMatrix))

if __name__ == "__main__":
    main()