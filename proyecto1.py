import sys

# prints all the rows of a 2d array
def printMatrix(matrix):
    for row in matrix: print(row)  

# Class for a calculator that generates an alignment matrix for the
# input sequences and generates an alignment for them
class AlignmentCalculator():
    # initializes all variables for the calculator and
    # generates a matrix with only the initial values.
    def __init__(self, seq1, seq2, matchScore=1, mismatchScore=-1, gapPenalty=-2):
        # setting up variables as well as width and height of alignment matrix
        self.gapPenalty = gapPenalty
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore
        self.gapPenalty = gapPenalty
        self.width = len(seq1)+1
        self.height = len(seq2)+1
        self.seq1 = seq1
        self.seq2 = seq2

        # initializing alignment matrix
        self.alignmentMatrix = [[0]*self.width for _ in range(self.height)]

        # setting initial values of matrix
        for y in range(self.height): self.alignmentMatrix[y][0] = gapPenalty * y
        for x in range(self.width): self.alignmentMatrix[0][x] = gapPenalty * x

        # automatically filling matrix upon initializaion
        self.fillMatrix()

    # populates the alignment matrix
    def fillMatrix(self):
        for y in range(1,self.height):
            for x in range(1,self.width):
                # setting the score value to 1 if seq1[current] = seq2[current]
                score = self.matchScore if self.seq1[x-1] == self.seq2[y-1] else self.mismatchScore

                topLeftScore = self.alignmentMatrix[y-1][x-1]
                leftScore = self.alignmentMatrix[y][x-1]
                topScore = self.alignmentMatrix[y-1][x]

                # deciding which value to put into the current cell
                value = max(topScore + self.gapPenalty, leftScore + self.gapPenalty, topLeftScore + score)
                self.alignmentMatrix[y][x] = value
    
    # returns a tuple containing two strings. seq1 aligned and seq2 aligned.
    def getAlignments(self):
        alignedSeq1 = []
        alignedSeq2 = []

        y = self.height - 1
        x = self.width - 1

        # iterating backwards through the matrix matrix
        while x > 0 or y > 0:
            currentScore = self.alignmentMatrix[y][x]

            # trying to see if the calculation using the left cell resulted in the max value
            if self.alignmentMatrix[y][x-1] + self.gapPenalty == currentScore:
                alignedSeq1.append(self.seq1[x-1])
                alignedSeq2.append('-')
                x -= 1
            # trying to see if the calculation using the top cell resulted in the max value
            elif self.alignmentMatrix[y-1][x] + self.gapPenalty == currentScore:
                alignedSeq2.append(self.seq2[y-1])
                alignedSeq1.append('-')
                y -= 1
            # when there is no gap, just append the next letter in the sequence
            else:
                alignedSeq1.append(self.seq1[x-1])
                alignedSeq2.append(self.seq2[y-1])
                x -= 1
                y -= 1
        
        # reversing the list to get the proper alignment and joining all the values into a string
        # it is done this way because it is more efficient than appending constantly to a string.
        # this is because strings are immutable and every time you "append" to a string, you're
        # actually copying all the characters over to a new string
        return (''.join(alignedSeq1[::-1]), ''.join(alignedSeq2[::-1]))
        
    # returns the score of the matrix
    def getScore(self): 
        return self.alignmentMatrix[-1][-1]

# returns a list of all the sequence pairs in tuples
def getSequencePairs():
    if len(sys.argv) < 2: return ()
    path = sys.argv[1]

    with open(path) as csv:
        # using line breaks as separators in order to get a list of all the rows in the csv.
        # the list is spliced from index 1 onwards in order to remove the csv titles/headers
        pairs = csv.read().split('\n')[1:]
        
        # creating a tuple for each row which contains both
        # sequences by using the comma as the separator.
        return (pair.split(',') for pair in pairs) 

def main():
    # getting list of sequence pairs from csv file
    sequencePairs = getSequencePairs()
    
    # computing and printing the alignment values for all the sequences
    for sequences in sequencePairs:
        if len(sequences) < 2: continue
        seq1, seq2 = sequences

        calculator = AlignmentCalculator(seq1, seq2)

        alignedSeq1, alignedSeq2 = calculator.getAlignments()
        print(alignedSeq1, alignedSeq2, calculator.getScore())

if __name__ == "__main__":
    main()