# Needleman–Wunsch Alignment Calculator

Using the Needleman–Wunsch algorithm, This module can generate an alignment matrix, get the alignment score, and generate an alignment for DNA sequences.

## Class Methods

### __init__(self, matchScore, mismatchScore, gapPenalty)
The **AlignmentCalculator** class initializes a calculator with the following parameters:
* _matchScore_ - score that is added to the cell if there is an alignment
* _mismatchScore_ - score that is added to the cell if there is a misalignment
* _gapPenalty_ - score that is added to the cell if there is an gap

### generateAlignmentMatrix(self, seq1, seq2)
Generates and returns an alignment matrix made from the following parameters:
* _seq1_ - sequence to be compared with _seq2_
* _seq2_ - sequence to be compared with _seq1_

### generateAlignments(self, seq1, seq2, matrix):
Returns a tuple containing two strings, the alignment for seq1 and for seq2 from the following parameters:
_seq1_ - sequence to be aligned with _seq2_
_seq2_ - sequence to be aligned with _seq1_
_matrix_ - alignment matrix generated using the _generateAlignmentMatrix_ method

### @classmethod getScore(cls, matrix)
Returns the score of a matrix from the following parameters:
_matrix_ - alignment matrix generated using the _generateAlignmentMatrix_ method

## Test File
This repository includes a file names _test.csv_ which can be used to test the functionality of the code.
Run **python AlignmentCalculator.py test.csv** in order to test the class methods using the test file.