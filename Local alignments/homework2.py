import numpy as np
import numba as nb
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo
from Bio.Alphabet import generic_dna
from Bio import SeqUtils

def initialStep(V0, V1, InSeq, In, M, Dir, isLocal = False):
    d = 8
    for i in range(V0.shape[0]):
        top = np.iinfo(np.int16).min
        if i > 0:
            top = V1[i-1] - d
        else:
            top = 0
        V1[i] = top
        if isLocal and top < 0:
            V1[i] = 0
            Dir[i] = 3
        else:
            Dir[i] = 2

def nextStep(V0, V1, InSeq, In, M, Dir, isLocal = False):
    d = 8
    for i in range(V0.shape[0]):
        left = V0[i] - d
        top = np.iinfo(np.int16).min
        corner = np.iinfo(np.int16).min
        if i > 0:
            top = V1[i-1] - d
            if (In,InSeq[i-1]) in M:
                corner = V0[i-1] + M[In,InSeq[i-1]]
            elif (InSeq[i-1],In) in M:
                corner = V0[i-1] + M[InSeq[i-1], In]
        V1[i] = max(left, top, corner)
        if isLocal:
            V1[i] = max(V1[i], 0)
        if V1[i] == top:
            Dir[i] = 2
        elif V1[i] == corner:
            Dir[i] = 1
        elif V1[i] == left:
            Dir[i] = 0
        else:
            Dir[i] = 3
    
def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))
    
def alignment_traceback(codonsA, codonsB, alignment, dnaARange = (0, -1), dnaBRange = (0, -1), ):
    (scoreMatrix, dirMatrix) = alignment
    
    colsCount = dirMatrix.shape[0]
    rowsCount = dirMatrix.shape[1]
    colIndex = dnaARange[1] + 1
    rowIndex = dnaBRange[1] + 1
    colIndexEnd = dnaARange[0]
    rowIndexEnd = dnaBRange[0]
    
    shift = 1
    
    if dnaARange[1] < 0:
        colIndex = colsCount - 1
        colIndexEnd = 0
        shift=2
    
    if dnaBRange[1] < 0:
        #rowIndex = np.argmax(scoreMatrix[colIndex])
        rowIndex = rowsCount - 1
        rowIndexEnd = 0
        shift=2
    
    score = scoreMatrix[colIndex, rowIndex];
    
    accA = np.chararray(colsCount + rowsCount)
    accB = np.chararray(colsCount + rowsCount)
    accAIter = colsCount + rowsCount - 1
    accBIter = colsCount + rowsCount - 1
    
    while True:
        if (colIndex < colIndexEnd or rowIndex < rowIndexEnd):
            break
    
        if dirMatrix[colIndex, rowIndex] == 0:
            colIndex = colIndex - 1
            accA[accAIter] = codonsA[colIndex]
            accAIter = accAIter - 1
            accB[accBIter] = '-'
            accBIter = accBIter - 1
        elif dirMatrix[colIndex, rowIndex] == 1:
            colIndex = colIndex - 1
            rowIndex = rowIndex - 1
            accA[accAIter] = codonsA[colIndex]
            accAIter = accAIter - 1
            accB[accBIter] = codonsB[rowIndex]
            accBIter = accBIter - 1
        elif dirMatrix[colIndex, rowIndex] == 2:
            rowIndex = rowIndex - 1
            accA[accAIter] = '-'
            accAIter = accAIter - 1
            accB[accBIter] = codonsB[rowIndex]
            accBIter = accBIter - 1
        else:
            break
    
    strA = np.chararray.tostring(accA[accAIter+shift:])
    strB = np.chararray.tostring(accB[accBIter+shift:])
    return (strA, strB, score, 0, len(strA))
    
def codons_alignment(codonsA, codonsB, isLocal = False):
    lenA = len(codonsA)
    lenB = len(codonsB)
    
    if lenA < lenB:
        lenA, codonsA, lenB, codonsB = lenB, codonsB, lenA, codonsA
    
    currentState = np.zeros(lenB+1, dtype=np.int16)
    nextState = np.zeros(lenB+1, dtype=np.int16)
    dirMatrix = np.zeros((lenA+1, lenB+1), dtype=np.int16)
    scoreMatrix = np.zeros((lenA+1, lenB+1), dtype=np.int16)
    
    initialStep(currentState, nextState, codonsB, codonsA[0], MatrixInfo.blosum50, dirMatrix[0], isLocal=isLocal)
    currentState, nextState = nextState, currentState
    
    for i in range(lenA):
        scoreMatrix[i] = currentState
        nextStep(currentState, nextState, codonsB, codonsA[i], MatrixInfo.blosum50, dirMatrix[i+1], isLocal=isLocal)
        currentState, nextState = nextState, currentState
    scoreMatrix[lenA] = currentState

    #print(scoreMatrix)
    return (scoreMatrix, dirMatrix)
    
def codons_align(codonsA, codonsB, isLocal = False, dnaARange = (0, -1), dnaBRange = (0, -1)):
    return (alignment_traceback(codonsA, codonsB, codons_alignment(codonsA, codonsB, isLocal=isLocal), dnaBRange=dnaBRange, dnaARange=dnaARange))
    
def dna_local_align3(dnaA, dnaB):

    results = []
    for i in range(3):
        for j in range(3):
            
            shiftA = i
            shiftB = j
            codonsA = Seq.translate(pad_seq(dnaA[shiftA:]))
            codonsB = Seq.translate(pad_seq(dnaB[shiftB:]))
            
            print(codonsA)
            print(codonsB)
            k = codons_align(codonsA, codonsB)
            
            if k[4] > 0:
                results.append(k)
        
    return sorted(results, key=lambda r: -r[2])
    
dnaA = Seq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", generic_dna)
dnaB = Seq("GTGGCCATTGTAATGGAAAGGGTGAAAGAT", generic_dna)

print(dna_local_align3(dnaA, dnaB))