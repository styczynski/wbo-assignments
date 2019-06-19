import click
from scripts.logger import createLogger, setLoggerLevel, logInfo, logDebug, logError
from scripts.files import getFile
import asyncio
import urllib
import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor
import time
import csv
from numpy import genfromtxt
import scipy.special
import scipy.stats
import matplotlib.pyplot
import plotly as py
import plotly.graph_objs as go
from scripts.logger import logInfo, logError, logDebug

DEFAULT_FILE_TMP_FOLDER = "./workdir"

def loadMatrixFromCSV(filename, logger=None):
    matrix = []
    domainNames = []
    
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=",", quotechar="\"", quoting=csv.QUOTE_MINIMAL)
        rows = [ row for row in reader ]
        domainNames = rows[0][1:]
        domainCounts = []
        rows.pop(0)
        nameIndex = 0
        for name in domainNames:
            domainCounts.append(sum([ int(row[nameIndex+1]) for row in rows ]))
            nameIndex = nameIndex + 1
        domainCounts.append(len(rows))
        domainNames.append('')
        matrix = domainCounts
    
    return (domainNames, matrix)
    
def normalizeTables(tables, logger=None):
    
    globNames = [ i for i in set([ domainName for (domainNames, _) in tables for domainName in domainNames ]) ]
    newTables = []
    for (domainNames, matrix) in tables:
        newMatrix = []
        for name in globNames:
            if name in domainNames:
                index = domainNames.index(name)
                newMatrix.append(matrix[index])
            else:
                newMatrix.append(0)
        newTables.append(newMatrix)  
    
    ret = (globNames, newTables)
    return ret
    
def customFisher(domainNames, tableA, tableB, logger=None):
    
    totalA = 0
    totalB = 0
    domainIndex = 0
    for domain in domainNames:
        if domain == '':
            totalA = tableA[domainIndex]
            totalB = tableB[domainIndex]
        domainIndex = domainIndex + 1
    
    domainIndex = 0
    domainPs = []
    for domain in domainNames:
        if domain != '':
            N = totalA + totalB
            K = totalA
            n = tableA[domainIndex] + tableB[domainIndex]
            x = tableA[domainIndex]
            
            p = 0
            for i in range(x, N+1):
                p = p + ( scipy.special.binom(K, i) * scipy.special.binom((N - K), (n - i)) ) / scipy.special.binom(N, n)
            domainPs.append(p)
        
        domainIndex = domainIndex + 1
    return domainPs
    
def scipyFisher(domainNames, tableA, tableB, logger=None):
    
    totalA = 0
    totalB = 0
    domainIndex = 0
    for domain in domainNames:
        if domain == '':
            totalA = tableA[domainIndex]
            totalB = tableB[domainIndex]
        domainIndex = domainIndex + 1
    
    domainIndex = 0
    domainPs = []
    for domain in domainNames:
       
        if domain != '':
            presentA = tableA[domainIndex]
            presentB = tableB[domainIndex]
            
            logInfo(logger, "\n\nDomain {}\n".format(domain))
            logInfo(logger, "  pA={}  pB={}  tA={}  tB={}\n\n".format(presentA, presentB, totalA, totalB))
            
            
            odds_ratio, p = scipy.stats.fisher_exact([[presentA, presentB], [totalA - presentA, totalB - presentB]], alternative='greater')
            domainPs.append(p)
        
        domainIndex = domainIndex + 1
        
    return domainPs
#
# MAIN CODE BLOCK
#
@click.command(short_help="Search for extensions of given FASTA sequences with BLAST web")
@click.argument('input', nargs=-1)
@click.option('--loggingLevel', '-l', default='INFO', type=click.Choice(['INFO', 'DEBUG']), help='Set logging level')
@click.option('--tmp', '-t', default=DEFAULT_FILE_TMP_FOLDER, help="Path to save temporary data. Default is {}".format(DEFAULT_FILE_TMP_FOLDER))
def main(logginglevel, input, tmp):
    logger = createLogger(__file__)
    logger = setLoggerLevel(logger, logginglevel)
    
    inputCount = len(input)
    if inputCount != 2:
        logError(logger, "Invalid files number was specified. Please provide exactly two input files.")
        exit(1)
    
    inputA = input[0]
    inputB = input[1]
    logInfo(logger, "Fisher test of {} compared to {}".format(inputA, inputB))
    
    tableB = loadMatrixFromCSV(inputA, logger=logger)
    tableA = loadMatrixFromCSV(inputB, logger=logger)   
    (domainNames, tables) = normalizeTables([ tableA, tableB ], logger=logger)
    
    domainPsCustom = [round(v, 2) for v in customFisher(domainNames, tables[0], tables[1], logger=logger)]
    domainPsFisher = [round(v, 2) for v in scipyFisher(domainNames, tables[0], tables[1], logger=logger)]
    
    results = [ [r[0], r[1], r[2]] for r in zip(domainNames, domainPsCustom, domainPsFisher) ]
    results.sort(key=lambda x: x[1])
    
    resultsT = [[row[i] for row in results] for i in range(3)]
    table = go.Table(
        header=dict(values=['Domain', 'Custom Fisher', 'SciPy Fisher']),
        cells=dict(values=resultsT)
    )
    
    bars = go.Bar(
            x=['giraffes', 'orangutans', 'monkeys'],
            y=[20, 14, 23]
    )

    data = [table, bars] 
    py.offline.plot({"data": data}, auto_open=False)
    
# Entry
if __name__ == "__main__":
    main()