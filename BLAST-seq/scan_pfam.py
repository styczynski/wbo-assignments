import click
import csv
import re
from scripts.logger import createLogger, setLoggerLevel, logInfo, logError, logDebug
from scripts.files import getFile
from Bio import SeqIO
import xml.etree.ElementTree as ET
import xml
import requests
import asyncio
from concurrent.futures import ThreadPoolExecutor
import urllib.parse

DEFAULT_FILE_TMP_FOLDER = "./workdir"

def makeHMMScanRequest(seq, seqID=None, logger=None, hmmUrl="https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan", queryParams={}, queryHeaders={}):
    if seqID:
        logDebug(logger, "Doing HMMScan on seqence \"{}\"...".format(seqID))
    
    defaultParams = {
        "hmmdb": (None, "pfam"),
        "seq": (None, seq),
        "threshold": (None, "cut_ga")
    }
    defaultHeaders = {
        "Accept": "text/xml",
        "User-Agent": "Mozilla/5.0",
    }
    
    s = requests.Session()
    r = requests.Request('POST', hmmUrl, files={ **defaultParams, **queryParams }, headers={ **defaultHeaders, **queryHeaders }).prepare()
    answer = s.send(r)
    answerContent = answer.content.decode("utf-8")
    
    if seqID:
        logDebug(logger, "Done HMMScan on seqence \"{}\"".format(seqID))
    
    # Replace invalid tags
    if answerContent:
        answerContent = re.sub(r'<\d+ H=.*\/>', '', answerContent)
    
    xmlTree = None 
    try:
        xmlTree = ET.fromstring(answerContent)
    except:
        logError(logger, "Failed HMMScan on sequence \"{}\" - Invalid response ".format(seqID))
        logError(logger, "xml [{}]".format(answerContent))
        raise Exception("XML Parse failed")
        
    if seqID:
        return {
            "seqID": seqID,
            "result": xmlTree
        }
    return xmlTree

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
    
    inputFiles = [ getFile(inputFile, logger=logger, tmp=tmp) for inputFile in input ]
    outputFiles = []
    
    for file in inputFiles:
        with ThreadPoolExecutor(max_workers=10) as executor:
            loop = asyncio.get_event_loop()
            seqTasks = []
        
            fileOut = file.replace(".fasta", "_hmm.csv")
            logInfo(logger, "Searching with sequences from {} saving to {}".format(file, fileOut))
            
            domainToSeqMap = {}
            seqToDomainMap = {}
            
            with open(file, "r") as handle:
            
                recordMap = {}
                for record in SeqIO.parse(handle, "fasta"):
                    recordMap[record.id] = record
                    if not (record.id in seqToDomainMap.keys()):
                        seqToDomainMap[record.id] = []
                    seqTasks.append(
                        loop.run_in_executor(
                            executor,
                            makeHMMScanRequest,
                            *("> {}\n{}".format(record.id, str(record.seq)), record.id, logger)
                        )
                    )
                    
                gatherFuture = asyncio.ensure_future(asyncio.gather(*seqTasks))
                loop.run_until_complete(gatherFuture)
                
                for answer in gatherFuture.result():
                    record = recordMap[answer["seqID"]]
                    root = answer["result"]
                    for hit in root.findall("./data/hits"):
                        nincluded = int(hit.attrib["nincluded"])
                        name = hit.attrib["name"]
                        if nincluded > 0:
                            seqToDomainMap[record.id].append(name)
                            if not (name in domainToSeqMap.keys()):
                                domainToSeqMap[name] = []
                            domainToSeqMap[name].append(record.id)
                
            logInfo(logger, "Searched in total {} sequences with {} returned unique domains".format(len(seqToDomainMap.keys()), len(domainToSeqMap.keys())))
            logInfo(logger, "Writing output to {}".format(file, fileOut))
            with open(fileOut, mode="w") as fileOutHandle:
                fileOutWriter = csv.writer(fileOutHandle, delimiter=",", quotechar="\"", quoting=csv.QUOTE_MINIMAL)
                headers = [ "seq_id" ] + [ key for key in domainToSeqMap.keys() ]
                fileOutWriter.writerow(headers)
                for seqID, seqDomains in seqToDomainMap.items():
                    rowResults = [ seqID ]
                    for domain in headers:
                        isInDomains = 0
                        if domain in seqDomains:
                           isInDomains = 1
                        rowResults.append(isInDomains)
                    fileOutWriter.writerow(rowResults)
            
            outputFiles.append(fileOut)
        
    logInfo(logger, "Written {} files in total. Returning them to stdout as plaintext paths list".format(len(outputFiles)))
    print("\n".join(outputFiles))
    
# Entry
if __name__ == "__main__":
    main()