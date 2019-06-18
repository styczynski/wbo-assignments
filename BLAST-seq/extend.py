import click
from scripts.logger import createLogger, setLoggerLevel, logInfo, logDebug
from scripts.files import getFile
from Bio.Blast.NCBIWWW import qblast
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import asyncio
import urllib
import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor
import xml.etree.ElementTree as ET
import time
from scripts.logger import logInfo, logError, logDebug

DEFAULT_FILE_TMP_FOLDER = "./workdir"
DEFAULT_ID_THRESHOLD = 90
DEFAULT_E_VALUE = "10e-10"

#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=1168090635&rettype=fasta&retmode=text
def entrezRetrieveSequence(accessionIDs, tmp='.', entrezUrl='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', dbName='protein', retries=5):
    
    if len(accessionIDs) <= 0:
        return []
    
    downloadedFilePath = os.path.join(tmp, "etrez.{}.{}.fasta".format(accessionIDs[0] + "_batch__" + str(len(accessionIDs)), dbName))
    
    retryNo = 0
    while retryNo < retries:
        retryNo = retryNo + 1
        try:
            urllib.request.urlretrieve('{}?db={}&id={}&rettype=fasta&retmode=text'.format(entrezUrl, dbName, ",".join(accessionIDs)), downloadedFilePath)
            break
        except urllib.error.HTTPError as err:
            if retryNo < retries-1:
                logError(logger, "Request failed retrying in 5 seconds...")
                time.sleep(5)
            else:
                raise err
    
    ret = []
    with open(downloadedFilePath, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ret.append(record)
        
    return ret

# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]
    
#
# MAIN CODE BLOCK
#
@click.command(short_help="Search for extensions of given FASTA sequences with BLAST web")
@click.argument('input', nargs=-1)
@click.option('--loggingLevel', '-l', default='INFO', type=click.Choice(['INFO', 'DEBUG']), help='Set logging level')
@click.option('--tmp', '-t', default=DEFAULT_FILE_TMP_FOLDER, help="Path to save temporary data. Default is {}".format(DEFAULT_FILE_TMP_FOLDER))
@click.option('--idThreshold', '-i', default=DEFAULT_ID_THRESHOLD, help="Set minimal identity threshold in percents. Default is {}".format(DEFAULT_ID_THRESHOLD))
@click.option('--eValue', '-e', default=DEFAULT_E_VALUE, help="Set minimal E-value using scientific notation. Default is {}".format(DEFAULT_E_VALUE))
def main(logginglevel, input, tmp, idthreshold, evalue):
    logger = createLogger(__file__)
    logger = setLoggerLevel(logger, logginglevel)
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        loop = asyncio.get_event_loop()
        
        inputFiles = [ getFile(inputFile, logger=logger, tmp=tmp) for inputFile in input ]
        outputFiles = []
        for file in inputFiles:

            fileOut = file.replace(".fasta", "_ext.fasta")
            logInfo(logger, "Extending file {} to {}".format(file, fileOut))
            seqs = []
            
            with open(file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seqs.append("> {}\n{}".format(record.id, str(record.seq)))
            queryStr = "\n".join(seqs)
            
            logInfo(logger, "Requesting web BLASTP search")
            requestOut = qblast("blastp", "nr_v5", queryStr, expect=float(evalue), perc_ident=float(idthreshold)).getvalue()
            
            root = ET.fromstring(requestOut)
            
            accessionIDs = []
            
            hits = root.findall("./BlastOutput_iterations/Iteration/Iteration_hits/Hit")
            logInfo(logger, "Got {} hits for given query".format(len(hits)))
            recordsOut = []
            recordCountOut = 0
            for hit in hits:
                seqs = hit.findall("./Hit_hsps/Hsp/Hsp_qseq")
                for seq in seqs:
                    accessionID = hit.find("./Hit_accession").text
                    accessionIDs.append(accessionID)
            
            gatherTasks = []
            
            requestsLimit = 3
            chunkSize = int(len(accessionIDs) / requestsLimit)+1
            accessionChunks = [ chunk for chunk in chunks(accessionIDs, chunkSize) ]
            logDebug(logger, "Chunked request, will request Entrez for {} batches of size {}.".format(len(accessionChunks), chunkSize))
            
            for accessionIDs in accessionChunks:
                gatherTasks.append(
                    loop.run_in_executor(
                        executor,
                        entrezRetrieveSequence,
                        *((accessionIDs, tmp))
                    )
                )
            
            gatherFuture = asyncio.ensure_future(asyncio.gather(*gatherTasks))
            loop.run_until_complete(gatherFuture)
            
            recordsOut = []
            for records in gatherFuture.result():
                for record in records:
                    recordsOut.append(record.format("fasta"))
            
            recordCountOut = len(recordsOut)
            logDebug(logger, "Writing to file {}".format(fileOut))
            with open(fileOut, "w") as fileOutHandle:
                fileOutHandle.write("\n".join(recordsOut))
                
            with open(fileOut, "r") as handle:
                recordCount = 0
                for record in SeqIO.parse(handle, "fasta"):
                    recordCount = recordCount + 1
                if recordCount != recordCountOut:
                    logError(logger, "Got mismatch records count after writing output file. {} records are present in {} and there should be {} reconds.".format(recordCount, fileOut, recordCountOut))
            outputFiles.append(fileOut)
            
        logInfo(logger, "Written {} files in total. Returning them to stdout as plaintext paths list".format(len(outputFiles)))
        print("\n".join(outputFiles))
    
# Entry
if __name__ == "__main__":
    main()