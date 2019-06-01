import click
from scripts.logger import createLogger, setLoggerLevel, logInfo, logDebug
from scripts.files import getFile
from Bio.Blast.NCBIWWW import qblast
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import xml.etree.ElementTree as ET

DEFAULT_FILE_TMP_FOLDER = "./workdir"
DEFAULT_ID_THRESHOLD = 90
DEFAULT_E_VALUE = "10e-10"

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
        
        hits = root.findall("./BlastOutput_iterations/Iteration/Iteration_hits/Hit")
        logInfo(logger, "Got {} hits for given query".format(len(hits)))
        recordsOut = []
        recordCountOut = 0
        for hit in hits:
            seqs = hit.findall("./Hit_hsps/Hsp/Hsp_qseq")
            seqNo = 0
            for seq in seqs:
                seqID = hit.find("./Hit_accession").text
                seqNo = seqNo + 1
                if len(seqs) > 1:
                    seqID = "{}_{}".format(seqID, seqNo)
                record = SeqRecord(Seq(seq.text, IUPAC.protein),
                   id=seqID, name=hit.find("./Hit_id").text,
                   description=hit.find("./Hit_def").text)
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