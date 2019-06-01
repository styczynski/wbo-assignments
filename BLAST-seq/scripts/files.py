#
# Helper module for handling files
#
import os
import re
import urllib
import urllib.request
from scripts.logger import logInfo, logError, logDebug

# Check https://regex101.com/r/A326u1/5 for reference
DOMAIN_FORMAT = re.compile(
    r"(?:^(\w{1,255}):(.{1,255})@|^)" # http basic authentication [optional]
    r"(?:(?:(?=\S{0,253}(?:$|:))" # check full domain length to be less than or equal to 253 (starting after http basic auth, stopping before port)
    r"((?:[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?\.)+" # check for at least one subdomain (maximum length per subdomain: 63 characters), dashes in between allowed
    r"(?:[a-z0-9]{1,63})))" # check for top level domain, no dashes allowed
    r"|localhost)" # accept also "localhost" only
    r"(:\d{1,5})?", # port [optional]
    re.IGNORECASE
)
SCHEME_FORMAT = re.compile(
    r"^(http|hxxp|ftp|fxp)s?$", # scheme: http(s) or ftp(s)
    re.IGNORECASE
)

def validateUrl(url):
    url = url.strip()

    if not url:
        raise Exception("No URL specified")
    if len(url) > 2048:
        raise Exception("URL exceeds its maximum length of 2048 characters (given length={})".format(len(url)))

    result = urllib.parse.urlparse(url)
    scheme = result.scheme
    domain = result.netloc

    if not scheme:
        raise Exception("No URL scheme specified")
    if not re.fullmatch(SCHEME_FORMAT, scheme):
        raise Exception("URL scheme must either be http(s) or ftp(s) (given scheme={})".format(scheme))
    if not domain:
        raise Exception("No URL domain specified")
    if not re.fullmatch(DOMAIN_FORMAT, domain):
        raise Exception("URL domain malformed (domain={})".format(domain))
    return url

def ensureTmpIsPresent(tmp, logger=None):
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    else:
        logInfo(logger, "Clearing temporary download directory")
        for the_file in os.listdir(tmp):
            file_path = os.path.join(tmp, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                    logError(logger, e)

#
# Get file by path.
# Path can be relative/absolute path or URL.
#
def getFile(path, logger=None, tmp=None):
    ensureTmpIsPresent(tmp, logger=logger)
    
    logDebug(logger, "Getting file {}".format(path))
    
    isValidFile = os.path.isfile(path)
    if isValidFile:
        logDebug(logger, "Path is a valid file.")
        return path
    else:
        logDebug(logger, "Path is not a valid file path")
    
    isValidUrl = True
    try:
        validateUrl(path)
    except:
        isValidUrl = False
        
    if isValidUrl:
        logDebug(logger, "Path seems like a fine url. Try to download the file")
        downloadPath = os.path.join(tmp, "download.fasta")
        logInfo(logger, "Downloading file {} from remote source".format(path))
        urllib.request.urlretrieve(path, downloadPath)
        return downloadPath
    else:
        logError(logger, "The given path {} is not a valid url nor points to any local file.".format(path))
        exit(1)
    
    