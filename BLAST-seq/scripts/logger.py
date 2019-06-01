import logging

def logInfo(logger, text):
    if logger:
        logger.info("[INF]   {}".format(text))

def logDebug(logger, text):
    if logger:
        logger.debug("[DBG]   {}".format(text))
        
def logError(logger, text):
    if logger:
        logger.error("[ERR]   {}".format(text))
        
def setLoggerLevel(logger, level):
    if level == 'INFO':
        logger.setLevel(logging.DEBUG)
    elif level == 'DEBUG':
        logger.setLevel(logging.DEBUG)
    else:
        raise Exception('Invalid logging level was specified: {}'.format(level))
    return logger

def createLogger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(name)s] %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger