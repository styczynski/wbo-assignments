import requests
import urllib.parse


def pretty_print_POST(req):
    """
    At this point it is completely built and ready
    to be fired; it is "prepared".

    However pay attention at the formatting used in 
    this function because it is programmed to be pretty 
    printed and may differ from the actual request.
    """
    print('{}\n{}\n{}\n\n{}'.format(
        '-----------START-----------',
        req.method + ' ' + req.url,
        '\n'.join('{}: {}'.format(k, v) for k, v in req.headers.items()),
        req.body,
    ))


def makeWebBlastRequest(program, database, queries, blastUrl="https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", queryParams={}, queryHeaders={}):
    defaultParams = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": "\n".join(queries), #[ urllib.parse.quote(query.encode("utf-8")) for query in queries ]
    }
    defaultHeaders = {
        "content-type": "application/x-www-form-urlencoded",
        "user-agent": "Mozilla/5.0"
    }
    #data = "&".join([ "{}={}".format(key, value) for key, value in ({ **defaultParams, **queryParams }).items() ])
    
    r = requests.Request('POST', blastUrl, data={ **defaultParams, **queryParams, "CMD": "request" }, headers={ **defaultHeaders, **queryHeaders }).prepare()
    pretty_print_POST(r)
    s = requests.Session()
    s.send(r)
    
    r = requests.Request('POST', blastUrl, data={ **defaultParams, **queryParams }, headers={ **defaultHeaders, **queryHeaders }).prepare()
    pretty_print_POST(r)
    return s.send(r)
    
    
r = makeWebBlastRequest("blastp", "nr_v5", [(
    ">s1"
    "MHEIKYITIDEADVLLTEEHEETTRFICQSANRDRQISLFSATTSERLDNFFDKVESSQQ"
    "IEVVAGEAKMPTTIDHIYIQVNPRDKVKTLYRLAQVENMRAIVFVNTIGRLNTVYEKLNH"
    "DGVKISALHGDLSKLQRQESVRDFKKGETSLLLATDVAARGIDLPNLPAIIQFDMAQSLT"
    "QYVHRSGRTGRMGEQGAAISLVTDREARELKQMVKENDVKMIEQIVKFGHLIDPQKTK"
)])

print(r.status_code)
print(r.text)