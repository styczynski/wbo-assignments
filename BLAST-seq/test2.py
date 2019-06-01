from Bio.Blast.NCBIWWW import qblast

query = (
    ">s1\n"
    "MHEIKYITIDEADVLLTEEHEETTRFICQSANRDRQISLFSATTSERLDNFFDKVESSQQ\n"
    "IEVVAGEAKMPTTIDHIYIQVNPRDKVKTLYRLAQVENMRAIVFVNTIGRLNTVYEKLNH\n"
    "DGVKISALHGDLSKLQRQESVRDFKKGETSLLLATDVAARGIDLPNLPAIIQFDMAQSLT\n"
    "QYVHRSGRTGRMGEQGAAISLVTDREARELKQMVKENDVKMIEQIVKFGHLIDPQKTK"
)

r = qblast("blastp", "nr_v5", query).getvalue()
print(r)