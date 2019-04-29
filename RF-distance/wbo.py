#
# [Script dependencies]
#    
#    Please install the following packages via pip:
#    $ pip install click urllib gunzip biopython networkx numpy dendropy
#
#
# [About]
#
# Python script to calculate RF metric between ClustalW and Muscle trees
#
# Usage: wbo.py [OPTIONS] URL
#
# Options:
#   -t, --tmp TEXT              Path to save temporary data. Default is
#                               ./downloads
#   -d, --disablecheck BOOLEAN  Disable DendroPy verification check
#   --help                      Show this message and exit.
#
# URL is a valid URL to download .gz archive with fasta files.
# Please use the following command to download example file:
#
#    $ python wbo.py "http://regulomics.mimuw.edu.pl/wp/wp-content/uploads/2018/05/homeo2.fa_.gz"
#
# The script calculates the RF metric in O(n) time using hashing.
# It also calculates the same metric using DendroPy to validate the result.
#
# [License]
#
# This code is licensed under MIT License.
# @Piotr Styczynski
#
import click
import urllib
from sh import gunzip
import os, shutil
from Bio import Phylo
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import random
import networkx as nx
import numpy as np
import itertools as it
from dendropy import Tree, TaxonNamespace
from dendropy.calculate import treecompare

FILE_TMP_FOLDER = "./downloads"
FILE_TO_DOWNLOAD = ""

#
# Downloads gz file referenced by url, then unpacks it to the temporary directory
# and returns list of relative paths pointing to all fasta (each file must end with .fa extension) files
# that were contained within the archive file.
#
def downloadFiles(url):
    
    if not os.path.exists(FILE_TMP_FOLDER):
        os.makedirs(FILE_TMP_FOLDER)
    else:
        print("[INFO] Clearing temporary download directory")
        for the_file in os.listdir(FILE_TMP_FOLDER):
            file_path = os.path.join(FILE_TMP_FOLDER, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)
    
    fileList = []
    filename = FILE_TMP_FOLDER+"/downloaded.fa.gz"
    
    print("[INFO] Downloading file {} from remote source".format(url))
    urllib.urlretrieve(FILE_TO_DOWNLOAD, filename)
    
    print("[INFO] Unpacking gz archive")
    gunzip(filename)
    
    for file in os.listdir(FILE_TMP_FOLDER):
        if file.endswith(".fa"):
            fileList.append(os.path.join(FILE_TMP_FOLDER, file))
            
    print("[INFO] Detected {} fasta files".format(len(fileList)))
    return fileList

#
# Runs clustalw on a given fasta file.
# Returns Bio.Phylo phylogenetic tree object.
#
def runClustal(filePath):
    print("[INFO] Running clustalw on {}".format(filePath))
    clustalw_cline = ClustalwCommandline("clustalw", infile=filePath)
    stdout, stderr = clustalw_cline()
    
    dndFilePath = os.path.splitext(filePath)[0]+".dnd"
    print("[INFO] Creating phylogenetic tree from {}".format(dndFilePath))
    tree = Phylo.read(dndFilePath, "newick")
    return tree

#
# Runs muscle on a given fasta file.
# Returns Bio.AlignIO object
#
def runMuscle(filePath):
    
    alnFilePath = os.path.splitext(filePath)[0]+".muscle.aln"
    print("[INFO] Running muscle on {}".format(filePath))
    
    muscle_cline = MuscleCommandline(input=filePath, out=alnFilePath, clw=True)
    stdout, stderr = muscle_cline()
    
    print("[INFO] Creating alignment from {}".format(alnFilePath))
    align = AlignIO.read(alnFilePath, "clustal")
    return align

#
# Creates NJ tree from alignment with the given distance calculator model
#
def createNJPhyloTree(align, distanceModel="identity", alignName="anonymous"):
    
    print("[INFO] Calculating distance matrix for {} alignment and {} distance model".format(alignName, distanceModel))
    calculator = DistanceCalculator(distanceModel)
    dm = calculator.get_distance(align)
    
    print("[INFO] Constructing NJ phylogenetic tree for {} alignment".format(alignName))
    constructor = DistanceTreeConstructor()
    njtree = constructor.nj(dm)
    return njtree

def saveTrees(trees, fileName, formatName="phyloxml"):
    print("[INFO] Saving phylogenetic tree/-s to file {}".format(fileName))
    Phylo.write(trees, fileName, formatName)

#
# Helper function from https://stackabuse.com/levenshtein-distance-and-text-similarity-in-python/
# to calculate Levenshtein distance (used by standarizeGraphLabeling function)
#
def levenshtein(seq1, seq2):  
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros ((size_x, size_y))
    for x in xrange(size_x):
        matrix [x, 0] = x
    for y in xrange(size_y):
        matrix [0, y] = y

    for x in xrange(1, size_x):
        for y in xrange(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix [x,y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
    return (matrix[size_x - 1, size_y - 1])

#
# Perform individual label standarization
#   * Characters "(" and ")" are replaced with "_"
#   * Inner nodes' labels are removed
#
def normalizeGraphLabeling(G, GLabel='G1'):
    print("[INFO] Normalizing labels on {}".format(GLabel))
    removedInnerNodeLabelsCnt = 0
    changedLeafLabelsCnt = 0
    
    for node in G.nodes():
        if G.out_degree(node)==0 and G.in_degree(node)==1:
            if node.name:
                name = node.name.replace("(", "_").replace(")", "_")
                if name != node.name:
                    changedLeafLabelsCnt = changedLeafLabelsCnt+1
                node.name = name
        else:
            if node.name:
                removedInnerNodeLabelsCnt = removedInnerNodeLabelsCnt+1
            node.name = None
    print("[INFO] Inner nodes labels removed: {}".format(removedInnerNodeLabelsCnt))
    print("[INFO] Leaf nodes labels changed: {}".format(changedLeafLabelsCnt))

#
# This function forks on Phylo trees (not the hashed graphs)
# And checks if the labelling of two trees exactly match if not then the exception is raised
#
def checkIfPhyloTreeLabelingMatch(T1, T2):
    startingNode1 = T1.clade
    dfsTree1 = nx.dfs_tree(Phylo.to_networkx(T1), source=startingNode1)
    nameSet1 = set([x.name for x in dfsTree1.nodes()])
    startingNode2 = T2.clade
    dfsTree2 = nx.dfs_tree(Phylo.to_networkx(T2), source=startingNode2)
    nameSet2 = set([x.name for x in dfsTree2.nodes()])
    if None in nameSet1:
        nameSet1.remove(None)
    if None in nameSet2:
        nameSet2.remove(None)
    leftDiff = nameSet1 - nameSet2
    rightDiff = nameSet2- nameSet1
    diffCount = len(leftDiff) + len(rightDiff)
    if diffCount > 0:
        raise Exception("Phylo tree labelling differs. Found {} mismatches.".format(diffCount))
    
#
# Muscle trees differs from custalw trees in labelling
# This function takes two hashed graphs (provided by createHashedPTreeGraph)
# and modifies labelling to standarized one.
#
# * Labels are fixed individually using normalizeGraphLabeling
# * Inner labels are removed
# * Labels for leafs from the second tree that are not present if the first one
#     are changed to the nearest labels existing in the first tree (using Levenshtein metric)
# * All labels left after those two steps that are in G1 and not in G2 (or in G2 and not in G1) are removed
#
#
# The Levenshtein metric is used because ClustalW and Muscle seem to provide labels that are cut off names
# without single letters or some other odd artifacts.
#
def standarizeGraphLabeling(G1, G2, G1Label='G1', G2Label='G2'):
    print("[INFO] Performing label standarization for graphs {} and {}".format(G1Label, G2Label))
    # Normalize labels individually due to standard rules
    normalizeGraphLabeling(G1, GLabel=G1Label)
    normalizeGraphLabeling(G2, GLabel=G2Label)
    
    # Determine sets of labels of nodes
    G1NameSet = set([x.name for x in G1.nodes()])
    G2NameSet = set([x.name for x in G2.nodes()])
    if None in G1NameSet:
        G1NameSet.remove(None)
    if None in G2NameSet:
        G2NameSet.remove(None)
    # Create differences sets
    leftDiff = G1NameSet - G2NameSet
    rightDiff = G2NameSet - G1NameSet
    diffCount = len(leftDiff) + len(rightDiff)

    print("[INFO] Found {} differences in labelling".format(diffCount))
    if diffCount <= 0:
        print("[INFO] Labelling is standarized now. Nothing to do.")
        return True
    
    # Use Levenshtein metric to store new node labels mapping inside labelMapping dict
    labelMapping = {}
    for name in leftDiff:
        if not (name in labelMapping):
            lowestScore = 1000000
            lowestLabel = None
            for name2 in rightDiff:
                score = levenshtein(name, name2)
                if score < lowestScore:
                    lowestScore = score
                    lowestLabel = name2
            if lowestLabel:
                labelMapping[lowestLabel] = name
    for name in rightDiff:
        if not (name in labelMapping):
            lowestScore = 1000000
            lowestLabel = None
            for name2 in leftDiff:
                score = levenshtein(name, name2)
                if score < lowestScore:
                    lowestScore = score
                    lowestLabel = name2
            if lowestLabel:
                labelMapping[lowestLabel] = name
    
    # Perform relabelling using labelMapping dict
    for node in G1.nodes():
        if node.name in labelMapping:
            node.name = labelMapping[node.name]
            
    for node in G2.nodes():
        if node.name in labelMapping:
            node.name = labelMapping[node.name]
            
    # Recalculate labels sets
    G1NameSet = set([x.name for x in G1.nodes()])
    G2NameSet = set([x.name for x in G2.nodes()])
    if None in G1NameSet:
        G1NameSet.remove(None)
    if None in G2NameSet:
        G2NameSet.remove(None)
    leftDiff = G1NameSet - G2NameSet
    rightDiff = G2NameSet - G1NameSet
    diffCount = len(leftDiff) + len(rightDiff)
    
    print("[INFO] Found {} differences in labelling after applied fixing with Levenshtein distance".format(diffCount))
    if diffCount <= 0:
        print("[INFO] Labelling is standarized now. Nothing to do.")
        return True
    else:
        print("[INFO] Will remove unfixed labels.")
        
    # If Levenshtein relabeling is not enough
    # then we remove rest of the unmatched labels
    for node in G1.nodes():
        if (node.name in leftDiff) or (node.name in rightDiff):
            node.name = None
            
    for node in G2.nodes():
        if (node.name in leftDiff) or (node.name in rightDiff):
            node.name = None
    
    # Recalculate labels sets
    G1NameSet = set([x.name for x in G1.nodes()])
    G2NameSet = set([x.name for x in G2.nodes()])
    if None in G1NameSet:
        G1NameSet.remove(None)
    if None in G2NameSet:
        G2NameSet.remove(None)
    leftDiff = G1NameSet - G2NameSet
    rightDiff = G2NameSet - G1NameSet
    diffCount = len(leftDiff) + len(rightDiff)
    
    # If the labelling is not matched then something bad happened, so we raise an error
    if diffCount <= 0:
        print("[INFO] Labelling is standarized now. Nothing to do.")
        return True
    else:
        raise Exception("Could not fix graph labelling")
        return False

#
# Creates hashed graph tree from Phylo tree.
# The hashed graph format is a graph in NetworkX format (directed DFS tree)
# with locally unique ids (unique for that graph)
# for nodes, standarized labelling for trees (if refTree param was provided) and hashes of
# bipartition graphs induced by the tree edges.
#
# The function acccepts the following parameters:
#   tree - Phylo tree to generate hashed tree graph from
#   maxLeafHashValue - Maximum number hash value for nodes
#   refTree - if provided it's the tree that will share labelling with the current one
#     (labelling standarization will be performed to make both trees share common labels)
#   treeName - provide tree debug name to be displayed in logs
#   hashAlgorithm - provide hash algorithm for nodes
#
# Supported hash algorithms are the following ones:
#   - pow - Hash is always 2^n where n are unique for each leaf
#        maxLeafHashValue value will be ignored and used only in case when node hash is missing
#       (exotic behaviour that should never happen)
#   - rand - Randomized hashes in range [0..maxLeafHashValue]
#
# Default hashing algorithm is "pow"
# Hashing is used to hash bipartitions of graph.
#
def createHashedPTreeGraph(tree, maxLeafHashValue=2000000000000000000000, refTree=None, treeName="anonymous", hashAlgorithm="pow"):
    
    print("[INFO] Constructing hashed PTree graph via networkx for {} tree".format(treeName))
    
    if refTree is None:
        print("[INFO] Using purely new node ids for nodes")
    else:
        print("[INFO] Using reference tree ids for nodes")
    
    # To each node we assign ids that are integers starting from 0
    def assignPhyloIds(tree, nextFreeId=0):            
        tree.id = nextFreeId
        nextFreeId = nextFreeId+1
        
        for child in tree:
            nextFreeId = assignPhyloIds(child, nextFreeId)
        return nextFreeId
    
    print("[INFO] Assigning random leafs hashes in range 0 to {}".format(maxLeafHashValue))
    assignPhyloIds(tree.clade)
    net = Phylo.to_networkx(tree)
    
    # Create directed DFS traversal graph for tree
    startingNode = tree.clade
    dfsTree = nx.dfs_tree(net, source=startingNode)
    if not (refTree is None):
        standarizeGraphLabeling(dfsTree, refTree["network"])
        
    nodesCount = len(dfsTree.nodes())
    
    # Assign hashes to all leafs
    leafId = 1
    nodeToHashMapping = {}
    hashToNodeMapping = {}
    for node in dfsTree.nodes():
        if len(node) <= 0:
            
            if hashAlgorithm == "pow":
                leafHash = 2 ** leafId
            elif hashAlgorithm == "rand":
                nodeHash = random.randint(1, maxLeafHashValue+1)
                
            leafId = leafId + 1
            if not (refTree is None):
                if not (node.name is None):
                    foundNodes = list(it.ifilter(lambda n: n.name == node.name, refTree["network"].nodes()))
                    if len(foundNodes) > 0:
                        leafHash = refTree["nodeToHashMapping"][foundNodes[0].id]["hash"]
            
            hashObj = {
                "clade": node,
                "hash": leafHash
            }
            nodeToHashMapping[node.id] = hashObj
            hashToNodeMapping[leafHash] = hashObj
    
    print("[INFO] Generating hashes for the rest of nodes")
    
    # Function that recursively generates hashes for nodes from
    # hashes of the children nodes
    def recAssignHashes(node, parent):
        
        for child in nx.neighbors(dfsTree, node):
            recAssignHashes(child, node)
        
        if not node.id in nodeToHashMapping:
            nodeHash = None
            okChild = 0
            for child in nx.neighbors(dfsTree, node):
                if child.id in nodeToHashMapping:
                    okChild = okChild + 1
                    if nodeHash is None:
                        nodeHash = nodeToHashMapping[child.id]["hash"]
                    else:
                        if hashAlgorithm == "pow":
                            nodeHash = nodeHash | nodeToHashMapping[child.id]["hash"]
                        elif hashAlgorithm == "rand":
                            nodeHash = nodeHash ^ nodeToHashMapping[child.id]["hash"]

            if nodeHash is None:
                nodeHash = random.randint(1, maxLeafHashValue+1)
            
            hashObj = {
                "clade": node,
                "hash": nodeHash
            }
            nodeToHashMapping[node.id] = hashObj
            hashToNodeMapping[nodeHash] = hashObj
        
    recAssignHashes(startingNode, None)
    return {
        "network": dfsTree,
        "nodeToHashMapping": nodeToHashMapping,
        "hashToNodeMapping": hashToNodeMapping
    }

#
# Function to calculate RF distance between two trees
# Returns the sum over A, B of (d(T1, A, B) - d(T2, A, B))^2
#
def RFDistance(t1, t2):
    print("[INFO] Calculating RF distance for the given trees")
    result = 0
    
    # Function to evaluate None to 0
    def defaultNone(val):
        if val is None:
            return 0
        return val
    
    lenDifs = []
    
    print("[INFO] Calculating RF-left distance for T1")
    for edge in t1["network"].edges():
        node = edge[1]
        hashNodeT1 = t1["nodeToHashMapping"][node.id]
        if (hashNodeT1["hash"] in t2["hashToNodeMapping"]):
            hashNodeT2 = t2["hashToNodeMapping"][hashNodeT1["hash"]]
            lenDifs.append((defaultNone(hashNodeT1["clade"].branch_length), defaultNone(hashNodeT2["clade"].branch_length)))
        else:
            lenDifs.append((defaultNone(hashNodeT1["clade"].branch_length), 0))

    print("[INFO] Calculating RF-right distance for T2")
    for edge in t2["network"].edges():
        node = edge[1]
        hashNodeT2 = t2["nodeToHashMapping"][node.id]
        if (hashNodeT2["hash"] in t1["hashToNodeMapping"]):
            hashNodeT1 = t1["hashToNodeMapping"][hashNodeT2["hash"]]
        else:
            lenDifs.append((0, defaultNone(hashNodeT2["clade"].branch_length)))
    
    lenDifs = [(round(i[0], 5), round(i[1], 5)) for i in lenDifs]
    for i in lenDifs:
        i = (round(i[0], 5), round(i[1], 5))
        result = result + (i[0] - i[1])*(i[0] - i[1])
    return result

#
# Helper function to calculate RF distance for trees using DendroPy
# DendroPy does not support (d1-d2)^2 metric - only |d1-d2| and sqrt(d1^2-d2^2) are supported.
# So we provide our own df lambda and use DendroPy internals
#
def custom_distance(tree1, tree2):
    df_helper = lambda length_diffs: (sum([(i[0] - i[1])*(i[0] - i[1]) for i in length_diffs]))
    return treecompare._bipartition_difference(tree1,tree2,dist_fn=df_helper,edge_weight_attr="length",value_type=float,is_bipartitions_updated=False)

#
# Function to calulate RF distance with help of DendroPy
# For given Nexus-format files
#
def refRFDistance(t1NexFilePath, t2NexFilePath):
    tns = TaxonNamespace()
    nexTree1 = Tree.get(unconstrained_taxa_accumulation_mode=True, path=t1NexFilePath, schema="nexus", taxon_namespace=tns)
    nexTree2 = Tree.get(unconstrained_taxa_accumulation_mode=True, path=t2NexFilePath, schema="nexus", taxon_namespace=tns)
    return (custom_distance(nexTree1, nexTree2))

#
# MAIN CODE BLOCK
#
@click.command()
@click.argument('url')
@click.option('--tmp', '-t', default=FILE_TMP_FOLDER, help="Path to save temporary data. Default is {}".format(FILE_TMP_FOLDER))
@click.option('--disablecheck', '-d', type=bool, default=False, help="Disable DendroPy verification check")
@click.option('--hashing', '-h', default="pow", help="Change default hashing strategy (can be 'pow' or 'rand')")
def main(url, tmp, disablecheck, hashing):
    global FILE_TMP_FOLDER
    global FILE_TO_DOWNLOAD
    FILE_TO_DOWNLOAD = url
    
    if tmp:
        FILE_TMP_FOLDER = tmp
        print("Using tmp folder {}".format(FILE_TMP_FOLDER))

    if (hashing <> "pow") and (hashing <> "rand"):
        raise Exception("Invalid hashing algorithm was specified: {}".format(hashing))
        
    for fastaFile in downloadFiles(FILE_TO_DOWNLOAD):
        # Create first tree
        T1tree = runClustal(fastaFile)
        # Create second tree
        muscleAlign = runMuscle(fastaFile)
        T2tree = createNJPhyloTree(muscleAlign, "blosum62", fastaFile)
        # Create hashed trees for T1 and T2
        T1Htree = createHashedPTreeGraph(T1tree, treeName="clustalTree", hashAlgorithm=hashing)
        T2Htree = createHashedPTreeGraph(T2tree, refTree=T1Htree, treeName="muscleTree", hashAlgorithm=hashing)
        # Check if the labelling is OK
        checkIfPhyloTreeLabelingMatch(T1tree, T2tree)
        # Save T1 and T2 to the xml file in case something bad would happen
        T1treeNexFile = os.path.splitext(fastaFile)[0] + ".t1.nex"
        T2treeNexFile = os.path.splitext(fastaFile)[0] + ".t2.nex"
        saveTrees(T1tree, T1treeNexFile, formatName="nexus")
        saveTrees(T2tree, T2treeNexFile, formatName="nexus")
        # Calculate RF distance
        dist = RFDistance(T1Htree, T2Htree)
        # Calculate reference DendroPY RF distance
        distRef = "<DISABLED BY --disable-check SWITCH>"
        if not disablecheck:
            distRef = refRFDistance(T1treeNexFile, T2treeNexFile)
        # Print result
        strResult = ("\n  | The calculated distance between trees is: {}\n  | The RF distance calculated by DendroPy is: {}\n".format(dist, distRef))
        print(strResult)
        # Save the results also to the text file
        outFile = os.path.splitext(fastaFile)[0] + ".answer.out"
        print("\n[INFO] Result was saved to {}".format(outFile))
        with open(outFile, 'w') as file:
            file.write(strResult)
    
# Entry
if __name__ == "__main__":
    main()
    