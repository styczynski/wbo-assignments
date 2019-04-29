The script description is provided in wbo.py file.
To run the script please use the following commands:
```

   $ pip install click urllib gunzip biopython networkx numpy dendropy
   $ python wbo.py "http://regulomics.mimuw.edu.pl/wp/wp-content/uploads/2018/05/homeo2.fa_.gz"
   
```
The result of running above commands is like follows:

```

    root@ip-172-31-25-39:/home/ubuntu/wbo# python wbo.py "http://regulomics.mimuw.edu.pl/wp/wp-content/uploads/2018/05/homeo2.fa_.gz"
    Using tmp folder ./downloads
    [INFO] Clearing temporary download directory
    [INFO] Downloading file http://regulomics.mimuw.edu.pl/wp/wp-content/uploads/2018/05/homeo2.fa_.gz from remote source
    [INFO] Unpacking gz archive
    [INFO] Detected 1 fasta files
    [INFO] Running clustalw on ./downloads/downloaded.fa
    [INFO] Creating phylogenetic tree from ./downloads/downloaded.dnd
    [INFO] Running muscle on ./downloads/downloaded.fa
    [INFO] Creating alignment from ./downloads/downloaded.muscle.aln
    [INFO] Calculating distance matrix for ./downloads/downloaded.fa alignment and blosum62 distance model
       [INFO] Constructing NJ phylogenetic tree for ./downloads/downloaded.fa alignment
    [INFO] Constructing hashed PTree graph via networkx for clustalTree tree
    [INFO] Using purely new node ids for nodes
    [INFO] Assigning random leafs hashes in range 0 to 2000000000000000000000
    [INFO] Generating hashes for the rest of nodes
    [INFO] Constructing hashed PTree graph via networkx for muscleTree tree
    [INFO] Using reference tree ids for nodes
    [INFO] Assigning random leafs hashes in range 0 to 2000000000000000000000
    [INFO] Performing label standarization for graphs G1 and G2
    [INFO] Normalizing labels on G1
    [INFO] Inner nodes labels removed: 837
    [INFO] Leaf nodes labels changed: 34
    [INFO] Normalizing labels on G2
    [INFO] Inner nodes labels removed: 0
    [INFO] Leaf nodes labels changed: 0
    [INFO] Found 4 differences in labelling
    [INFO] Found 0 differences in labelling after applied fixing with Levenshtein distance
    [INFO] Labelling is standarized now. Nothing to do.
    [INFO] Generating hashes for the rest of nodes
    [INFO] Saving phylogenetic tree/-s to file ./downloads/downloaded.t1.nex
    [INFO] Saving phylogenetic tree/-s to file ./downloads/downloaded.t2.nex
    [INFO] Calculating RF distance for the given trees
    [INFO] Calculating RF-left distance for T1
    [INFO] Calculating RF-right distance for T2

      | The calculated distance between trees is: 0.3524015274
      | The RF distance calculated by DendroPy is: 0.3524015274


    [INFO] Result was saved to ./downloads/downloaded.answer.out
    
```  
   
Piotr Styczyński (ps386038)