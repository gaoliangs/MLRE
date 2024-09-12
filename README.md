# MLRE(Maximum Likelihood REtrotransposon)

1. preprocessing.py
   input file: .nex
   thresthold: 
   output file: newseq1.csv

2. triplet-joining.py
   input file: newseq1.csv
   output file: tjtree_newick.txt

3. treevalue.py/NNI.py/SPR.py
   input file: tjtree_newick.txt
   input file: newseq1.csv
   option: parameter c (yes/no)
   output: treetopo & F-value
   
