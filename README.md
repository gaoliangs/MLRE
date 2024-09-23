# MLRE(Maximum Likelihood REtrotransposon)


---

### 1. Data Preprocessing

The input data for the algorithm is in `.nex` format. To preprocess the data, run the `preprocessing.py` script. This will generate a `.csv` file containing a presence/absence matrix for retrotransposons (i.e., a 0-1 matrix).

During execution, the script will display the weights of all Buneman clusters. You can input a threshold based on these weights to filter the clusters. 

If you do not wish to use Buneman clustering, simply set the threshold to a value larger than the maximum Buneman weight to bypass the clustering step.


You can adjust this as needed based on the specifics of your project!

   ```
   input file: .nex  
   thresthold: buneman weight  
   output file: newseq1.csv  
   ```




---

### 2. Generating the Initial Tree

To generate an initial tree for the subsequent heuristic search, use the presence/absence matrix of retrotransposons as input data and run the `triplet-joining.py` script. This will output a triplet-joining tree in Newick format.

Alternatively, you can use other software such as ASTRAL to generate the initial tree.


This provides a clear explanation of the steps for generating the initial tree.
   ```
   input file: newseq1.csv  
   output file: tjtree_newick.txt
   ```

---

### 3. Tree Search Using NNI or SPR

To perform a tree search, provide an initial tree in Newick format along with the retrotransposon marker matrix in `.csv` format as input. You can also choose whether to use parameter `c` and/or parameter `q`. 

- **NNI Search**: Run the `NNI.py` script. This will first output the logF value of the initial tree, then perform the NNI (Nearest Neighbor Interchange) search. The script will output the Newick format trees from the NNI search along with their logF values, and finally, it will output the optimal tree and its logF value after NNI search.
  
- **SPR Search**: Run the `SPR.py` script. This will perform an NNI search first, followed by an SPR (Subtree Pruning and Regrafting) search. You can choose which search method best fits your needs.


This section clearly explains the process of running tree searches using both NNI and SPR methods.

   ```
   input file: tjtree_newick.txt  
   input file: newseq1.csv  
   option: parameter c (yes/no)  
            parameter q (yes/no)    
   output: treetopo & F-value
   ```



---

### 4. Calculating the Maximum Likelihood Value for a Specific Tree

If you only want to calculate the maximum likelihood value for a specific tree topology without performing a tree search, you can run the `treevalue.py` script. 

Provide one or more trees in Newick format along with the retrotransposon marker matrix in `.csv` format as input. You can also choose whether to use parameter `c` and/or parameter `q`. The script will output the branch lengths (plus mutation rates/character deletion rates) and the F value for each input tree.


This section explains how to calculate the maximum likelihood value for a given tree topology.
   ```
   starting tree: tree(newick)  
   input file: newseq1.csv  
   option: parameter c (yes/no)  
            parameter p (yes/no)     
   output: treetopo & F-value
   ```
     
   
