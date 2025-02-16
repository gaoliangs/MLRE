# MLRE(Maximum Likelihood REtrotransposon)

**Algorithm for Phylogenetic Tree Construction from Low-homoplasy Markers**  
This project introduces an advanced algorithm for constructing phylogenetic trees using retrotransposon markers. By combining maximum-likelihood estimation with heuristic search techniques, it provides high accuracy in tree reconstruction.


## Features

- **Handling ILS and Homoplasy**: Optimized to manage the challenges posed by both Incomplete Lineage Sorting (ILS) and homoplasy, ensuring more accurate tree estimation.  
- **Initial Tree Generation**: Supports methods like Triplet-joining and linear programming for generating initial trees or using trees from existing research.
- **Maximum-Likelihood Model**: A probabilistic model is used to compute maximum likelihood values for various tree topologies, ensuring a reliable fit to the data.
- **Heuristic Search Optimization**: Efficient heuristic search techniques are employed to find the optimal tree topology, even for large and complex datasets.

## Usage

### Dependencies
- **Python**: Required for running the core scripts and integrating with Matlab.
- **Matlab**: Matlab is used via the **Matlab engine** for Python. Please refer to the official [Matlab documentation](https://www.mathworks.com/help/matlab/matlab-engine-for-python.html) for installation instructions.
- The project also utilizes Matlab MEX files for computational efficiency. Please ensure a compatible C++ compiler is installed

### Algorithm Workflow


#### 1. Data Preprocessing

The input data for the algorithm is in `.nex` format. To preprocess the data, run the `preprocessing.py` script. This will generate a `.csv` file containing a presence/absence matrix for retrotransposon markers (i.e., a 0-1 matrix).

During execution, the script will display the weights of all Buneman clusters. You may choose to apply a threshold based on these weights to filter the clusters. Clusters with weights exceeding the threshold will be retained, resulting in a partially unresolved tree. For internal nodes with a degree greater than 3, a separate CSV file will be generated for each node to facilitate further analysis.

If you do not intend to use Buneman clustering, you may skip the threshold input by pressing Enter without providing a value. In this case, the Buneman clustering step will be bypassed, and no filtering based on cluster weights will occur.

```
Input file:  .nex  
Threshold:   A number (Buneman weight) / Press Enter to skip clustering
Output files: newseq1.csv, newseq2.csv, ...
```



#### 2. Generating the Initial Tree

To generate an initial tree for the subsequent heuristic search, use the presence/absence matrix of retrotransposons as input data and run the `triplet-joining.py` script. This will output a triplet-joining tree in Newick format.

Alternatively, you can use other software such as ASTRAL to generate the initial tree.


   ```
   input file: .csv  
   output file: tjtree_newick.txt
   ```



#### 3. Tree Search Using NNI and SPR

To perform a tree search, provide an initial tree in Newick format along with the retrotransposon marker matrix in `.csv` format as input. You can also choose whether to use parameter `c` and/or parameter `q`. (Hint: The parameter `q` is not recommended for cases with a large number of taxa.)

Run the `treesearch.py` script. This will perform a NNI search first, followed by a SPR (Subtree Pruning and Regrafting) search. 
User can input a constraint tree to limit the search space. Clusters within the input tree remain fixed, with adjustments made only to other edges. Typically, a Buneman tree（`bunemantree.py`） can be used for this purpose, while a star tree can be input for a complete search.

   ```
   input file: tjtree/...
   input file: newseq1.csv  
   option: parameter c (yes/no)  
           parameter q (yes/no)
   constraint tree: buneman tree/star tree/...   
   output: treetopo & F-value
   ```




#### 4. Calculating the Maximum Likelihood Value for a Specific Tree

If you only want to calculate the maximum likelihood value for a specific tree topology without performing a tree search, you can run the `treevalue.py` script. 

Provide one or more trees in Newick format along with the retrotransposon marker matrix in `.csv` format as input. You can also choose whether to use parameter `c` and/or parameter `q`. The script will output the branch lengths (plus mutation rates/character deletion rates) and the F value for each input tree.


   ```
   starting tree: tree(newick)  
   input file: newseq1.csv  
   option: parameter c (yes/no)  
           parameter q (yes/no)     
   output: treetopo with parameters & F-value
   ```


#### 5. Support value
This method uses the aBayes approach to calculate the support values for each edge in the tree. To compute these values, run the `supportvalue.py` script, which will evaluate the support for each edge in the phylogenetic tree.
```
input tree: tree topology(newick)
input file: .csv
option: parameter c (yes/no)
        parameter q (yes/no)
output: support value list
```



