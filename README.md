# MLRE(Maximum Likelihood REtrotransposon)

**Algorithm for Phylogenetic Tree Construction from Low-homoplasy Markers**  
This project introduces an advanced algorithm for constructing phylogenetic trees using retrotransposon markers. By combining maximum-likelihood estimation with heuristic search techniques, it provides high accuracy in tree reconstruction.


## Features

- **Handling ILS and Homoplasy**: Optimized to manage the challenges posed by both Incomplete Lineage Sorting (ILS) and homoplasy, ensuring more accurate tree estimation.  
- **Maximum-Likelihood Model**: A probabilistic model is used to compute maximum likelihood values for various tree topologies, ensuring a reliable fit to the data.
- **Heuristic Search Optimization**: Efficient heuristic search techniques are employed to find the optimal tree topology, even for large and complex datasets.

## Usage

### Dependencies
- **Python**: Required for running the core scripts and integrating with Matlab.
- **Matlab**: Matlab is used via the **Matlab engine** for Python. Please refer to the official [Matlab documentation](https://www.mathworks.com/help/matlab/matlab-engine-for-python.html) for installation instructions.
- The project also utilizes MATLAB MEX files for computational efficiency. Please ensure a compatible C++ compiler is installed, and be mindful of version compatibility. All tests were conducted using MATLAB 2023b and Xcode 16.2.

### Algorithm Workflow


#### 1. Data Preprocessing

The input data for the algorithm is in `.nex` format. To preprocess the data, run the `preprocessing.py` script. This will generate a `.csv` file containing a presence/absence matrix for retrotransposon markers (i.e., a 0-1 matrix).

During execution, the script will display the weights of all Buneman clusters. You may choose to apply a threshold based on these weights to filter the clusters. Clusters with weights exceeding the threshold will be retained, resulting in a partially unresolved tree. For internal nodes with a degree greater than 3, a separate CSV file will be generated for each node to facilitate further analysis.

If you do not intend to use Buneman clustering, you may skip the threshold input by pressing Enter without providing a value. In this case, the Buneman clustering step will be bypassed, and no filtering based on cluster weights will occur.

```
Input file: .nex  
Threshold: buneman weight/ Press Enter to skip clustering
Output files: .csv
```


#### 2. Tree Search Using NNI and SPR

Run the `treesearch.py` script. This will perform a NNI search first, followed by a SPR (Subtree Pruning and Regrafting) search. 

```
Initial tree:
   The default initial tree is a Triplet-Joining tree. Simply press Enter to use the default.
   If you wish to use a custom tree, input a valid Newick format tree.
Constraint tree:
   If you do not require a constraint tree, simply press Enter.
   If you need to use a constraint tree, you can: Input buneman to use the default Buneman tree. Alternatively, input a custom Newick format tree.
Marker file: .csv 
Option:
  parameter c (yes/no)  
  parameter q (yes/no) (Hint: The parameter q is not recommended for cases with a large number of taxa.)
```


#### 3. Calculating the Maximum Likelihood Value for a Specific Tree

If you only want to calculate the maximum likelihood value for a specific tree topology without performing a tree search, you can run the `treevalue.py` script. 

Provide one or more trees in Newick format along with the retrotransposon marker matrix in `.csv` format as input. You can also choose whether to use parameter `c` and/or parameter `q`. The script will output the branch lengths (plus insertion rates/missing rates) and the Likelihood value for each input tree.

```
input tree: Newick string or filename(.txt)
input marker file: .csv 
option: parameter c (yes/no)  
        parameter q (yes/no)     
output: treetopo with parameters & Likelihood value
```


#### 4. Support values
This method uses the aBayes approach to calculate the support values for each edge in the tree. To compute these values, run the `supportvalue.py` script, which will evaluate the support for each edge in the phylogenetic tree.
```
input tree: Newick string or filename(.txt)
input file: .csv
option: parameter c (yes/no)
        parameter q (yes/no)
output: support value list
```

#### 5. Forward Simulation
The Forward Simulation tool allows users to simulate markers based on phylogenetic trees in Newick format. 
To run the simulation, execute the `forward_simulation.py` script with the following parameters:
```
input samples number: the number of replications
input markers number: the number of markers per replication
input tree: Newick string or filename
optional input marker file: CSV file for replicating missing data pattern (?).  If you do not need to simulate missing data, press Enter to skip.
output: Simulated markers in CSV format
```



