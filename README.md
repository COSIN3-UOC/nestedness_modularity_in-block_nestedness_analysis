# Structural Analysis Nestedness, Modularity and In-block

This code performs structural analysis in binary unipartite and bipartite networks by means of in-block nestedness
as defined by Solé-Ribalta et al, PRE 2018 (https://doi.org/10.1103/PhysRevE.97.062302), as well as global nestedness and modularity.
        
## Inputs:
       
1) `path-folder` =  directory where the network data files are. The data files with .csv file extension.
It can be either in edge list or adjacency matrix format with no headers. 
2) `bipartite` =  boolean value to indicate if "filename" is a bipartite (True) or unipartite (False) network.
3) `edge_data` = boolean value indicating the format of the data file. Three-column/edge list (True) or matrix format (False).

## Outputs:

1) A summary file, `data_structures_NQI_results.csv`, containing as many rows as input files found in `path-folder`. 
 Its columns report for each network analysed the name of its input file, and
the values of ![formula](https://render.githubusercontent.com/render/math?math=\mathcal{N}),
![formula](https://render.githubusercontent.com/render/math?math=Q), 
and ![formula](https://render.githubusercontent.com/render/math?math=\mathcal{I}).

2) As many files as input files with general naming `in-block_partition_<filename>.csv`, where `<filename>` is the extensionless name of the respective input file.
 Each files are formed by one or two columns, depending whether the network is unipartite or bipartite.
A column contains the indices of the block each node is assigned to by optimizing the in-block nestedness, following zero indexing.
Therefore, columns' lengths equal the number of nodes present in the corresponding graph or subgraph.
As an example for a bipartite network, the i-th value of the column labelled by `rows` (`columns`) 
indicates the block of the `i-th` node of the subgraph associated to the rows (columns) of the input adjacency matrix.
**N.B.** block indices are generally not contiguous, namely, the number of blocks does not in principle coincide with the largest block index+1.   

3) Like above, but the general naming is `modularity_partition_<filename>.csv` and partitions are calculated by optimising the modularity.
	
If for example we pass bipartite=True and edge_data=True, all the networks we want to analyze have to fulfill such conditions.

### example: 
```
python structural_analysis.py home/User/data/ True False

```
# Nestedness and In-block nestedness

Both metrics are computed considering the condition ![formula](https://render.githubusercontent.com/render/math?math=k_i>=k_j) to compare the paired overlap between pairs of nodes, in contrast with the NODF metrics that considers ![formula](https://render.githubusercontent.com/render/math?math=k_i>k_j)

# Modularity and in-block nestedness optimization

The optimization of modularity and in-block nestedness is performed by employing the Extremal optimization algorithm (https://doi.org/10.1103/PhysRevE.72.027104).
The main code for this function was written in c++ and should be compiled as a file with a .so extension for Python 3.x. This file will be imported in Python as a library. 

This will be possible for MacOS or Linux.


## System Requirements 	
### Compilers 

1) Clang/LLVM 3.3 or newer (for Apple Xcode's clang, this is 5.0.0 or newer) or
2) GCC 4.8 or newer

Python 3.x.

pybind11 (via pip or conda)

### To compile the c++ files 

#### Compilation on Linux: 
```
g++ -O3 -Wall -shared -std=c++11 -fPIC `python -m pybind11 --includes` EO_functions_bipartite.cpp -o extremal_bi.so
```
	
#### Compilation on MacOS: 
```
g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python -m pybind11 --includes` EO_functions_bipartite.cpp -o extremal_bi.so
```

# Citations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4557009.svg)](https://doi.org/10.5281/zenodo.4557009)

MJ Palazzi, J Borge-Holthoefer, CJ Tessone and A Solé-Ribalta. Macro- and mesoscale pattern interdependencies in complex networks. J. R. Soc. Interface, 16, 159, 20190553 (2019). DOI: [10.1098/rsif.2019.0553](https://doi.org/10.1098/rsif.2019.0553)

MJ Palazzi, A Solé-Ribalta, V Calleja-Solanas, CA Plata, S Meloni, S Suweis and J Borge-Holthoefer. An ecological approach to structural flexibility in online communication systems. Nature Communications, 12, 1941 (2021). DOI: [10.1038/s41467-021-22184-2](https://doi.org/10.1038/s41467-021-22184-2)

