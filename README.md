# Structural Analysis Nestedness, Modularity and In-block

To perform structural analysis in binary unipartite and bipartite networks by means of nestedness, as defined by Solé-Ribalta et al, PRE 2018 (https://doi.org/10.1103/PhysRevE.97.062302), in-block nested and modularity.
    
        
## Inputs:
       
1) path-folder =  directory where the network data files are. The data files with .csv file extension. It can be either in edge list or adjacency matrix format with no headers. 
2) bipartite =  boolean to indicate if "filename" is a bipartite (True) or unipartite (False) network.
3) edge_data = boolean indicating the format of the data file. Three-column/edge list (True) or matrix format (False).
## output:
1) The function return  two .csv files for each generated network containing the label partitions from the Modularity and in-block nestednesss analysis. 
2) An additional .csv file containing a table with name, ![formula](https://render.githubusercontent.com/render/math?math=\mathcal{N}) value, ![formula](https://render.githubusercontent.com/render/math?math=Q) value, and ![formula](https://render.githubusercontent.com/render/math?math=\mathcal{I}) of value the networks
	
If for example we pass bipartite=True and edge_data=True, all the networks we want to analyze have to fulfill such conditions.

### example: 
```
python structural_analysis.py home/User/data/ True False

```
# Nestedness and In-block nestedness

Both metrics are computed considering the condition ![formula](https://render.githubusercontent.com/render/math?math=k_i>=k_j) to compare the paired overlap between pairs of nodes, in contrast with the NODF metrics that considers ![formula](https://render.githubusercontent.com/render/math?math=k_i>k_j)

# Modularity and in-block nestedness optimization

The optimization of modularity and in-block nestedness was perform by employing the Extremal optimization algorithm (https://doi.org/10.1103/PhysRevE.72.027104).
The main code for this function was written in c++ and should be compiled as a file with a .so extension for Python 3.x. This file will be imported to python as a library. 

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
Palazzi M. J., Borge-Holthoefer J., Tessone C. J. and Solé-Ribalta A. 
Macro- and mesoscale pattern interdependencies in complex networksJ. R. Soc. Interface, 16, 159, 20190553, 2019. 
[[DOI]: 10.1098/rsif.2019.0553](https://doi.org/10.1098/rsif.2019.0553)

MJ Palazzi, A Solé-Ribalta, V Calleja-Solanas, CA Plata, S Meloni, S Suweis, J Borge-Holthoefer. An Ecological Approach to Structural Flexibility in Online Communication Networks. Nature Communications (accepted, 2021).
