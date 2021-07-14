# DCQ
Efficient Graph Isomorphism Query Processing using Degree Sequences and Color-Label Distributions

## Binary File
- DCQ: The binary file of DCQ for linux (CentOS release 6.5)

## Usage
Options:
```
-d,   specify the file containing data graphs
-q,   specify the file containing query graphs
```

To run the program:
```
./DCQ -d [data file] -q [query file]
```

## Running Example
```
./DCQ -d data/amazon_D4000_d2.igraph -q query/amazon_d2.igraph
./run_amazon.sh   //a batch script for amazon
./run_hprd.sh     //a batch script for HPRD
```

## Results
DCQ outputs indexing time and query processing time in milliseconds, 
and the answer file for a set of data graphs and a query set (answer.txt).

Each line of the answer file is as follows.
```
[query ID]: [data graph ID] [data graph ID] ...
```

## Input File Format
The file format is a text format to store an undirected graph. 
A file contans a set of data graphs (or a set of query graphs).
- The first line of the file should be "t ID #vertices" which means the start of a graph with graph id=ID
- Following lines of "v [vertex ID] [vertex label]" indicate the vertices in the graph.
- The vertices should be written in the file in ascending order of their IDs, and a vertex ID should be in range [0, #vertices - 1].
- Following lines of "e [vertex ID] [vertex ID] 0" after the vertices indicate the undirected edges in the graph.

For example:
```
Line "t 1 3112" means that the start of a graph with ID=1 and #vertices=3112.
Line "v 0 1" means that there is a vertex with ID=0 and vertex-label=1 in the graph.
Line "e 1 133 0" means that there is an undirected edge between vertices with IDs 1 and 133.
```

## Licensing
This project is provided under the Apache License 2.0.
