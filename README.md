# Quantum_LDPC_Codes
Collection of quantum LDPC codes

Hyperbolic Codes Planar: a collection of planar hyperbolic codes ordered by Schläfli symbol {f,d} of the corresponding infinite graph.  These are CSS codes whose Z and X generators are edge-site incidence matrices of a pair of mutually dual locally planar vertex transitive graphs, finite quotients of the regular hyperbolic tiling.  
Hyperbolic Codes Q_ary: q-ary version of the hyperbolic codes, where the incidence matrix elements are 0, 1 and -1, such that the incidence matrix of a graph Gx is orthogonal to that of the dual graph Gz over Z, instead of only Z/2Z.

qudit_hyperbolic_1.cpp: To create q-ary incidence matrices based on the binary ones.
Input: Put the incidence matrices in the same directory, write the filenames of the same {p,q} tilings and the dual tilings in "latticelist" and "duallatticelist". "d" is the vertex degree of the graph, and "dd" is the degree of the dual graph. The graphs must be regular.
Output: The q-ary incidence matrices, filename is "q_" followed by the original filename. If the graph are non-orientable, e.g. on a Klein bottle, it returns -1 with an error.
