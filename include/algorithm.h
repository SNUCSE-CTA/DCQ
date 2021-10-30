//***************************************************************************
// This file is part of the graph isomorphism algorithm.
// Copyright by Geonmo Gu, Yehyun Nam, and Kunsoo Park
//
// Name: algorithm.h
// Author: Geonmo Gu, Yehyun Nam
// Version
//     August 20, 2020: the first stable version. (version 1.0)
//***************************************************************************

#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <iostream>

#include "graph.h"
#include "refine.h"
#include "backtrack.h"
#include "global.h"

// extern int num_check_colorings;
// extern int num_backtracks;

class Algorithm
{
public:
	Algorithm();
	~Algorithm();

	//parameters: [aG1], [aG2]
	//VERIFY whether aG1 and aG2 are isomorphic
	//RETURN true if aG1 and aG2 are isomorphic, false otherwise
	bool run(Graph*, Graph*);

private:
	//parameters: [aG1], [aG2]
	//VERIFY whether aG1 and aG2 has same number of vertices for each degree and label
	//RETURN true if the condition is satisfied, false otherwise
	bool checkSimpleInvariants(Graph*, Graph*);


////////////////////////////////
//FOR SINGLE GRAPH
public:
	Backtrack bt;

	//parameters: [aG1], [aG2], [aColoring1], [aColoring2], [numTreeNode]
	//VERIFY whether aG1 and aG2 are isomorphic
	//RETURN true if aG1 and aG2 are isomorphic, false otherwise
	bool run(Graph*, Graph*, Coloring*, Coloring*, int32_t);

	bool checkColorings(Graph*, Graph*, Coloring*, Coloring*);
	
	//parameters: [aG1], [aG2], [aColoring1], [aColoring2], [u], [v]
	//DETERMINE wheter the trees rooted in aG1 and aG2 rooted at u and v, respectively, are isomorphic
	//RETURN true if isomorphic, false otherwise
	bool treeIsomorphism(Graph*, Graph*, Coloring*, Coloring*, int32_t, int32_t);
};

#endif

