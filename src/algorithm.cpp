//***************************************************************************
// This file is part of the graph isomorphism algorithm.
// Copyright by Geonmo Gu, Yehyun Nam, and Kunsoo Park
// 
// Name: algorithm.cpp
// Author: Geonmo Gu, Yehyun Nam
// Version
//     August 20, 2020: the first stable version. (version 1.0)
//***************************************************************************

#include "algorithm.h"

#include <algorithm>

// int num_check_colorings = 0;
// int num_backtracks = 0;

Algorithm::Algorithm() {}

Algorithm::~Algorithm() {}

//VERIFY whether aG1 and aG2 are isomorphic
//RETURN true if aG1 and aG2 are isomorphic, false otherwise
bool Algorithm::run(Graph* aG1, Graph* aG2)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__ << endl;
	#endif

	initGlobal(aG1->numNode);

	// bool result = checkSimpleInvariants(aG1, aG2);
	// if (result == false) {
	// 	return false;
	// }

	Refinement cr;
	bool result = cr.run(aG1, aG2);
	if (result == false) {
		return false;
	}

	Coloring* coloring = cr.getStableColoring();

	Backtrack bt;
	result = bt.run(coloring, aG1, aG2, cr.getNumTreeNode());

	clearGlobal();

	return result;
}

//VERIFY whether aG1 and aG2 has same number of vertices for each degree and label
//RETURN true if the condition is satisfied, false otherwise
bool Algorithm::checkSimpleInvariants(Graph* aG1, Graph* aG2)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__ << endl;
	#endif

	if (aG1->numNode != aG2->numNode)
		return false;
	if (aG1->numEdge != aG2->numEdge)
		return false;

	const int32_t numNode = aG1->numNode;

	int32_t* nodes1 = global_memory.getLLArray(numNode);
	int32_t* nodes2 = global_memory.getLLArray(numNode);

	for (int32_t i = 0; i < numNode; ++i)
		nodes1[i] = nodes2[i] = i;

	sort(nodes1, nodes1 + numNode,
			[aG1](const int32_t& i, const int32_t& j) -> bool {
				if (aG1->d[i] == aG1->d[j])
					return aG1->l[i] < aG1->l[j];

				return aG1->d[i] < aG1->d[j];
			});

	sort(nodes2, nodes2 + numNode,
			[aG2](const int32_t& i, const int32_t& j) -> bool {
				if (aG2->d[i] == aG2->d[j])
					return aG2->l[i] < aG2->l[j];

				return aG2->d[i] < aG2->d[j];
			});

	for (int32_t i = 0; i < numNode; ++i) {
		if (aG1->d[nodes1[i]] != aG2->d[nodes2[i]])
			return false;
		if (aG1->l[nodes1[i]] != aG2->l[nodes2[i]])
			return false;
	}

	global_memory.returnLLArray(nodes1, numNode);
	global_memory.returnLLArray(nodes2, numNode);

	return true;
}


////////////////////////////////
//FOR SINGLE GRAPH


//VERIFY whether aG1 and aG2 are isomorphic
//RETURN true if aG1 and aG2 are isomorphic, false otherwise
bool Algorithm::run(Graph* aG1, Graph* aG2, Coloring* aColoring1, Coloring* aColoring2, int32_t aNumTreeNode)
{
	//simple invariant check
	if(aG1->numNode != aG2->numNode || aG1->numEdge != aG2->numEdge) {
		return false;
	}

    // ++num_check_colorings;

	bool result = checkColorings(aG1, aG2, aColoring1, aColoring2);
	if( result == false ) {
//cout << "checkColoring FALSE!" << endl;
		return false;
	}

    // ++num_backtracks;

	//allocate memory for  CS, etc
	bt.initSingle(aG1, aG2, aColoring1, aColoring2, aNumTreeNode);

	//build DAG only once
	if(bt.dag == NULL) {
		bt.dag = new DAG;
		bt.dag->numNode = bt.n;
		bt.dag->numEdge = bt.e;
		bt.dag->childSize = new int32_t[bt.n];
		bt.dag->parentSize = new int32_t[bt.n];
		bt.dag->dagArr = new int32_t[bt.e2];
		bt.buildDAGSingle();
	}
	bt.buildCSSingle();

	int32_t numMatching = 0;
	numMatching = bt.mapBinaryCellSingle();

	if( bt.failingset == NULL )
		bt.failingset = new vector<int32_t>[aG1->numNode - numMatching - aNumTreeNode];
	result = bt.backtrack(numMatching + aNumTreeNode);
	if( bt.failingset != NULL) {
		delete[] bt.failingset;
		bt.failingset = NULL;
	}

	bt.clearSingle();

	if(result == false)
		cout << "backtrack FALSE!" << endl;

	return result;
}


//check if the union of the colorings is the coarsest stable coloring of the union of the graphs
bool Algorithm::checkColorings(Graph* aG1, Graph* aG2, Coloring* aColoring1, Coloring* aColoring2)
{
	if(aColoring1->numCell != aColoring2->numCell) {
//cout << "numCell different: " << aColoring1->numCell << " " << aColoring2->numCell << endl;
		return false;
	}
/*
cout << "COLORING 1" << endl;
cout << "PERM:";
for(int32_t i = 0; i < aColoring1->numNode; ++i)
	cout << aColoring1->perm[i] << " ";
cout << endl;
cout << "COLOR:";
for(int32_t i = 0; i < aColoring1->numNode; ++i)
	cout << aColoring1->color[i] << " ";
cout << endl;

cout << "COLORING 2" << endl;
cout << "PERM:";
for(int32_t i = 0; i < aColoring2->numNode; ++i)
	cout << aColoring2->perm[i] << " ";
cout << endl;
cout << "COLOR:";
for(int32_t i = 0; i < aColoring2->numNode; ++i)
	cout << aColoring2->color[i] << " ";
cout << endl;
*/
	//for each cell of aC1 and aC2, check that the numbers of vertices are the same
	for(int32_t i = 0; i < aColoring1->numNode; i += aColoring1->cellSize[i]) {
		if(aColoring1->cellSize[i] != aColoring2->cellSize[i]) {
//cout << "cellSize different: " << aColoring1->cellSize[i] << " " << aColoring2->cellSize[i] << endl;
			return false;
		}
		if(aG1->one[ aColoring1->perm[i] ] != aG2->one[ aColoring2->perm[i] ])
			return false;
	}
	for(int32_t i = 0; i < aColoring2->numNode; i += aColoring2->cellSize[i]) {
		if(aColoring2->cellSize[i] != aColoring1->cellSize[i]) {
//cout << "cellSize different: " << aColoring1->cellSize[i] << " " << aColoring2->cellSize[i] << endl;
			return false;
		}
		if(aG2->one[ aColoring2->perm[i] ] != aG1->one[ aColoring1->perm[i] ])
			return false;
	}

	int32_t* count = new int32_t[2 * aColoring1->numNode];
	memset(count, 0, sizeof(int32_t) * 2 * aColoring1->numNode);
	char* visited = new char[aColoring1->numNode];
	memset(visited, 0, sizeof(char) * aColoring1->numNode);
	vector<int32_t> checkList;
	int padn = aColoring1->numNode;

	//for each cell c
	for(int32_t c = 0; c < aColoring1->numNode; c += aColoring1->cellSize[c]) {
//cout << "c: " << c << endl;
		//for each vertex u in c (aG1)
		for(int32_t i = 0; i < aColoring1->cellSize[c]; ++i) {
			int32_t u = aColoring1->perm[c + i];
//cout << "\t" << "u: " << u << endl;
			//for each neighbor v of u
			for(int32_t j = 0; j < aG1->d[u]; ++j) {
				int32_t v = aG1->e[u][j];
				if(aG1->one[v]) //we check direction from non-tree nodes to tree nodes (not the reverse)
					continue;
//cout << "\t\t" << "v: " << v << endl;
				int32_t cv = aColoring1->color[ aColoring1->inv[v] ];
				++(count[v]);
				//if cv first visited
				if(visited[cv] == 0) {
					visited[cv] = 1;
					checkList.push_back(cv);
				} //if
			} //for (j)
		} //for (i) aG1
		
		//for each vertex u in c (aG2)
		for(int32_t i = 0; i < aColoring2->cellSize[c]; ++i) {
			int32_t u = aColoring2->perm[c + i];
//cout << "\t" << "u: " << u << endl;
			//for each neighbor v of u
			for(int32_t j = 0; j < aG2->d[u]; ++j) {
				int32_t v = aG2->e[u][j];
				if(aG2->one[v]) //we check direction from non-tree nodes to tree nodes (not the reverse)
					continue;
//cout << "\t\t" << "v: " << v << endl;
				int32_t cv = aColoring2->color[ aColoring2->inv[v] ];
				++(count[v + padn]);
				//if cv first visited
				if(visited[cv] == 0) {
					visited[cv] = 1;
					checkList.push_back(cv);
				} // if
			} //for (j)
		} //for (i) aG2

		//for each cell nc in checkList
		for(int32_t i = 0; i < checkList.size(); ++i) {
			int32_t nc = checkList[i];
			//check if there exists u such that count(u) is not the same as others
			int32_t u = aColoring1->perm[nc]; //the first vertex (aG1)
			int32_t theCount = count[u];
//cout << "\t" << "G1 nc: " << nc << ", theCount: " << theCount << endl;
			count[u] = 0; //reset
			//for each vertex u in nc (aG1)
			for(int32_t j = 1; j < aColoring1->cellSize[nc]; ++j) {
				u = aColoring1->perm[nc + j];
//cout << "\t\t" << "count of " << u << ": " << count[u] << endl;
				if( theCount != count[u] ) {
					delete[] count;
					delete[] visited;
//cout << "\t\t" << " count different: " << endl;
					return false;
				}
				else {
					count[u] = 0; //reset
				}
			}
//cout << "\t" << "G2 nc: " << nc << ", theCount: " << theCount << endl;
			//for each vertex u in nc (aG2)
			for(int32_t j = 0; j < aColoring2->cellSize[nc]; ++j) {
				u = aColoring2->perm[nc + j] + padn;
//cout << "\t\t" << "count of " << u << ": " << count[u] << endl;
				if( theCount != count[u] ) {
					delete[] count;
					delete[] visited;
//cout << "\t\t" << "count different: " << endl;
					return false;
				}
				else {
					count[u] = 0; //reset
				}
			}
			visited[nc] = 0; //reset
		}
		checkList.clear();
	} //for(c)


///check tree isomorphism for roots of the subtrees in G1 and G2///

	for(int32_t c = 0; c < aColoring1->numNode; c += aColoring1->cellSize[c]) {
		int32_t u = aColoring1->perm[c]; //a vertex in G1
		int32_t v = aColoring2->perm[c]; //a vertex in G2
		if(aG1->one[u]) //we are assured that aG1->one[u] == aG2->one[v]
			continue; //we skip non-roots nodes
		if(aG1->od[u] == aG1->d[u] && aG2->od[v] == aG2->d[v])
			continue; //we skip non-tree nodes

//cout << "call treeIsomorphism(" << c << ");" << endl;
		//at this point, u or v (or both of them) is a root of a subtree
		if(treeIsomorphism(aG1, aG2, aColoring1, aColoring2, u, v) == false) //check tree isomorphism between trees rooted at u and v
			return false;
	} //for(c)

//cout << " return true" << endl;
	return true;
}


//check if the subtrees in aG1 and aG2 that are rooted by a vertex colored aRootColor are isomorphic or not
bool Algorithm::treeIsomorphism(Graph* aG1, Graph* aG2, Coloring* aColoring1, Coloring* aColoring2, int32_t aRoot1, int32_t aRoot2)
{
	//NOTE: stack should contain vertices NOT COLORS 
	//(vertices with the same color may have parents with different colors, which leads to infinite loop)
	struct Quad{int32_t u; int32_t v; int32_t up; int32_t vp; };
	vector<Quad> stack;
	vector<int32_t> list;

	int32_t* count = new int32_t[aColoring1->numNode];
	memset(count, 0, sizeof(int32_t) * aColoring1->numNode);
	int32_t* rep1 = new int32_t[aColoring1->numNode];
	int32_t* rep2 = new int32_t[aColoring1->numNode];

	Quad quad = {aRoot1, aRoot2, -1, -1};
	stack.push_back(quad);

	while(stack.size() > 0) {
		list.clear();
		quad = stack.back();
		stack.pop_back();

		int32_t u = quad.u; //assured that color(u)=color(v)
		int32_t v = quad.v; 
		int32_t up = quad.up;
		int32_t vp = quad.vp;
//cout << "u: " << u << ", v: " << v << endl;

		int32_t i = aG1->d[u] < 0 ? 0 : aG1->d[u]; // d == -1 if (one[u])
		for(; i < aG1->od[u]; ++i) {
			int32_t w = aG1->e[u][i]; //for each one-neighbor w in N_G1(u)
			if(w == up) //except parent
				continue;
//cout << "one-neighbor " << w << " of N_G1(" << u << ")" << endl;
			int32_t wc = aColoring1->color[ aColoring1->inv[w] ];
			if(count[wc] == 0)
				list.push_back(wc);
			++(count[wc]);
			rep1[wc] = w;
//cout << "LINE " << __LINE__ << ", count[" << wc << "]: " << count[wc] << endl;
		} //for i in G1

		i = aG2->d[v] < 0 ? 0 : aG2->d[v];
		for(; i < aG2->od[v]; ++i) {
			int32_t w = aG2->e[v][i]; //for each one-neighbor w in N_G2(v)
			if(w == vp) //except parent
				continue;
//cout << "one-neighbor " << w << " of N_G2(" << v << ")" << endl;
			int32_t wc = aColoring2->color[ aColoring2->inv[w] ];
//cout << "LINE " << __LINE__ << ", count[" << wc << "]: " << count[wc] << endl;
			if(count[wc] == 0) {
				delete[] count;
				delete[] rep1;
				delete[] rep2;
				return false; //color distribution is not the same
			}
			--(count[wc]);
			rep2[wc] = w;
		} //for i in G2

		for(int32_t j = 0; j < list.size(); ++j) {
			int32_t cc = list[j];
			if(count[cc] != 0) {
				delete[] count;
				delete[] rep1;
				delete[] rep2;
//cout << __LINE__ << " color distribution is not the same" << endl;
				return false; //color distribution is not the same
			}
		
			quad = {rep1[cc], rep2[cc], u, v};
			stack.push_back(quad);
		} //for j

	} //while stack is not empty

	delete[] count;
	delete[] rep1;
	delete[] rep2;
	return true;
}
