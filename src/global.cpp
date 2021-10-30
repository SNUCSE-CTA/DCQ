//***************************************************************************
// This file is part of the graph isomorphism algorithm.
// Copyright by Geonmo Gu, Yehyun Nam, and Kunsoo Park
// 
// Name: global.cpp
// Author: Geonmo Gu
// Version
//     August 20, 2020: the first stable version. (version 1.0)
//***************************************************************************
#include "global.h"

vector<int32_t> global_temp_vector;
int32_t* markCell = NULL;
int32_t* markNode = NULL;
int32_t global_mark = 0;
Memory global_memory;

int32_t markLen = 0;

//note that aNumNode is the number of vertices in g1.
//ALLOCATE global variables
void initGlobal(int32_t aNumNode)
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__ << endl;
	#endif
	int32_t n2 = aNumNode * 2;
	
	if( markLen < n2 ) {
		clearGlobal(); //avoid double allocation
		markCell = new int32_t[n2];
		markNode = new int32_t[n2];
		markLen = n2;
	}

	memset(markCell, 0, sizeof(int32_t) * n2);
	memset(markNode, 0, sizeof(int32_t) * n2);
	global_mark = 0;
}

//DEALLOCATE global variables
void clearGlobal()
{
	#ifdef DEBUG
	cout << __PRETTY_FUNCTION__ << endl;
	#endif

	if( markCell != NULL ) {
		delete[] markCell;
		markCell = NULL;
	}
	if( markNode != NULL ) {
		delete[] markNode;
		markNode = NULL;
	}
	markLen = 0;
}
