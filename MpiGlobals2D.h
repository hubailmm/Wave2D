// This file is part of Swept2D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu AND Maitham Alhubail, hubailmm@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) an later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT An WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef H_MPIGLOBALS2D
#define H_MPIGLOBALS2D

#include <mpi.h>
#include <algorithm>
#include <string>
#include "ProcessGraph.h"
#include <string>
#include <iostream>
#include <vector>
#include <cctype>
#ifdef _WIN32
#include <Windows.h>
#else
#include <pthread.h>
#endif

ProcessGraph pg;
string constantsArray;
string resultsArray;
vector<string> remoteCommands;
vector<int> remoteCodes;
unsigned char* constantsArrayBytes;

#ifdef _WIN32
HANDLE threadHd;
DWORD  threadId;
HANDLE remoteCommandsThreadMutexLock;
#else
pthread_t remoteCommandsThread;
pthread_mutex_t remoteCommandsThreadMutexLock;
#endif

static inline std::string &ltrim(std::string &s) 
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

static inline std::string &rtrim(std::string &s) 
{
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

static inline std::string &trim(std::string &s) 
{
	return ltrim(rtrim(s));
}

int ijToIndex(int i,int j,int xNodes,int yNodes)
{
	int ii,jj;
	ii = i;
	jj = j;
	if(i==-1)      ii = xNodes-1;
	else if(i==xNodes) ii = 0;
	if(j==-1)      jj = yNodes-1;
	else if(j==yNodes) jj = 0;

	int index = jj * xNodes + ii;
	return index;
}

ProcessGraph getProcessGraph(int myrank,int xNodes,int yNodes,int size)
{
	ProcessGraph pg;
	pg.rank = myrank;
	int j = (myrank % (xNodes*yNodes)) / xNodes;
	int i = myrank % xNodes;
	
	pg.iIndex  = i;
	pg.jIndex  = j;
	pg.xNodes  = xNodes;
	pg.yNodes  = yNodes;
	pg.mpiSize = size;
	pg.Wrank   = ijToIndex(i-1,j,xNodes,yNodes);
	pg.Erank   = ijToIndex(i+1,j,xNodes,yNodes);
	pg.Srank   = ijToIndex(i,j+1,xNodes,yNodes);
	pg.Nrank   = ijToIndex(i,j-1,xNodes,yNodes);				
	return pg;
}

void initMPI(int argc, char* argv[],int xNodes,int yNodes)
{
	int myrank,size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	pg = getProcessGraph(myrank,xNodes,yNodes,size);
}

void finMPI()
{
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

#endif
