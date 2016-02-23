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

#ifndef H_CLASSICDISC2D
#define H_CLASSICDISC2D

#include "Swept2DGlobals.h"
#include "DebugUtils.h"
#include "Substeps2D.h"
#include "mpi.h"
#include "ProcessGraph.h"
#include "string.h"
#include "stdio.h"

extern ProcessGraph pg;

class ClassicDiscretization2D
{
private:
	int n;
	int panelSize;
	int substeps;
	int dataPointSize;
	int constants;
	int communicationSize;
	double *foundation ;
	double *staging    ;
	double *northRegion;
	double *southRegion;
	double *eastRegion;
	double *westRegion;
	double *northRegionRecv;
	double *southRegionRecv;
	double *eastRegionRecv;
	double *westRegionRecv;
	int totalConservedQuantitie;
	bool firstRun;
	MPI_Request reqs[4];
	int ijToIndex(int i,int j);
	void initGhostRegion(double *dataSource,double *destData);
	void sendGhostRegion(double *dataSource);
	void recvGhostRegion(double *destData);
	
public:
	ClassicDiscretization2D(int n,int substeps,int dataPointSize,int constants,int totalConservedQuantitie);

	void setInitCondition(void initFnc(int,int,InitPoint2D*));
	void updateSystem(void updateFnc(int,int,InitPoint2D*));
	void getFoundation(double *result,int dataPoint = 0);
	void printFoundation(int dataPoint = 0);
	void printFullFoundation(int dataPoint = 0);
	void printFoundationHex(int dataPoint = 0);
	void printFullFoundationHex(int dataPoint = 0);
	double calculate(int steps);	
	~ClassicDiscretization2D();
};

#endif

