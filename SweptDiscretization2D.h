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

#ifndef H_SWEPTDISC2D
#define H_SWEPTDISC2D

#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)

#include "Swept2DGlobals.h"
#include "DebugUtils.h"
#include "mpi.h"
#include "UpPyramid.h"
#include "DownPyramid.h"
#include "HorizontalBridge.h"
#include "VerticalBridge.h"
#include "ProcessGraph.h"
#include <string>

using namespace std;
extern ProcessGraph pg;
extern string getRemoteCommand(int *code);
extern unsigned char* constantsArrayBytes;
extern void printHex(double a);

class SweptDiscretization2D
{
	friend class UpPyramid;
	friend class DownPyramid;
	friend class VerticalBridge;
	friend class HorizontalBridge;

private:
	int n;
	int panelSize;
	int substeps;
	int dataPointSize;
	int constants;
	int communicationSize;
	int foundationSize;
	double *foundation ;
	double *remoteConstants;
	int totalConservedQuantities;
	int conserveCheckFreq;
	double *staging    ;
	double *northPanel ;
	double *southPanel ;
	double *eastPanel  ;
	double *westPanel  ;
	double *resultArray;
	string sendString;
	string outputDirectory;
	int remoteCommandCode,checkForRemoteCycles;
	int outputToFile;
	int remoteConstantsCount;
	UpPyramid *up;
	DownPyramid *dp;
	HorizontalBridge *hb;
	VerticalBridge *vb;
	double *convervedVariables;	
	MPI_Win foundationWindow,constantsWindow;
	MPI_Request reqs[4];
	bool firstConservedCheck;
	int ijToIndex(int i,int j);
	int ijToConstantIndex(int i,int j);
	void printConstantFromWindow(int constNum = 0);
	void printConstantFromFoundation(int constNum = 0);
	void updateRemoteConstants(unsigned char *buffer);
	void allGatherOutputToFile(int dataPoint,string filename);
	void allGatherAllOutputToFile(string filename);
	void reportConveredQuantities();
	bool firstRun;
	int outputLength;
	
public:
	SweptDiscretization2D(int n,int substeps,int dataPointSize,int totalConstants,int totalConservedQuantities,int remoteConstantsCount = 0);
	void setInitCondition(void initFnc(int,int,InitPoint2D*));
	void updateSystem(void initFnc(int,int,InitPoint2D*));
	void getFoundation(double *result,int dataPoint = 0);
	void printFoundation(int dataPoint = 0);
	void printFullFoundation(int dataPoint = 0);
	void printFoundationHex(int dataPoint = 0);
	void printFullFoundationHex(int dataPoint = 0);
	double calculate(int cycles);	
	void setCheckForRemoteCommands(int cycles);
	void setFileOutput(int cycles);
	void setOutputDirectory(string outputDirectory);
	void setConservationCheck(int cycles);
	void setOutputLength(int outputLength);
    ~SweptDiscretization2D();
};

#endif

