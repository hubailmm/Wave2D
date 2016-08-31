// This file is part of Swept2D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu AND Maitham Alhubail, hubailmm@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) an later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT Any WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef H_SUBSTEPS2D
#define H_SUBSTEPS2D

#include <math.h>
#include "Swept2DGlobals.h"

extern void printHex(double a);

class Substeps2D
{

public:

	static inline void executeStepFnc(int executeFnc,PointStruct2D *point)
	{
		const double CFL = 0.3;
		double laplacian = point->E_input[0] +  point->W_input[0] +  point->N_input[0] +  point->S_input[0] - 4 * point->C_input[0];
		point->output[1] = point->C_input[0];
		point->output[0] = 2 * point->C_input[0] - point->C_input[1] + laplacian * CFL;
	}

	static void conservationCheck(double *previous,double *current)
	{
		
	}
	static void linearSystemUpdate(int i,int j,InitPoint2D *point)
	{
		printf("Updating(%d,%d)...\n",i,j);
	}


	static void init(int i,int j,InitPoint2D *point)
	{
		int nx = globals[0];
		int ny = globals[1];
		if(i == nx/2 && j == ny/2)
			point->U_input[0] = 1;
		else if(i == (nx/2)+1 && j == (ny/2)+1)
			point->U_input[0] = 1;
		else if(i == (nx/2) && j == (ny/2)+1)
			point->U_input[0] = 1;
		else if(i == (nx/2)+1 && j == (ny/2))
			point->U_input[0] = 1;
		else
			point->U_input[0] = 0;
		
	}

	static void setGlobalVariables(int nx,int ny)
	{
		
		globals[0] = nx;
		globals[1] = ny;
		if(pg.rank == 0)
			printf("Global Variables Setup DONE!\n");
	}

	static void setInitParameters(int nx,int ny)
	{
		//PARAMETERS SETUP
		totalGlobals = 2;
		substeps = 1;
		dataPointSize = 2;
	}	
};

#endif
