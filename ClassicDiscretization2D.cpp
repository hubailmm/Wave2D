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

#include "ClassicDiscretization2D.h"

ClassicDiscretization2D::ClassicDiscretization2D(int n,int substeps,int dataPointSize,int constants,int totalConservedQuantities)
{
	this->n = n+2;
	this->substeps = substeps;
	this->dataPointSize = dataPointSize;
	this->totalConservedQuantitie = totalConservedQuantities;
	this->constants = constants + totalConservedQuantities;
	this->communicationSize = n * substeps * dataPointSize + n * this->constants;
	int foundationSize = (n+2) * (n+2) * substeps * dataPointSize + ((n+2) * (n+2) * this->constants);
	foundation  = new double[(n+2) * (n+2) * substeps * dataPointSize + ((n+2) * (n+2) * this->constants)];
	staging     = new double[(n+2) * (n+2) * substeps * dataPointSize + ((n+2) * (n+2) * this->constants)];
	memset(foundation,'0',sizeof(double) * foundationSize);
	memset(staging   ,'0',sizeof(double) * foundationSize);
	northRegion = new double[n * substeps * dataPointSize + n * this->constants];
	southRegion = new double[n * substeps * dataPointSize + n * this->constants];
	eastRegion  = new double[n * substeps * dataPointSize + n * this->constants];
	westRegion  = new double[n * substeps * dataPointSize + n * this->constants];

	northRegionRecv = new double[n * substeps * dataPointSize + n * this->constants];
	southRegionRecv = new double[n * substeps * dataPointSize + n * this->constants];
	eastRegionRecv  = new double[n * substeps * dataPointSize + n * this->constants];
	westRegionRecv  = new double[n * substeps * dataPointSize + n * this->constants];

	firstRun = true;

}

int ClassicDiscretization2D::ijToIndex(int i,int j)
{
	int index = j * n + i;
	int constToAdd = (index * constants);
	index = index * this->substeps * this->dataPointSize;
	index += constToAdd;
	return index;
}

void ClassicDiscretization2D::setInitCondition(void initFnc(int,int,InitPoint2D*))
{
	for(int i=1;i<n-1;i++)
	{
		for(int j=1;j<n-1;j++)
		{
			InitPoint2D point;
			point.U_input = &foundation[this->ijToIndex(i,j) + constants];
			point.U_constants  = &foundation[this->ijToIndex(i,j)];
			point.U_conserved  = &foundation[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
			initFnc((pg.iIndex*(n-2))+i-1,(pg.jIndex*(n-2))+j-1,&point);			
		}
	}
}

void ClassicDiscretization2D::updateSystem(void updateFnc(int,int,InitPoint2D*))
{
	for(int i=1;i<n-1;i++)
	{
		for(int j=1;j<n-1;j++)
		{
			InitPoint2D point;
			point.U_input = &foundation[this->ijToIndex(i,j) + constants];
			point.U_constants  = &foundation[this->ijToIndex(i,j)];
			point.U_conserved  = &foundation[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
			updateFnc((pg.iIndex*(n-2))+i-1,(pg.jIndex*(n-2))+j-1,&point);			
		}
	}
}

void ClassicDiscretization2D::sendGhostRegion(double *dataSource)
{
	//West Region
	int index = 0;	
	for(int j=1;j<n-1;j++)
	{
		memcpy(&this->westRegion[index],&dataSource[this->ijToIndex(1,j)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->westRegion,communicationSize,MPI_DOUBLE,pg.Wrank,1,MPI_COMM_WORLD,&reqs[0]);

	//East Region
	index = 0;
	for(int j=1;j<n-1;j++)
	{
		memcpy(&this->eastRegion[index],&dataSource[this->ijToIndex(n-2,j)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->eastRegion,communicationSize,MPI_DOUBLE,pg.Erank,2,MPI_COMM_WORLD,&reqs[1]);

	//South Region
	index = 0;
	for(int i=1;i<n-1;i++)
	{
		memcpy(&this->southRegion[index],&dataSource[this->ijToIndex(i,n-2)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->southRegion,communicationSize,MPI_DOUBLE,pg.Srank,3,MPI_COMM_WORLD,&reqs[2]);

	//North Region
	index = 0;
	for(int i=1;i<n-1;i++)
	{
		memcpy(&this->northRegion[index],&dataSource[this->ijToIndex(i,1)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->northRegion,communicationSize,MPI_DOUBLE,pg.Nrank,4,MPI_COMM_WORLD,&reqs[3]);	
}

void ClassicDiscretization2D::recvGhostRegion(double *destData)
{
	//EastRegion
	int index = 0;
	MPI_Recv(this->eastRegionRecv,communicationSize,MPI_DOUBLE,pg.Erank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int j=1;j<n-1;j++)
	{
		memcpy(&destData[this->ijToIndex(n-1,j)],&this->eastRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	//WestRegion
	index = 0;
	MPI_Recv(this->westRegionRecv,communicationSize,MPI_DOUBLE,pg.Wrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int j=1;j<n-1;j++)
	{
		memcpy(&destData[this->ijToIndex(0,j)],&this->westRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	//NorthRegion
	index = 0;
	MPI_Recv(this->northRegionRecv,communicationSize,MPI_DOUBLE,pg.Nrank,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int i=1;i<n-1;i++)
	{
		memcpy(&destData[this->ijToIndex(i,0)],&this->northRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	//SouthRegion
	index = 0;	
	MPI_Recv(this->southRegionRecv,communicationSize,MPI_DOUBLE,pg.Srank,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int i=1;i<n-1;i++)
	{
		memcpy(&destData[this->ijToIndex(i,n-1)],&this->southRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	MPI_Request_free(&reqs[0]);
	MPI_Request_free(&reqs[1]);
	MPI_Request_free(&reqs[2]);
	MPI_Request_free(&reqs[3]);
}
void ClassicDiscretization2D::initGhostRegion(double *dataSource,double *destData)
{
	//West Region
	int index = 0;	
	for(int j=1;j<n-1;j++)
	{
		memcpy(&this->westRegion[index],&dataSource[this->ijToIndex(1,j)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->westRegion,communicationSize,MPI_DOUBLE,pg.Wrank,1,MPI_COMM_WORLD,&reqs[0]);

	//East Region
	index = 0;
	for(int j=1;j<n-1;j++)
	{
		memcpy(&this->eastRegion[index],&dataSource[this->ijToIndex(n-2,j)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->eastRegion,communicationSize,MPI_DOUBLE,pg.Erank,2,MPI_COMM_WORLD,&reqs[1]);

	//South Region
	index = 0;
	for(int i=1;i<n-1;i++)
	{
		memcpy(&this->southRegion[index],&dataSource[this->ijToIndex(i,n-2)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->southRegion,communicationSize,MPI_DOUBLE,pg.Srank,3,MPI_COMM_WORLD,&reqs[2]);

	//North Region
	index = 0;
	for(int i=1;i<n-1;i++)
	{
		memcpy(&this->northRegion[index],&dataSource[this->ijToIndex(i,1)],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}
	MPI_Isend(this->northRegion,communicationSize,MPI_DOUBLE,pg.Nrank,4,MPI_COMM_WORLD,&reqs[3]);

	//EastRegion
	index = 0;
	MPI_Recv(this->eastRegionRecv,communicationSize,MPI_DOUBLE,pg.Erank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int j=1;j<n-1;j++)
	{
		memcpy(&destData[this->ijToIndex(n-1,j)],&this->eastRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	//WestRegion
	index = 0;
	MPI_Recv(this->westRegionRecv,communicationSize,MPI_DOUBLE,pg.Wrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int j=1;j<n-1;j++)
	{
		memcpy(&destData[this->ijToIndex(0,j)],&this->westRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	//NorthRegion
	index = 0;
	MPI_Recv(this->northRegionRecv,communicationSize,MPI_DOUBLE,pg.Nrank,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int i=1;i<n-1;i++)
	{
		memcpy(&destData[this->ijToIndex(i,0)],&this->northRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	//SouthRegion
	index = 0;	
	MPI_Recv(this->southRegionRecv,communicationSize,MPI_DOUBLE,pg.Srank,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	for(int i=1;i<n-1;i++)
	{
		memcpy(&destData[this->ijToIndex(i,n-1)],&this->southRegionRecv[index],((substeps * dataPointSize) + constants) * sizeof(double));		
		index += ((substeps * dataPointSize) + constants);
	}

	MPI_Request_free(&reqs[0]);
	MPI_Request_free(&reqs[1]);
	MPI_Request_free(&reqs[2]);
	MPI_Request_free(&reqs[3]);
}

double ClassicDiscretization2D::calculate(int steps)
{
	if(!firstRun)
		this->updateSystem(Substeps2D::linearSystemUpdate);

	firstRun = false;

	double *sourceData;
	double *destData;
	int execute = 0;
	this->initGhostRegion(foundation,foundation);
	double startTime = MPI_Wtime();
	for(int s=1;s<=steps;s++)
	{
		if(s%2 != 0)
		{
			destData   = staging;
			sourceData = foundation;
		}
		else
		{
			destData   = foundation;
			sourceData = staging;
		}
				
		//West Region
		for(int j=1;j<n-2;j++)
		{
			int i = 1;
			PointStruct2D point;
			point.ConstantsOutput = &destData[this->ijToIndex(i,j)];
			memcpy(&destData[this->ijToIndex(i,j)],&sourceData[this->ijToIndex(i,j)],constants * sizeof(double));
			point.output = &destData[this->ijToIndex(i,j) + constants];
			point.C_input = &sourceData[this->ijToIndex(i,j) + constants];
			point.C_constants  = &sourceData[this->ijToIndex(i,j)];
			point.C_conserved  = &destData[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
			point.W_input = &sourceData[this->ijToIndex(i-1,j) + constants];
			point.W_constants  = &sourceData[this->ijToIndex(i-1,j)];
			point.E_input = &sourceData[this->ijToIndex(i+1,j) + constants];
			point.E_constants  = &sourceData[this->ijToIndex(i+1,j)];
			point.S_input = &sourceData[this->ijToIndex(i,j+1) + constants];
			point.S_constants  = &sourceData[this->ijToIndex(i,j+1)];
			point.N_input = &sourceData[this->ijToIndex(i,j-1) + constants];
			point.N_constants  = &sourceData[this->ijToIndex(i,j-1)];				
			Substeps2D::executeStepFnc(execute,&point);
		}
		
		//East Region
		for(int j=2;j<n-1;j++)
		{
			int i = n-2;
			PointStruct2D point;
			point.ConstantsOutput = &destData[this->ijToIndex(i,j)];
			memcpy(&destData[this->ijToIndex(i,j)],&sourceData[this->ijToIndex(i,j)],constants * sizeof(double));
			point.output = &destData[this->ijToIndex(i,j) + constants];
			point.C_input = &sourceData[this->ijToIndex(i,j) + constants];
			point.C_constants  = &sourceData[this->ijToIndex(i,j)];
			point.C_conserved  = &destData[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
			point.W_input = &sourceData[this->ijToIndex(i-1,j) + constants];
			point.W_constants  = &sourceData[this->ijToIndex(i-1,j)];
			point.E_input = &sourceData[this->ijToIndex(i+1,j) + constants];
			point.E_constants  = &sourceData[this->ijToIndex(i+1,j)];
			point.S_input = &sourceData[this->ijToIndex(i,j+1) + constants];
			point.S_constants  = &sourceData[this->ijToIndex(i,j+1)];
			point.N_input = &sourceData[this->ijToIndex(i,j-1) + constants];
			point.N_constants  = &sourceData[this->ijToIndex(i,j-1)];
			Substeps2D::executeStepFnc(execute,&point);								
		}
		
		//South Region
		for(int i=1;i<n-2;i++)
		{
			int j = n-2;
			PointStruct2D point;
			point.ConstantsOutput = &destData[this->ijToIndex(i,j)];
			memcpy(&destData[this->ijToIndex(i,j)],&sourceData[this->ijToIndex(i,j)],constants * sizeof(double));
			point.output = &destData[this->ijToIndex(i,j) + constants];
			point.C_input = &sourceData[this->ijToIndex(i,j) + constants];
			point.C_constants  = &sourceData[this->ijToIndex(i,j)];
			point.C_conserved  = &destData[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
			point.W_input = &sourceData[this->ijToIndex(i-1,j) + constants];
			point.W_constants  = &sourceData[this->ijToIndex(i-1,j)];
			point.E_input = &sourceData[this->ijToIndex(i+1,j) + constants];
			point.E_constants  = &sourceData[this->ijToIndex(i+1,j)];
			point.S_input = &sourceData[this->ijToIndex(i,j+1) + constants];
			point.S_constants  = &sourceData[this->ijToIndex(i,j+1)];
			point.N_input = &sourceData[this->ijToIndex(i,j-1) + constants];
			point.N_constants  = &sourceData[this->ijToIndex(i,j-1)];
			Substeps2D::executeStepFnc(execute,&point);							
		}
		
		//North Region
		for(int i=2;i<n-1;i++)
		{
			int j = 1;
			PointStruct2D point;
			point.ConstantsOutput = &destData[this->ijToIndex(i,j)];
			memcpy(&destData[this->ijToIndex(i,j)],&sourceData[this->ijToIndex(i,j)],constants * sizeof(double));
			point.output = &destData[this->ijToIndex(i,j) + constants];
			point.C_input = &sourceData[this->ijToIndex(i,j) + constants];
			point.C_constants  = &sourceData[this->ijToIndex(i,j)];
			point.C_conserved  = &destData[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
			point.W_input = &sourceData[this->ijToIndex(i-1,j) + constants];
			point.W_constants  = &sourceData[this->ijToIndex(i-1,j)];
			point.E_input = &sourceData[this->ijToIndex(i+1,j) + constants];
			point.E_constants  = &sourceData[this->ijToIndex(i+1,j)];
			point.S_input = &sourceData[this->ijToIndex(i,j+1) + constants];
			point.S_constants  = &sourceData[this->ijToIndex(i,j+1)];
			point.N_input = &sourceData[this->ijToIndex(i,j-1) + constants];
			point.N_constants  = &sourceData[this->ijToIndex(i,j-1)];
			Substeps2D::executeStepFnc(execute,&point);				
		}		
		this->sendGhostRegion(destData);
		//
		for(int i=2;i<(n-2);i++)
		{
			for(int j=2;j<(n-2);j++)
			{
				PointStruct2D point;
				point.ConstantsOutput = &destData[this->ijToIndex(i,j)];
				memcpy(&destData[this->ijToIndex(i,j)],&sourceData[this->ijToIndex(i,j)],constants * sizeof(double));
				point.output = &destData[this->ijToIndex(i,j) + constants];
				point.C_input = &sourceData[this->ijToIndex(i,j) + constants];
				point.C_constants  = &sourceData[this->ijToIndex(i,j)];
				point.C_conserved  = &destData[(this->ijToIndex(i,j)) + (this->constants-this->totalConservedQuantitie)];
				point.W_input = &sourceData[this->ijToIndex(i-1,j) + constants];
				point.W_constants  = &sourceData[this->ijToIndex(i-1,j)];
				point.E_input = &sourceData[this->ijToIndex(i+1,j) + constants];
				point.E_constants  = &sourceData[this->ijToIndex(i+1,j)];
				point.S_input = &sourceData[this->ijToIndex(i,j+1) + constants];
				point.S_constants  = &sourceData[this->ijToIndex(i,j+1)];
				point.N_input = &sourceData[this->ijToIndex(i,j-1) + constants];
				point.N_constants  = &sourceData[this->ijToIndex(i,j-1)];				
				Substeps2D::executeStepFnc(execute,&point);				
			}
		}
		this->recvGhostRegion(destData);
		execute++;
		if(execute==substeps)execute = 0;
	}
	if(this->foundation != destData)
		memcpy(foundation,destData,sizeof(double) * ((n*n) * substeps * dataPointSize + ((n*n) * constants)));	
	return MPI_Wtime() - startTime;
}

void ClassicDiscretization2D::printFullFoundation(int dataPoint)
{
	for(int j=0;j<=n-1;j++)
	{
		for(int i=0;i<=n-1;i++)
		{
			printf("%.1f ",foundation[ijToIndex(i,j) + constants + dataPoint]);
		}
		printf("\n");
	}
	printf("\n\n");
}

void ClassicDiscretization2D::printFoundation(int dataPoint)
{
	for(int j=1;j<n-1;j++)
	{
		for(int i=1;i<n-1;i++)
		{
			printf("%.8f ",foundation[ijToIndex(i,j) + constants + dataPoint]);
		}
		printf("\n");
	}
	printf("\n\n");
}

void ClassicDiscretization2D::printFoundationHex(int dataPoint)
{
	for(int j=1;j<n-1;j++)
	{
		for(int i=1;i<n-1;i++)
		{
			DebugUtils::printHex(foundation[ijToIndex(i,j) + constants + dataPoint]);
			printf(" ");
		}
		printf("\n");
	}
	printf("\n\n");
}

void ClassicDiscretization2D::printFullFoundationHex(int dataPoint)
{
	for(int j=0;j<=n-1;j++)
	{
		for(int i=0;i<=n-1;i++)
		{
			DebugUtils::printHex(foundation[ijToIndex(i,j) + constants + dataPoint]);
			printf(" ");
		}
		printf("\n");
	}
	printf("\n\n");
}
ClassicDiscretization2D::~ClassicDiscretization2D()
{
	delete[] foundation;
	delete[] staging;
	delete[] northRegion;
	delete[] southRegion;
	delete[] eastRegion;
	delete[] westRegion;
	delete[] northRegionRecv;
	delete[] southRegionRecv;
	delete[] eastRegionRecv;
	delete[] westRegionRecv;
}
