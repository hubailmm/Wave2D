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

#include "SweptDiscretization2D.h"

SweptDiscretization2D::SweptDiscretization2D(int n,int substeps,int dataPointSize,int totalConstants,int totalConservedQuantities,int remoteConstantsCount)
{
	if(totalConstants<remoteConstantsCount)
	{
		printf("Remote Constants Cannot be more than the total constants!!\n");
		exit(1);
	}
	this->n = n;
	this->substeps = substeps;
	this->dataPointSize = dataPointSize;
	//this->constants = totalConstants;
	this->constants = totalConstants + totalConservedQuantities;
	this->remoteConstantsCount = remoteConstantsCount;
	this->totalConservedQuantities = totalConservedQuantities;
	this->checkForRemoteCycles = -1;
	this->outputToFile = -1;
	this->conserveCheckFreq = -1;
	this->firstConservedCheck = true;
	this->outputLength = 1;
	panelSize = 2;
	for(int i=n;i>2;i=i-2)panelSize += i;panelSize *= 2;
	int constantsToAdd = panelSize * this->constants;
	panelSize = panelSize * substeps * dataPointSize;
	communicationSize = panelSize + constantsToAdd;
	foundationSize = ((n+2) * (n+2) * substeps * dataPointSize + ((n+2) * (n+2) * this->constants));
	resultArray = NULL;		
	
	MPI_Alloc_mem(foundationSize * sizeof(double), MPI_INFO_NULL, &foundation);
	MPI_Win_create(foundation,foundationSize * sizeof(double),1, MPI_INFO_NULL,MPI_COMM_WORLD, &foundationWindow);
	if(remoteConstantsCount > 0)
	{
		MPI_Alloc_mem(n * n * remoteConstantsCount * sizeof(double), MPI_INFO_NULL, &remoteConstants);
		MPI_Win_create(remoteConstants,n * n * remoteConstantsCount * sizeof(double),1, MPI_INFO_NULL,MPI_COMM_WORLD, &constantsWindow);
		constantsArrayBytes = (unsigned char*) malloc(n * n * remoteConstantsCount * sizeof(double) * pg.mpiSize);		
	}

	if(totalConservedQuantities > 0)
	{
		convervedVariables = new double[totalConservedQuantities];
		memset(convervedVariables,'0',sizeof(double) * totalConservedQuantities);
	}
	staging    = new double[foundationSize];
	northPanel = new double[panelSize + constantsToAdd];
	southPanel = new double[panelSize + constantsToAdd];
	eastPanel  = new double[panelSize + constantsToAdd];
	westPanel  = new double[panelSize + constantsToAdd];

	//Define Swept2D components
	up = new UpPyramid(n,foundation,staging,substeps,dataPointSize,this->constants,totalConservedQuantities);
	dp = new DownPyramid(n,foundation,staging,substeps,dataPointSize,this->constants,totalConservedQuantities);
	hb = new HorizontalBridge(n,foundation,staging,substeps,dataPointSize,this->constants,totalConservedQuantities);
	vb = new VerticalBridge(n,foundation,staging,substeps,dataPointSize,this->constants,totalConservedQuantities);

	firstRun = true;
}

int SweptDiscretization2D::ijToIndex(int i,int j)
{
	return Swept2DUtils::ijToIndex(n+2,i,j,this->substeps,this->dataPointSize,this->constants);
}

int SweptDiscretization2D::ijToConstantIndex(int i,int j)
{
	return this->remoteConstantsCount * (i + j * n);	
}

void SweptDiscretization2D::setConservationCheck(int cycles)
{
	this->conserveCheckFreq = cycles;
}

void SweptDiscretization2D::setOutputLength(int outputLength)
{
	this->outputLength = outputLength;
}

void SweptDiscretization2D::reportConveredQuantities()
{
	double *newConservedLocal = new double[this->totalConservedQuantities];
	double *newConservedAll   = new double[this->totalConservedQuantities];
	for(int c=0;c<this->totalConservedQuantities;c++)newConservedLocal[c] = 0;
	for(int i=1;i<n+1;i++)
	{
		for(int j=1;j<n+1;j++)
		{
			for(int c=0;c<this->totalConservedQuantities;c++)
			{
				int cellIndex = ijToIndex(i,j);
				int shift = this->constants-this->totalConservedQuantities + c;
				double value = this->foundation[cellIndex + shift];
				newConservedLocal[c] += value;
			}
		}		
	}
	MPI_Allreduce(newConservedLocal,newConservedAll,this->totalConservedQuantities,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	/*
	printf("\n\n");
	if(this->firstConservedCheck)
	{
		for(int c=0;c<this->totalConservedQuantities;c++)
			printf("Conserved Quantity %d: = %f\n",c+1,newConservedAll[c]);			
		this->firstConservedCheck = false;
	}
	else
	{
		for(int c=0;c<this->totalConservedQuantities;c++)
			printf("Conserved Quantity %d: OLD = %f , NEW = %f\n",c+1,convervedVariables[c],newConservedAll[c]);
	}
	printf("\n");
	*/
	Substeps2D::conservationCheck(this->convervedVariables,newConservedAll);
	memcpy(this->convervedVariables,newConservedAll,sizeof(double) * this->totalConservedQuantities);
	delete[] newConservedLocal;
	delete[] newConservedAll;
}
void SweptDiscretization2D::updateRemoteConstants(unsigned char *buffer)
{
	void *sendingBuffer = NULL;
	FILE *inFile = NULL;

	if(pg.rank == 0)
	{
		int bufferSize = this->remoteConstantsCount * n * n * pg.mpiSize * sizeof(double);
		MPI_Alloc_mem(bufferSize, MPI_INFO_NULL, &sendingBuffer);
		
		for(int r=0;r<pg.mpiSize;r++)
		{
			double *processing = (double*)sendingBuffer + (this->remoteConstantsCount * n * n * r);
			int jIndex = (r % (pg.xNodes*pg.yNodes)) / pg.xNodes;
			int iIndex = r % pg.xNodes;
			for(int j=0;j<n;j++)
			{
				for(int i=0;i<n;i++)
				{
					int iGlobal = n*iIndex + (i);
					int jGlobal = n*jIndex + (j);
					int index   = this->ijToConstantIndex(i,j);
					int globalIndex = this->remoteConstantsCount * (iGlobal + jGlobal * n * pg.xNodes);
					for(int k=0;k<this->remoteConstantsCount;k++)
					{
						processing[index + k] = ((double*)buffer)[k + globalIndex];
					}					
				}
			}		
		}
	}
	
	MPI_Win_fence(MPI_MODE_NOPRECEDE, this->constantsWindow);
	if(pg.rank == 0)
	{
		for(int r=0;r<pg.mpiSize;r++)
		{
			MPI_Put((unsigned char*)sendingBuffer + (r * remoteConstantsCount * n * n * sizeof(double)), remoteConstantsCount * n * n * sizeof(double), MPI_BYTE, r, 0, remoteConstantsCount * n * n * sizeof(double), MPI_BYTE, constantsWindow);			
		}
	}
	MPI_Win_fence((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), this->constantsWindow);

	if(pg.rank == 0)
	{
		MPI_Free_mem(sendingBuffer);					
	}
	for(int i=1;i<n+1;i++)
	{
		for(int j=1;j<n+1;j++)
		{
			for(int k=0;k<this->remoteConstantsCount;k++)
			{
				int windowIndex    = this->ijToConstantIndex(i-1,j-1);
				int foundationIndex = this->ijToIndex(i,j);
				this->foundation[foundationIndex + k] = this->remoteConstants[windowIndex + k];
			}
		}
	}
}
void SweptDiscretization2D::allGatherOutputToFile(int dataPoint,string filename)
{

	void *buffer = NULL;
	FILE *output;
	if(pg.rank == 0)
	{
		MPI_Alloc_mem(foundationSize * pg.mpiSize * sizeof(double), MPI_INFO_NULL, &buffer);		
		output = fopen(filename.c_str(),"wb");
	}

	MPI_Win_fence((MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE), foundationWindow);
	if(pg.rank == 0)
	{
		for(int r=0;r<pg.mpiSize;r++)
		{
			MPI_Get((char*)buffer + (r * foundationSize * sizeof(double)), foundationSize * sizeof(double), MPI_BYTE, r, 0, foundationSize * sizeof(double), MPI_BYTE, foundationWindow);			
		}
	}
	MPI_Win_fence(MPI_MODE_NOSUCCEED, foundationWindow);
	
	if(pg.rank == 0)
	{
		int w = (n * pg.xNodes);
		int h = (n * pg.yNodes);
		int resultArraySize = w * h;
		if(resultArray == NULL)
			resultArray = (double*) malloc(resultArraySize * sizeof(double));
		
		for(int r=0;r<pg.mpiSize;r++)
		{
			double *processing = (double*)buffer + (foundationSize * r);			
			int jIndex = (r % (pg.xNodes*pg.yNodes)) / pg.xNodes;
			int iIndex = r % pg.xNodes;
			for(int j=1;j<n+1;j++)
			{
				for(int i=1;i<n+1;i++)
				{
					int iGlobal = n*iIndex + (i-1);
					int jGlobal = n*jIndex + (j-1);
					int index   = this->ijToIndex(i,j);
					double val  = processing[index + constants + dataPoint];
					int resultIndex = (iGlobal + jGlobal * n * pg.xNodes);
					resultArray[resultIndex] = val;
					//printf("Result Index = %d -- Global[%d][%d] = %f\n",resultIndex,iGlobal,jGlobal,val);					
				}
			}		
		}
		fwrite((const void*)resultArray,sizeof(double),resultArraySize,output);
		fclose(output);
		MPI_Free_mem(buffer);
				
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void SweptDiscretization2D::allGatherAllOutputToFile(string filename)
{

	void *buffer = NULL;
	FILE *output;
	if(pg.rank == 0)
	{
		MPI_Alloc_mem(foundationSize * pg.mpiSize * sizeof(double), MPI_INFO_NULL, &buffer);		
		output = fopen(filename.c_str(),"wb");
	}

	MPI_Win_fence((MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE), foundationWindow);
	if(pg.rank == 0)
	{
		for(int r=0;r<pg.mpiSize;r++)
		{
			MPI_Get((char*)buffer + (r * foundationSize * sizeof(double)), foundationSize * sizeof(double), MPI_BYTE, r, 0, foundationSize * sizeof(double), MPI_BYTE, foundationWindow);			
		}
	}
	MPI_Win_fence(MPI_MODE_NOSUCCEED, foundationWindow);
	
	if(pg.rank == 0)
	{
		int w = (n * pg.xNodes);
		int h = (n * pg.yNodes);
		int resultArraySize = w * h;
		if(resultArray == NULL)
			resultArray = (double*) malloc(resultArraySize * sizeof(double) * outputLength);
		
		for(int r=0;r<pg.mpiSize;r++)
		{
			double *processing = (double*)buffer + (foundationSize * r);			
			int jIndex = (r % (pg.xNodes*pg.yNodes)) / pg.xNodes;
			int iIndex = r % pg.xNodes;
			for(int j=1;j<n+1;j++)
			{
				for(int i=1;i<n+1;i++)
				{
					int iGlobal = n*iIndex + (i-1);
					int jGlobal = n*jIndex + (j-1);
					int index   = this->ijToIndex(i,j);
					for(int point=0;point<outputLength;point++)
					{
						double val  = processing[index + constants + point];
						int resultIndex = (iGlobal + jGlobal * n * pg.xNodes) * outputLength + point;
						resultArray[resultIndex] = val;
					}					
				}
			}		
		}
		fwrite((const void*)resultArray,sizeof(double),resultArraySize,output);
		fclose(output);
		MPI_Free_mem(buffer);
				
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


void SweptDiscretization2D::setInitCondition(void initFnc(int,int,InitPoint2D*))
{
	up->setInitCondition(initFnc);	
	//up->printFoundation();
}

void SweptDiscretization2D::updateSystem(void updateFnc(int,int,InitPoint2D*))
{
	up->updateSystem(updateFnc);	
	//up->printFoundation();
}
void SweptDiscretization2D::setOutputDirectory(string outputDirectory)
{
	this->outputDirectory.clear();
	this->outputDirectory += outputDirectory;		
}

/*
void SweptDiscretization2D::printFoundation(int dataPoint)
{
	dp->printFoundation(foundation,dataPoint);		
}
*/

void SweptDiscretization2D::printConstantFromWindow(int constNum)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			int index  = this->ijToConstantIndex(i,j);
			double value = this->remoteConstants[index + constNum];
			printf("%.1f ",value);
		}
		printf("\n");
	}
}

void SweptDiscretization2D::printConstantFromFoundation(int constNum)
{
	for(int i=1;i<n+1;i++)
	{
		for(int j=1;j<n+1;j++)
		{
			int index  = this->ijToIndex(i,j);
			double value = this->foundation[index + constNum];
			printf("%.1f ",value);
		}
		printf("\n");
	}
}

void SweptDiscretization2D::getFoundation(double *result,int dataPoint)
{
	dp->getFoundation(foundation,result,dataPoint);
}

void SweptDiscretization2D::setCheckForRemoteCommands(int cycles)
{
	this->checkForRemoteCycles = cycles;
}

void SweptDiscretization2D::setFileOutput(int cycles)
{
	this->outputToFile = cycles;
}

double SweptDiscretization2D::calculate(int cycles)
{
	double outputTime = 0;
	if(!firstRun)
		this->updateSystem(Substeps2D::linearSystemUpdate);
	else
	{
		firstRun = false;
		/*
		for(int i=0;i<this->dataPointSize;i++)
		{
			char filename[80];
			memset(filename,'\0',80);
			sprintf(filename,"output%d_0.bin",i);		
			string file(outputDirectory + "/" + filename);
			printf("Generating Output File: %s\n",file.c_str());
			this->allGatherOutputToFile(i,file);
		}
		*/
	}
	
	
	//printf("\n");
	this->reportConveredQuantities();
	int execFnc = 0;
	double startTime = MPI_Wtime();
	for(int c=1;c<=cycles;c++)
	{
		if(c % 10000 == 0 && pg.rank == 0)	
			printf("Performing Swept Cycle: %d\n",c);
		//Sleep(100);
		//Build the first Top Pyramid	
		up->build(execFnc);
		MPI_Isend(up->getWestPanel() ,communicationSize,MPI_DOUBLE,pg.Nrank ,2,MPI_COMM_WORLD,&reqs[1]);		
		MPI_Isend(up->getNorthPanel(),communicationSize,MPI_DOUBLE,pg.Wrank,1,MPI_COMM_WORLD,&reqs[0]);
		MPI_Irecv(eastPanel,communicationSize,MPI_DOUBLE,pg.Srank,2,MPI_COMM_WORLD,&reqs[2]);	
		MPI_Irecv(southPanel,communicationSize,MPI_DOUBLE,pg.Erank,1,MPI_COMM_WORLD,&reqs[3]);
		
		//Build the fist pair of Bridges
		MPI_Wait(&reqs[2],MPI_STATUS_IGNORE);		
		hb->setEastWestPanels(eastPanel,up->getEastPanel());
		hb->build(execFnc);
		MPI_Isend(hb->getNorthPanel(),communicationSize,MPI_DOUBLE,pg.Wrank ,1,MPI_COMM_WORLD,&reqs[2]);

		MPI_Wait(&reqs[3],MPI_STATUS_IGNORE);
		vb->setNorthSouthPanels(up->getSouthPanel(),southPanel);
		vb->build(execFnc);
		MPI_Isend(vb->getWestPanel(),communicationSize,MPI_DOUBLE,pg.Nrank,2,MPI_COMM_WORLD,&reqs[3]);
		
		MPI_Wait(&reqs[0],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[1],MPI_STATUS_IGNORE);

		MPI_Irecv(southPanel,communicationSize,MPI_DOUBLE,pg.Erank,1,MPI_COMM_WORLD,&reqs[0]);
		MPI_Irecv(eastPanel,communicationSize,MPI_DOUBLE,pg.Srank,2,MPI_COMM_WORLD,&reqs[1]);

		MPI_Wait(&reqs[0],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[1],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[2],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[3],MPI_STATUS_IGNORE);
		
		//Build the firse Down Pyramid
		dp->setPanels(hb->getSouthPanel(),southPanel,eastPanel,vb->getEastPanel());
		dp->build(execFnc);
		execFnc += (n)/2;
		execFnc = execFnc % substeps;
		
		//Build the second Up Pyramid
		up->build(execFnc);
		MPI_Isend(up->getEastPanel() ,communicationSize,MPI_DOUBLE,pg.Srank ,1,MPI_COMM_WORLD,&reqs[0]);
		MPI_Isend(up->getSouthPanel(),communicationSize,MPI_DOUBLE,pg.Erank,2,MPI_COMM_WORLD,&reqs[1]);
		MPI_Irecv(eastPanel,communicationSize,MPI_DOUBLE,pg.Nrank,1,MPI_COMM_WORLD,&reqs[2]);
		MPI_Irecv(southPanel,communicationSize,MPI_DOUBLE,pg.Wrank,2,MPI_COMM_WORLD,&reqs[3]);
		
		//Build the second pair of Bridges
		MPI_Wait(&reqs[2],MPI_STATUS_IGNORE);
		hb->setEastWestPanels(up->getWestPanel(),eastPanel);
		hb->build(execFnc);
		MPI_Isend(hb->getSouthPanel(),communicationSize,MPI_DOUBLE,pg.Erank,1,MPI_COMM_WORLD,&reqs[2]);

		MPI_Wait(&reqs[3],MPI_STATUS_IGNORE);
		vb->setNorthSouthPanels(southPanel,up->getNorthPanel());
		vb->build(execFnc);
		MPI_Isend(vb->getEastPanel(),communicationSize,MPI_DOUBLE,pg.Srank,2,MPI_COMM_WORLD,&reqs[3]);
		
		MPI_Wait(&reqs[0],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[1],MPI_STATUS_IGNORE);
				
		MPI_Irecv(northPanel,communicationSize,MPI_DOUBLE,pg.Wrank,1,MPI_COMM_WORLD,&reqs[0]);
		MPI_Irecv(westPanel,communicationSize,MPI_DOUBLE,pg.Nrank,2,MPI_COMM_WORLD,&reqs[1]);

		MPI_Wait(&reqs[0],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[1],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[2],MPI_STATUS_IGNORE);
		MPI_Wait(&reqs[3],MPI_STATUS_IGNORE);

		//Build the second Down Pyramid, which completes the cycle
		dp->setPanels(northPanel,hb->getNorthPanel(),vb->getWestPanel(),westPanel);
		dp->build(execFnc);
		execFnc += (n)/2;
		execFnc = execFnc % substeps;

		if((outputToFile != -1) && (c % outputToFile == 0))
		{
			double outputStart = MPI_Wtime();
			for(int i=0;i<this->outputLength;i++)
			{
				char filename[80];
				memset(filename,'\0',80);
				sprintf(filename,"output%d_%d.bin",i,c);				
				string file(outputDirectory + "/" + filename);
				printf("Generating Output File: %s\n",file.c_str());
				this->allGatherOutputToFile(i,file);				
			}
			printf("\n");
			outputTime += MPI_Wtime()-outputStart;
		}

		if((this->conserveCheckFreq != -1) && (c % this->conserveCheckFreq == 0))
		{
			this->reportConveredQuantities();
		}
	}
	//printf("Total Output Time: %f\n",outputTime*1000000);
	return (MPI_Wtime() - startTime) - outputTime;
}

void SweptDiscretization2D::printFullFoundation(int dataPoint)
{
	for(int j=0;j<=n-1+2;j++)
	{
		for(int i=0;i<=n-1+2;i++)
		{
			printf("%.1f ",foundation[ijToIndex(i,j) + constants + dataPoint]);
		}
		printf("\n");
	}
	printf("\n\n");
}


void SweptDiscretization2D::printFoundation(int dataPoint)
{
	for(int j=1;j<n-1+2;j++)
	{
		for(int i=1;i<n-1+2;i++)
		{
			printf("%.8f ",foundation[ijToIndex(i,j) + constants + dataPoint]);
		}
		printf("\n");
	}
	printf("\n\n");
}


void SweptDiscretization2D::printFoundationHex(int dataPoint)
{
	for(int j=1;j<n-1+2;j++)
	{
		for(int i=1;i<n-1+2;i++)
		{
			DebugUtils::printHex(foundation[ijToIndex(i,j) + constants + dataPoint]);
			printf(" ");
		}
		printf("\n");
	}
	printf("\n\n");
}

void SweptDiscretization2D::printFullFoundationHex(int dataPoint)
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

SweptDiscretization2D::~SweptDiscretization2D(void)
{
	delete up;
	delete dp;
	delete hb;
	delete vb;
	delete[] staging;
	delete[] northPanel;
	delete[] southPanel;
	delete[] eastPanel;
	delete[] westPanel;
}
