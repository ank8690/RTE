/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Today is 14 Jaunary 2018. Sir asked me to write the code for the solution of radiative transfer equation from Scratch.
I need to complete it in 10 days. Let's see! I need to plot the graph between heat flux and length of the square domain of unit length!!

\*---------------------------------------------------------------------------*/

#include "radiativeTEqn.H"
//#include"cell.H"
//#include"face.H"
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>  // for atoi
using namespace std;
float pi= 3.14159f;
float DEGS_TO_RAD = 3.14159f/180.0f;
// first input is absorption coefficient
	double dkappa = 0.1;
double dSBConst = 5.67e-08; // Stefan-Boltzmann Constant
int numVertices = 0;    // Tallies the number of vertex points added.
float fTemp = 1000.0;
int ni = 3;
int nj = 3;
float fEmiss_Bound = 0.1;
// I am trying to construct the grid the way Chai has made using Fortran for RAT
	// X is the cell-centre value at the location
	// XU is the face-centre value  
	//float X[ni],Y[nj],XU[ni],YV[nj], XCV[ni], YCV[nj], Vol[ni][nj];




void radiativeTransferEquation::printFaceAndCellInformation(Face *pFaces, Cell *pCells,long l_nbCells,
															 long l_nbFaces)
{
	// first check for inputs
	if(!l_nbCells || !l_nbFaces || !pFaces || !pCells)
	{
		return;
	}

	//Now, print face information
	long i = 0;
	long j = 0;
	printf("Face Information:\n");
	for(i = 0; i < l_nbFaces; ++i)
	{
		printf("FaceIndex: %d, Number of rays = %d, Temp = %lf, Area = %lf\n",i,pFaces[i].m_l_nbRays, 																					pFaces[i].m_f_Temp,
																				pFaces[i].m_f_Area);
	} 

	printf("Cell Information:\n");
	for(i = 0; i < l_nbCells; ++i)
	{
		printf("CellIndex: %d, Number of rays = %d, Temp = %lf, Volume = %lf\n",i,pCells[i].m_l_nbRays, 																					pCells[i].m_f_Temp,
																				pCells[i].m_f_Vol);

		// Face Indicies to all the sides
		printf("EastFaceIndex: %d, WestFaceIndex = %d, NorthFaceIndex = %d, SouthFaceIndex = %d\n",
																			pCells[i].m_l_FaceIndex_East,
																			pCells[i].m_l_FaceIndex_West, 																				pCells[i].m_l_FaceIndex_North,
																			pCells[i].m_l_FaceIndex_South);
		// Cell Indicies to all the sides
		printf("EastCellIndex: %d, WestCellIndex = %d, NorthCellIndex = %d, SouthCellIndex = %d\n",
																		pCells[i].m_l_CellIndex_East,
																		pCells[i].m_l_CellIndex_West, 																			pCells[i].m_l_CellIndex_North,
																		pCells[i].m_l_CellIndex_South);	
		printf("\n");
	} 
}

// This function is written to accomodate face and cell values seperately
void radiativeTransferEquation::updateIntensity()
{
	int inbRays = 0;
	long i = 0;
	extractDataAndGetRayParameters(inbRays);
	if(!inbRays)
	{
		return;
	}
	Face *pFaces = NULL;
	Cell *pCells = NULL;
	long l_nbCells = 0;
	long l_nbFaces = 0;
	constructGridUsingCellAndFaces(inbRays, pFaces,pCells, l_nbCells, l_nbFaces);
	printFaceAndCellInformation(pFaces,pCells,l_nbCells, l_nbFaces);
	Point pnfe,pnfw,pnfn,pnfs;
	pnfe.x = 1.0;
	pnfe.y = 0.0;
	pnfe.z = 0.0;
	
	pnfw.x = -1.0;
	pnfw.y = 0.0;
	pnfw.z = 0.0;

	pnfn.x = 0.0;
	pnfn.y = 1.0;
	pnfn.z = 0.0;

	pnfs.x = 0.0;
	pnfs.y = -1.0;
	pnfs.z = 0.0;

	float fConvergence = 1e-06; 
	bool *pbRayConverged = NULL;
	pbRayConverged = new bool[inbRays];
	int iRayCt = 0;
	for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
	{
		pbRayConverged[iRayCt] = false;
	}
	int iConvergedRayCt = 0;
	int iNbIterations = 0;
	float f_SMALL = 1e-10;
	float f_DMAX = -1000;
	do
	{
		iNbIterations++;
		iConvergedRayCt = 0;
		f_DMAX = -1000;
		float f_Sum_Old = 0.0;
		float f_Sum_New = 0.0;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(pbRayConverged[iRayCt])
			{
				iConvergedRayCt++;
				continue;
			}
			float fOmega = vec_omega_[iRayCt];
			Point p_cent = vec_Cent_[iRayCt];
			Point p_dAve = vec_Dfm_[iRayCt];
			float fVol = 1.0;
			float fAbsCoeff = 1.0;
			float fdp_e = 0.0;
			float fdp_w = 0.0;
			float fdp_n = 0.0;
			float fdp_s = 0.0;
			dotProduct(p_dAve,pnfe,fdp_e);
			dotProduct(p_dAve,pnfw,fdp_w);
			dotProduct(p_dAve,pnfn,fdp_n);
			dotProduct(p_dAve,pnfs,fdp_s);
			
			// Loop for each cell, where I will write the equation 
			for(int iCellCt = 0; iCellCt < l_nbCells; ++iCellCt)
			{
				long lFaceIndex_East = pCells[iCellCt].m_l_FaceIndex_East;
				long lFaceIndex_West = pCells[iCellCt].m_l_FaceIndex_West;
				long lFaceIndex_North = pCells[iCellCt].m_l_FaceIndex_North;
				long lFaceIndex_South = pCells[iCellCt].m_l_FaceIndex_South;

				long l_Cell_East = pCells[iCellCt].m_l_CellIndex_East;
				long l_Cell_West = pCells[iCellCt].m_l_CellIndex_West;
				long l_Cell_North = pCells[iCellCt].m_l_CellIndex_North;
				long l_Cell_South = pCells[iCellCt].m_l_CellIndex_South;
				float f_Ip_Old = pCells[iCellCt].m_d_vec_cell_ray_Inten[iRayCt];
				// Apply Boundary condition first
					
				// Start from the east face
				if(fdp_e > 0)
				{
					pFaces[lFaceIndex_East].m_d_vec_face_ray_Inten[iRayCt] = f_Ip_Old;
				}
				else
				{
					// otherwise find the cell east to the face of the current cell;
					if(l_Cell_East != -1)
					{
						pFaces[lFaceIndex_East].m_d_vec_face_ray_Inten[iRayCt] = 
													pCells[l_Cell_East].m_d_vec_cell_ray_Inten[iRayCt];
					}
					// Also, we need to use the emissivity boundary condition
					
				}
				// Now , west face
				if(fdp_w > 0)
				{	
					if(l_Cell_West != -1)
					{
						pFaces[lFaceIndex_West].m_d_vec_face_ray_Inten[iRayCt] = 
													pCells[l_Cell_West].m_d_vec_cell_ray_Inten[iRayCt];
					}
				}
				else
				{
					pFaces[lFaceIndex_East].m_d_vec_face_ray_Inten[iRayCt] = f_Ip_Old;
				}
				// North Face
				if(fdp_n > 0)
				{
					pFaces[lFaceIndex_North].m_d_vec_face_ray_Inten[iRayCt] = f_Ip_Old;

				}
				else
				{
					// otherwise find the cell east to the face of the current cell;
					if(l_Cell_North != -1)
					{
						pFaces[lFaceIndex_North].m_d_vec_face_ray_Inten[iRayCt] = 
													pCells[l_Cell_North].m_d_vec_cell_ray_Inten[iRayCt];
					}
				}
				// Now , south face
				if(fdp_s > 0)
				{	
					if(l_Cell_South != -1)
					{
						pFaces[lFaceIndex_South].m_d_vec_face_ray_Inten[iRayCt] = 
												pCells[l_Cell_South].m_d_vec_cell_ray_Inten[iRayCt];
					}
				}
				else
				{
					pFaces[lFaceIndex_South].m_d_vec_face_ray_Inten[iRayCt] = f_Ip_Old;
				}

				// After assigning the intensities on the faces , now we can write the general equation to 							find the intensity on the current cell
				float I_e = pFaces[lFaceIndex_East].m_d_vec_face_ray_Inten[iRayCt];
				float I_w = pFaces[lFaceIndex_West].m_d_vec_face_ray_Inten[iRayCt];
				float I_n = pFaces[lFaceIndex_North].m_d_vec_face_ray_Inten[iRayCt];
				float I_s = pFaces[lFaceIndex_South].m_d_vec_face_ray_Inten[iRayCt];
				float A_e = pFaces[lFaceIndex_East].m_f_Area;
				float A_w = pFaces[lFaceIndex_West].m_f_Area;
				float A_n = pFaces[lFaceIndex_North].m_f_Area;
				float A_s = pFaces[lFaceIndex_South].m_f_Area;
				
				float f_Vol_Cell = pCells[iCellCt].m_f_Vol;
				float fP_V_Omega = f_Vol_Cell*fOmega; 
				float f_Source = fP_V_Omega*fAbsCoeff*dSBConst*pow(pCells[iCellCt].m_f_Temp,4)/pi;
				float f_Num = f_Source -  ( (I_e*A_e*fdp_e) + (I_w*A_w*fdp_w) + (I_n*A_n*fdp_n) + 									(I_s*A_s*fdp_s) );
				float f_DEN = fAbsCoeff*fP_V_Omega;
				
				f_Sum_Old += f_Ip_Old;
				float f_Ip_New = f_Num/(f_DEN + f_SMALL);
				f_Sum_New += f_Ip_New;
				pCells[iCellCt].m_d_vec_cell_ray_Inten[iRayCt] = f_Ip_New;
				float f_diff = fabs(f_Ip_New - f_Ip_Old)/(f_Ip_New + f_SMALL);
				f_DMAX = max(f_DMAX,f_diff);
			}
				/*for(i = 1; i < ni - 1; ++i)
				{
					for(j = 1; j < nj - 1; ++j)
					{
						float fP_V_Omega = ppf_Vol[i-1][j-1]*fOmega; // product of volume and solid angle
						float f_Num = fabs(p_dAve.x*pfYCV[j-1])*I[iRayCt][i-1][j] + 
										fabs(p_dAve.y*pfXCV[i-1])*I[iRayCt][i][j-1] + 											(fP_V_Omega*fAbsCoeff*dSBConst*pow(T[i][j],4)/pi);
							
						float f_DEN = fabs(p_dAve.x*pfYCV[j-1]) + fabs(p_dAve.y*pfXCV[i-1]) + 											fAbsCoeff*fP_V_Omega;
						float fI_OLD = I[iRayCt][i][j];
						I[iRayCt][i][j] = f_Num/(f_DEN + f_SMALL);
						float f_diff = fabs(I[iRayCt][i][j] - fI_OLD)/(I[iRayCt][i][j] + f_SMALL);
						f_DMAX = max(f_DMAX,f_diff);
						

						//float fResidual = (Ip_New[iRayCt][i][j] - Ip[iRayCt][i][j] ) / 
											Ip_New[iRayCt][i][j];
						/*if(fResidual < fConvergence)
						{	
							pbRayConverged[iRayCt] = true;
							iConvergedRayCt++;
						}
						else
						{
							Ip[iRayCt][i][j] = Ip_New[iRayCt][i][j];
						}
					}
				}*/
			if(f_DMAX < fConvergence)
			{
					pbRayConverged[iRayCt] = true;
					iConvergedRayCt++;
			}
		}
	}while( (inbRays != iConvergedRayCt) /*|| (iNbIterations < 5)*/);


	// to calculate wall heat flux at bottom
	for( i = l_nbFaces/2; i < l_nbFaces; i+=nj)
	{
		float q_w = 0.0;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			Point p_dAve = vec_Dfm_[iRayCt];
			float fdp_s = 0.0;
			dotProduct(p_dAve,pnfs,fdp_s);
			q_w += pFaces[i].m_d_vec_face_ray_Inten[iRayCt]*fdp_s;
		}
		q_w = q_w/(dSBConst*pow(fTemp,4));
		printf("%lf\t",q_w);		
	}
	if(pFaces)
	{
		delete pFaces; pFaces = NULL; 
	}
	if(pCells)
	{
		delete pCells; pCells = NULL;
	}
// the end
}

// This function is written to accomodate only cell values , it will include face values inside
void radiativeTransferEquation::setBoundaryCondition(float **&T, float ***&I)
{
	int inbRays = 0;
	extractDataAndGetRayParameters(inbRays);
	if(!inbRays)
	{
		return;
	}
	int i = 0;
	int j = 0;
	float *pfXCV = NULL;
	pfXCV = new float [ni];

	float *pfYCV = NULL;
	pfYCV = new float [nj];

	float **ppf_Vol = NULL;
	ppf_Vol = new float*[ni];
	for(j = 0; j < ni; ++j)
	{
		ppf_Vol[j] = new float[nj];
	}
	Face *pFaces = NULL;
	Cell *pCells = NULL;
	//constructGridUsingCellAndFaces(inbRays, pFaces,pCells);
	constructGridForControlVolume(pfXCV, pfYCV, ppf_Vol);
	//float ***Ip = NULL;
	//float ***Ip_New = NULL;
	
	// Allocate Memory
	T = new float*[ni];
	for(i = 0; i < ni; ++i)
	{
		T[i] = new float[nj];
	}

	// Allocate Memory for Intensity
	I = new float**[inbRays];
	for(i = 0; i < inbRays; ++i)
	{
		I[i] = new float*[ni];
		for(j = 0; j < ni; ++j)
		{
			I[i][j] = new float[nj];
		}
	}

	// Allocate Memory for Intensity as internal field
	/*Ip = new float**[inbRays];
	for(i = 0; i < inbRays; ++i)
	{
		Ip[i] = new float*[ni];
		for(j = 0; j < ni; ++j)
		{
			Ip[i][j] = new float[nj];
		}
	}
	// Allocate Memory for New Intensity as internal field
	Ip_New = new float**[inbRays];
	for(i = 0; i < inbRays; ++i)
	{
		Ip_New[i] = new float*[ni];
		for(j = 0; j < ni; ++j)
		{
			Ip_New[i][j] = new float[nj];
		}
	}

	// calculate all face centre values for the intensity

	float ***pfI_face_ij = NULL;
	pfI_face_ij = new float**[inbRays];
	for(i = 0; i < inbRays; ++i)
	{
		pfI_face_ij[i] = new float*[ni-1];
		for(j = 0; j < ni - 1; ++j)
		{
			pfI_face_ij[i][j] = new float[nj];
		}
	}
	
	float ***pfI_face_ji = NULL;
	pfI_face_ji = new float**[inbRays];
	for(i = 0; i < inbRays; ++i)
	{
		pfI_face_ji[i] = new float*[ni];
		for(j = 0; j < ni; ++j)
		{
			pfI_face_ji[i][j] = new float[nj-1];
		}
	}*/
	

	int iRayCt = 0;
	// Initialize T everywhere in the domain
	for(i = 0; i < ni; ++i)
	{
		if( (i == 0) || (i == ni - 1) ) // Boundary field
		{
			for(j = 0; j < nj; ++j)
			{
				T[i][j] = 0.0;
				for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
				{
					I[iRayCt][i][j] = 0.0;
					//Ip[iRayCt][i][j] = 0.0;
				}
			}
		}
		else
		{
			for(j = 0; j < nj; ++j) 
			{
				if( (j == 0) || (j == nj-1) ) // Boundary field
				{
					T[i][j] = 0.0;
					for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
					{
						I[iRayCt][i][j] = 0.0;
						//Ip[iRayCt][i][j] = 0.0;
					}
				}
				else // internal field
				{
					T[i][j] = fTemp;
					for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
					{
						I[iRayCt][i][j] = 0.0;
						//Ip[iRayCt][i][j] = 0.0;
					}
				}
			}
		}
		
	}
	
	// Print T for the whole domain
	/*for(i = 0; i < ni; ++i)
	{
		for(j = 0; j < nj; ++j)
		{
			printf(" %lf\t",T[i][j]);
			for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
			{
				printf(" %lf\t",I[iRayCt][i][j]);
			}
			//std::cout<<T[i][j]<<"\t";
		}
		printf(" \n");
	}
	printf(" T is initialized properly");
	// Delete pointers to release memory 
	
	if(T)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; i++)
		{
			if(T[i])
			{
				delete []T[i]; T[i] = NULL;
			}
		}
		delete []T; T = NULL;
	}
	if(I)
	{
		//delete []T; T = NULL;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(I[iRayCt])
			{
				for(i = 0; i < ni; i++)
				{
					if(I[iRayCt][i])
					{
						delete []I[iRayCt][i]; I[iRayCt][i] = NULL;
					}
				}
				delete []I[iRayCt]; I[iRayCt] = NULL;
			}
			
		}
		delete []I; I = NULL;
	}


	if(Ip)
	{
		//delete []T; T = NULL;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(Ip[iRayCt])
			{
				for(i = 0; i < ni; i++)
				{
					if(Ip[iRayCt][i])
					{
						delete []Ip[iRayCt][i]; Ip[iRayCt][i] = NULL;
					}
				}
				delete []Ip[iRayCt]; Ip[iRayCt] = NULL;
			}
			
		}
		delete []Ip; Ip = NULL;
	}*/
		
	Point pnfe,pnfw,pnfn,pnfs;
	pnfe.x = 1.0;
	pnfe.y = 0.0;
	pnfe.z = 0.0;
	
	pnfw.x = -1.0;
	pnfw.y = 0.0;
	pnfw.z = 0.0;

	pnfn.x = 0.0;
	pnfn.y = 1.0;
	pnfn.z = 0.0;

	pnfs.x = 0.0;
	pnfs.y = -1.0;
	pnfs.z = 0.0;

	float fConvergence = 1e-06; 
	bool *pbRayConverged = NULL;
	pbRayConverged = new bool[inbRays];
	for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
	{
		pbRayConverged[iRayCt] = false;
	}
	int iConvergedRayCt = 0;
	int iNbIterations = 0;
	float f_SMALL = 1e-10;
	float f_DMAX = -1000;
	do
	{
		iNbIterations++;
		iConvergedRayCt = 0;
		f_DMAX = -1000;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(pbRayConverged[iRayCt])
			{
				iConvergedRayCt++;
				continue;
			}
			float fOmega = vec_omega_[iRayCt];
			Point p_cent = vec_Cent_[iRayCt];
			Point p_dAve = vec_Dfm_[iRayCt];
			float fVol = 1.0;
			float fAbsCoeff = 1.0;
			float fdp_e = 0.0;
			float fdp_w = 0.0;
			float fdp_n = 0.0;
			float fdp_s = 0.0;
			dotProduct(p_dAve,pnfe,fdp_e);
			dotProduct(p_dAve,pnfw,fdp_w);
			dotProduct(p_dAve,pnfn,fdp_n);
			dotProduct(p_dAve,pnfs,fdp_s);
			
				for(i = 1; i < ni - 1; ++i)
				{
					for(j = 1; j < nj - 1; ++j)
					{
						float fP_V_Omega = ppf_Vol[i-1][j-1]*fOmega; // product of volume and solid angle
						float f_Num = fabs(p_dAve.x*pfYCV[j-1])*I[iRayCt][i-1][j] + fabs(p_dAve.y*pfXCV[i-1])*I[iRayCt][i][j-1] + (fP_V_Omega*fAbsCoeff*dSBConst*pow(T[i][j],4)/pi);
						//Ip_New[iRayCt][i][j] = (fAbsCoeff*dSBConst*pow(T[i][j],4)/pi) - ( ( (I[iRayCt][i+1][j]*fdp_e) + (I[iRayCt][i-1][j]*fdp_w) + (I[iRayCt][i][j+1]*fdp_n) +
							//				 (I[iRayCt][i][j-1]*fdp_s) ) / (fOmega*fVol) );
						
						float f_DEN = fabs(p_dAve.x*pfYCV[j-1]) + fabs(p_dAve.y*pfXCV[i-1]) + fAbsCoeff*fP_V_Omega;
						float fI_OLD = I[iRayCt][i][j];
						I[iRayCt][i][j] = f_Num/(f_DEN + f_SMALL);
						float f_diff = fabs(I[iRayCt][i][j] - fI_OLD)/(I[iRayCt][i][j] + f_SMALL);
						f_DMAX = max(f_DMAX,f_diff);
						

						//float fResidual = (Ip_New[iRayCt][i][j] - Ip[iRayCt][i][j] ) / Ip_New[iRayCt][i][j];
						/*if(fResidual < fConvergence)
						{	
							pbRayConverged[iRayCt] = true;
							iConvergedRayCt++;
						}
						else
						{
							Ip[iRayCt][i][j] = Ip_New[iRayCt][i][j];
						}*/
					}
				}
				if(f_DMAX < fConvergence)
				{
						pbRayConverged[iRayCt] = true;
						iConvergedRayCt++;
				}
		}

		// Update the boundary field using internal field of intensity

		/*for(j = 0; j < nj; ++j) // this will update the left and right boundaries
		{
		
		
			for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
			{
				float fdp_e = 0.0;
				float fdp_w = 0.0;
				Point p_cent = vec_Cent_[iRayCt];
				dotProduct(p_cent,pnfe,fdp_e);
				dotProduct(p_cent,pnfw,fdp_w);
				if(fdp_w < 0)
				{
					I[iRayCt][0][j] = I[iRayCt][1][j];
				}
				if(fdp_e < 0)
				{
					I[iRayCt][ni-1][j] = I[iRayCt][ni-2][j];
				}
			}
		}


		for(i = 0; i < ni; ++i) // this will update the top and bottom boundaries
		{
			for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
			{
				float fdp_n = 0.0;
				float fdp_s = 0.0;
				Point p_cent = vec_Cent_[iRayCt];
				dotProduct(p_cent,pnfn,fdp_n);
				dotProduct(p_cent,pnfs,fdp_s);
				if(fdp_s < 0)
				{
					I[iRayCt][i][0] = I[iRayCt][i][1];
				}
				if(fdp_n < 0)
				{
					I[iRayCt][i][nj-1] = I[iRayCt][i][nj-2];
				}
			}
		}*/
		if(inbRays == iConvergedRayCt)
		{
			printf("All rays are converged successfully in %d number of iterations ",iNbIterations);
			break;
		}
	}while( (inbRays != iConvergedRayCt) /*|| (iNbIterations < 5)*/);
	

	// Update the boundary field using internal field of intensity

	
	for(j = 0; j < nj; ++j) // this will update the left and right boundaries
	{
		
		
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			float fdp_e = 0.0;
			float fdp_w = 0.0;
			Point p_cent = vec_Cent_[iRayCt];
			dotProduct(p_cent,pnfe,fdp_e);
			dotProduct(p_cent,pnfw,fdp_w);
			if(fdp_w < 0)
			{
				I[iRayCt][0][j] = I[iRayCt][1][j];
			}
			if(fdp_e < 0)
			{
				I[iRayCt][ni-1][j] = I[iRayCt][ni-2][j];
			}
		}
	}


	for(i = 0; i < ni; ++i) // this will update the top and bottom boundaries
	{
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			float fdp_n = 0.0;
			float fdp_s = 0.0;
			Point p_cent = vec_Cent_[iRayCt];
			dotProduct(p_cent,pnfn,fdp_n);
			dotProduct(p_cent,pnfs,fdp_s);
			if(fdp_s < 0)
			{
				I[iRayCt][i][0] = I[iRayCt][i][1];
			}
			if(fdp_n < 0)
			{
				I[iRayCt][i][nj-1] = I[iRayCt][i][nj-2];
			}
		}
	}
	//printf("total number of iterations:%d ",iNbIterations);
	
	/*for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
	{
		for(i = 0; i < ni - 1; ++i)
		{
			for(j = 0; j < nj; ++j)
			{
				pfI_face_ij[iRayCt][i][j] =  Ip[iRayCt][i+1][j] - Ip[iRayCt][i][j] ;
			}
		}
	}

	for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
	{
		for(i = 0; i < ni; ++i)
		{
			for(j = 0; j < nj - 1; ++j)
			{
				pfI_face_ji[iRayCt][i][j] =  Ip[iRayCt][i][j+1] - Ip[iRayCt][i][j] ;
			}
		}
	}*/


	// CALCULATION OF INCIDENT RADIATION ENERGY(G) 
	float **fG_total = NULL;
	fG_total = new float*[ni];
	for(i = 0; i < ni; ++i)
	{
		fG_total[i] = new float[nj];
	}
	
	//float fI_total = 0.0;
	//printf("for ray ID = %d\n",iRayCt);
	printf("\n");
	printf(" CALCULATION OF INCIDENT RADIATION INTENSITY(I) \n");
	for(i = 0; i < ni; ++i)
	{
		for(j = 0; j < nj; ++j)
		{
			printf("I for i = %d,j = %d \n",i,j);
			//float fG_sum = 0.0;
			for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
			{
				printf("%lf\t",I[iRayCt][i][j]);
				//float fOmega = vec_omega_[iRayCt];
				//fG_sum += I[iRayCt][i][j]*fOmega;
				//fG_total[i][j] += I[iRayCt][i][j]*fOmega;
			}
			//fG_total[i][j] = fG_sum;
			//fI_total += Ip[iRayCt][i][j];
			//printf("%lf\t",fG_total[i][j]);
			printf("\n");
		}
		printf("\n");
	}
	
	printf(" CALCULATION OF INCIDENT RADIATION ENERGY(G) \n");

	for(i = 0; i < ni; ++i)
	{
		for(j = 0; j < nj; ++j)
		{
			float fG_sum = 0.0;
			for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
			{
				float fOmega = vec_omega_[iRayCt];
				fG_sum += I[iRayCt][i][j]*fOmega;
				//fG_total[i][j] += I[iRayCt][i][j]*fOmega;
			}
			fG_total[i][j] = fG_sum;
			//fI_total += Ip[iRayCt][i][j];
			printf("%lf\t",fG_total[i][j]);
		}
		printf("\n");
	}


	// CALCULATION OF RADIATIVE FLUX(W/m2) 
	float **ppfQ_rad_PY = NULL;
	ppfQ_rad_PY = new float*[ni];
	for(i = 0; i < ni; ++i)
	{
		ppfQ_rad_PY[i] = new float[nj];
	}

	float **ppfQ_rad_MY = NULL;
	ppfQ_rad_MY = new float*[ni];
	for(i = 0; i < ni; ++i)
	{
		ppfQ_rad_MY[i] = new float[nj];
	}

	float **ppfQ_rad_PX = NULL;
	ppfQ_rad_PX = new float*[ni];
	for(i = 0; i < ni; ++i)
	{
		ppfQ_rad_PX[i] = new float[nj];
	}

	float **ppfQ_rad_MX = NULL;
	ppfQ_rad_MX = new float*[ni];
	for(i = 0; i < ni; ++i)
	{
		ppfQ_rad_MX[i] = new float[nj];
	}

	float *pfQ_rad_bottom = NULL;
	pfQ_rad_bottom = new float[ni];

	float *pfQ_rad_top = NULL;
	pfQ_rad_top = new float[ni];

	float *pfQ_rad_left = NULL;
	pfQ_rad_left = new float[nj];

	float *pfQ_rad_right = NULL;
	pfQ_rad_right = new float[nj];
	
	
	//float fI_total = 0.0;
	//printf("for ray ID = %d\n",iRayCt);
	
	printf(" CALCULATION OF RADIATIVE FLUX(W/m2) \n");

	// This is done in the same way HFLUX Subroutine has been called in fortran code, but then i found a different way to do it, because finally 
	// they have to update the heat flux at the boundaries only
	/*for(i = 0; i < ni; ++i)
	{
		for(j = 0; j < nj; ++j)
		{
			float f_SQPY = 0.0;
			float f_SQMY = 0.0;
			float f_SQPX = 0.0;
			float f_SQMX = 0.0;
			for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
			{
				Point p_dfm = vec_Dfm_[iRayCt];
				float f_ADCY = fabs(p_dfm.y);
				float f_ADCX = fabs(p_dfm.x);
				if(i > 0 && i < ni - 1)
				{
					if(p_dfm.y > 0)
					{
						f_SQPY += f_ADCY*I[iRayCt][i][j];
					}
					else
					{
						f_SQMY += f_ADCY*I[iRayCt][i][j];
					}
				}
				if(j > 0 && j < nj - 1)
				{
					if(p_dfm.x > 0)
					{
						f_SQPX += f_ADCX*I[iRayCt][i][j];
					}
					else
					{
						f_SQMX += f_ADCX*I[iRayCt][i][j];
					}
				}
				//float fOmega = vec_omega_[iRayCt];
				//fG_total[i][j] += Ip[iRayCt][i][j]*fOmega;
			}
			//fI_total += Ip[iRayCt][i][j];
			
			ppfQ_rad_PY[i][j] = f_SQPY;
			ppfQ_rad_MY[i][j] = f_SQMY;
			ppfQ_rad_PX[i][j] = f_SQPX;
			ppfQ_rad_MX[i][j] = f_SQMX;

			//printf("%lf\t",ppfQ_rad_PY[i][j]);
			//printf("%lf\t",ppfQ_rad_PY[i][j]);
		}
		//printf("\n");
	}*/

	// pRINT THE HEAT FLUX AT BOTTOM BOUNDARY
	// Update the boundary field using internal field of intensity

	for(j = 0; j < nj; ++j) // this will update the left and right boundaries
	{
		float f_SQPY_first = 0.0;
		//float f_SQMY_first = 0.0;
		float f_SQPY_last = 0.0;
		//float f_SQMY_last = 0.0;
		//	float f_SQPX = 0.0;
		//	float f_SQMX = 0.0;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			Point p_dfm = vec_Dfm_[iRayCt];
			float f_ADCY = /*fabs*/(p_dfm.y);
			Point p_cent = vec_Cent_[iRayCt];
			float fOmega = vec_omega_[iRayCt];
			//float f_ADCX = fabs(p_dfm.x);
			if(-(p_cent.x) < 0)
			{
				f_SQPY_first += -(p_cent.x)*I[iRayCt][0][j]*fOmega;
			}
			if((p_cent.x) < 0)
			{
				f_SQPY_last += (p_cent.x)*I[iRayCt][ni-1][j]*fOmega;
			}
			/*else
			{
				f_SQMY_first += f_ADCY*I[iRayCt][0][j];
				f_SQMY_last += f_ADCY*I[iRayCt][ni-1][j];
			}*/
			//I[iRayCt][0][j] = I[iRayCt][1][j];
			//I[iRayCt][ni-1][j] = I[iRayCt][ni-2][j];
		}
		pfQ_rad_left[j] = fabs(f_SQPY_first/* - f_SQMY_first*/);
		pfQ_rad_right[j] = fabs(f_SQPY_last/* - f_SQMY_last*/);
	}


	for(i = 0; i < ni; ++i) // this will update the top and bottom boundaries
	{
		//_last = 0.0;
		float f_SQPX_first = 0.0;
		//float f_SQMX_first = 0.0;
		float f_SQPX_last = 0.0;
		//float f_SQMX_last = 0.0;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			Point p_dfm = vec_Dfm_[iRayCt];
			//float f_ADCY = fabs(p_dfm.y);
			float f_ADCX = /*fabs*/(p_dfm.x);
			Point p_cent = vec_Cent_[iRayCt];
			float fOmega = vec_omega_[iRayCt];
			if(-(p_cent.y) < 0)
			{
				f_SQPX_first += -(p_cent.y)*I[iRayCt][i][0]*fOmega;
			}
			if((p_cent.y) < 0)
			{	
				f_SQPX_last += (p_cent.y)*I[iRayCt][i][nj-1]*fOmega;
			}
			/*else
			{
				f_SQMX_first += f_ADCX*I[iRayCt][i][0];
				f_SQMX_last += f_ADCX*I[iRayCt][i][nj-1];
			}*/
			//I[iRayCt][0][j] = I[iRayCt][1][j];
			//I[iRayCt][ni-1][j] = I[iRayCt][ni-2][j];
		}
		//printf("f_SQPX_first = %lf\t ,f_SQPX_last=  %lf\t ,f_SQMX_first=  %lf\t ,f_SQMX_last=  %lf \n ",f_SQPX_first,f_SQPX_last ,f_SQMX_first, f_SQMX_last);
		pfQ_rad_bottom[i] = fabs(f_SQPX_first /*- f_SQMX_first*/);
		
		pfQ_rad_top[i] = fabs(f_SQPX_last/* - f_SQMX_last*/);

		printf("Q_bottom = %lf\t ,Q_top=  %lf\n ",pfQ_rad_bottom[i],pfQ_rad_top[i] );
	}


	/*for(i = 0; i < ni; ++i) 
	{
		ppfQ_rad_bottom[i] = fabs(ppfQ_rad_PY[i+1][1] - ppfQ_rad_MY[i+1][1]);
	}*/

	
	if(pfXCV)
	{
		delete []pfXCV; pfXCV = NULL;
	}
	if(pfYCV)
	{
		delete []pfYCV; pfYCV = NULL;
	}
	if(ppf_Vol)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(ppf_Vol[i])
			{
				delete []ppf_Vol[i]; ppf_Vol[i] = NULL;
			}
		}
		delete []ppf_Vol; ppf_Vol = NULL;
	}
	if(ppfQ_rad_PY)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(ppfQ_rad_PY[i])
			{
				delete []ppfQ_rad_PY[i]; ppfQ_rad_PY[i] = NULL;
			}
		}
		delete []ppfQ_rad_PY; ppfQ_rad_PY = NULL;
	}
	if(ppfQ_rad_MY)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(ppfQ_rad_MY[i])
			{
				delete []ppfQ_rad_MY[i]; ppfQ_rad_MY[i] = NULL;
			}
		}
		delete []ppfQ_rad_MY; ppfQ_rad_MY = NULL;
	}
	if(ppfQ_rad_PX)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(ppfQ_rad_PX[i])
			{
				delete []ppfQ_rad_PX[i]; ppfQ_rad_PX[i] = NULL;
			}
		}
		delete []ppfQ_rad_PX; ppfQ_rad_PX = NULL;
	}
	if(ppfQ_rad_MX)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(ppfQ_rad_MX[i])
			{
				delete []ppfQ_rad_MX[i]; ppfQ_rad_MX[i] = NULL;
			}
		}
		delete []ppfQ_rad_MX; ppfQ_rad_MX = NULL;
	}
	if(fG_total)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(fG_total[i])
			{
				delete []fG_total[i]; fG_total[i] = NULL;
			}
		}
		delete []fG_total; fG_total = NULL;
	}
	if(T)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; ++i)
		{
			if(T[i])
			{
				delete []T[i]; T[i] = NULL;
			}
		}
		delete []T; T = NULL;
	}
	if(I)
	{
		//delete []T; T = NULL;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(I[iRayCt])
			{
				for(i = 0; i < ni; ++i)
				{
					if(I[iRayCt][i])
					{
						delete []I[iRayCt][i]; I[iRayCt][i] = NULL;
					}
				}
				delete []I[iRayCt]; I[iRayCt] = NULL;
			}
			
		}
		delete []I; I = NULL;
	}


	/*if(Ip)
	{
		//delete []T; T = NULL;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(Ip[iRayCt])
			{
				for(i = 0; i < ni; ++i)
				{
					if(Ip[iRayCt][i])
					{
						delete []Ip[iRayCt][i]; Ip[iRayCt][i] = NULL;
					}
				}
				delete []Ip[iRayCt]; Ip[iRayCt] = NULL;
			}
			
		}
		delete []Ip; Ip = NULL;
	}

	if(Ip_New)
	{
		//delete []T; T = NULL;
		for(iRayCt = 0; iRayCt < inbRays; ++iRayCt)
		{
			if(Ip_New[iRayCt])
			{
				for(i = 0; i < ni; ++i)
				{
					if(Ip_New[iRayCt][i])
					{
						delete []Ip_New[iRayCt][i]; Ip_New[iRayCt][i] = NULL;
					}
				}
				delete []Ip_New[iRayCt]; Ip_New[iRayCt] = NULL;
			}
			
		}
		delete []Ip_New; Ip_New = NULL;
	}*/
}
void radiativeTransferEquation::constructGridForControlVolume(float *&pfXCV, float *&pfYCV, float **&ppf_Vol)
{
	// I am trying to construct the grid the way Chai has made using Fortran for RAT
	// X is the cell-centre value at the location
	// XU is the face-centre value  
	float X[ni],Y[nj],XU[ni],YV[nj], XCV[ni], YCV[nj], Vol[ni][nj];
	// XU is the value of X at control volume face
	//int L1 = ni + 2; // ni is the number of control volumes, but we need to start the grid formation from the face
	//int L2 = L1 -1;
	int XL = 1; // Total length in x-direction
 	int YL = 1;// Total length in y-direction
	//XU[0] = 0.0; // value at starting face
	//XU[L1] = XL;
	float dX = (float)XL/ni;
	float dY = (float)YL/nj;
	int i = 0;
	int j  =0;	

	// CALCULATION OF X-DIRECTION GRID
	for(i = 0; i < ni; ++i)
	{
		XU[i] = (i+1)*dX;
		if(i == 0 )
		{
			X[i] = XU[i]/2.0;
			XCV[i] = XU[i];
		}
		else
		{
			X[i] = 0.5*(XU[i] + XU[i-1]);
			XCV[i] = XU[i] - XU[i-1];
		}
		pfXCV[i] = XCV[i];
	}

	//int M1 = nj + 2;
	//int M2 = M1 -1;
	//YV[0] = 0.0; // value at starting face
	//YV[M1] = YL;

	// CALCULATION OF Y-DIRECTION GRID
	for(j = 0; j < nj; ++j)
	{
		YV[j] = (j+1)*dY;
		if(j == 0 )
		{
			Y[j] = YV[j]/2.0;
			YCV[j] = YV[j];
		}
		else
		{
			Y[j] = 0.5*(YV[j] + YV[j-1]);
			YCV[j] = YV[j] - YV[j-1];
		}
		pfYCV[j] = YCV[j];
	}
		
	// CALCULATIONS OF CONTROL VOLUME VOLUMES
	for(i = 0; i < ni; ++i)
	{
		for(j = 0; j < nj; ++j)
		{
			Vol[i][j] = XCV[i]*YCV[j];
			ppf_Vol[i][j] = XCV[i]*YCV[j];
		}
	}
	

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
radiativeTransferEquation::radiativeTransferEquation
(
    const scalar nPhi,
    const scalar nTheta,
    const float radius,
    const Point centrePt
)
:
    nPhi_(nPhi),
    nTheta_(nTheta),
    radius_(radius),
    centrePt_(centrePt)
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

radiativeTransferEquation::~radiativeTransferEquation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void radiativeTransferEquation::crossProduct(Point a, Point b, Point &cp)
{
  cp.x = 0.0;
  cp.y = 0.0;
  cp.z = 0.0;
  cp.x = a.y*b.z - b.y*a.z;
  cp.y = b.x*a.z - a.x*b.z;
  cp.z = a.x*b.y - b.x*a.y;
}
void radiativeTransferEquation::getSubtractionVector(Point a, Point b, Point &sv)
{
 	sv.x = 0.0;
	sv.y = 0.0;
	sv.z = 0.0;
	sv.x = a.x - b.x;
	sv.y = a.y - b.y;
	sv.z = a.z - b.z;
}
//------------------------
//-- Entry point. This main() function demonstrates how you can
//-- use "printStandardSphere()", but you probably won't
//-- want/need to copy it in your own code.
 
int main(int argc, char *argv[])
{  

	float **T = NULL;
	float ***I = NULL;
	 //setBoundaryConditions(T);
	
  int nVertical  = 4;                  // Number vertical lines.
  int nHorizontal = 2;      // Number horizontal lines.
  // NOTE: for a good sphere use ~half the number of longitude lines than latitude.
  Point centerPt;            // Position the center of out sphere at (0,0,0).
  centerPt.x =0.0;
  centerPt.y = 0.0;
  centerPt.z = 0.0;
  /*if (argc < 2) {
    fprintf(stderr, "Must enter: './programname outputfile.obj'\n");
    return (-1);
  }*/
 
  FILE *fout = fopen(argv[1] , "w");
  /*if (fout == NULL) {
     printf("Couldn't open output file %s.\n", argv[1]);
     return (-1);
  }*/

  radiativeTransferEquation pnewSphere(nHorizontal,nVertical,1,centerPt);
	//pnewSphere.extractDataAndGetRayParameters();
	//pnewSphere.setBoundaryCondition(T, I);
	pnewSphere.updateIntensity();
	// Delete pointers to release memory 
	int i = 0;
	if(T)
	{
		//delete []T; T = NULL;
		for(i = 0; i < ni; i++)
		{
			if(T[i])
			{
				delete []T[i]; T[i] = NULL;
			}
		}
		delete []T; T = NULL;
	}
  //pnewSphere.createSphere(fout);
  //fclose(fout);
  //fprintf(stdout, "  # vertices:   %d\n", numVertices);
  return (0);
}


// ************************************************************************* //
void radiativeTransferEquation::extractDataAndGetRayParameters(int &inbRays)
{
	//string mainPath = path();
	//Info<<"mainPath : "<<mainPath<<endl;
	
	
	//string filePath = "out.obj";
	//string filePath = path()/"out_rec1.obj";
	//Info<<"filePath : "<<filePath<<endl;
	// Open your file
	ifstream someStream("out.obj");
    

    // Set up a place to store our data read from the file
    string line;

    // Read and throw away the first line simply by doing
    // nothing with it and reading again
    getline( someStream, line );
    getline( someStream, line );
 
	string s = line;
	string delimiter = ",";

	size_t pos = 0;
	string token;
	int total = 0;
	int nodes=0;
	int elements=0;
	int count=0;
	int iElemType = 0;
	std::vector<Point>vec_UpperTrianglePoints;
  	std::vector<Point>vec_LowerTrianglePoints;
  	std::vector<Point>vec_QuadPoints;

	std::vector<std::vector<int> >vec_UpperTriangleNodeIndicies;
	std::vector<std::vector<int> >vec_LowerTriangleNodeIndicies;
	std::vector<std::vector<int> >vec_QuadNodeIndicies;

	
	std::vector<std::pair<Point, float> >vec_pair_CentroidDfm_;

	//std::vector<vector> vec_lower_Dfm;

	//std::vector<vector> vec_lower_Cent;

	//std::vector<scalar> vec_lower_omega;

	//Foam::label iFieldSize = 0;
	//label lRayCt = 0;
	double dSum = 0.0;
	while (getline(someStream,s))
	{
    	if(count==total)
		{
       	 	total=0;
        	while ((pos = s.find(delimiter)) != string::npos)
			{	
            	token = s.substr(0, pos);
            	if (token.find("NODES") != string::npos)
				{
            	    nodes=atoi(token.substr(token.find("NODES=")+6 ).c_str());
               		total+=nodes;
            	}
            	if (token.find("ELEMENTS") != string::npos)
				{
                	elements=atoi(token.substr(token.find("ELEMENTS=") + 9).c_str());
                	total+=elements;
            	}	
            	s.erase(0, pos + delimiter.length());
        	}
        	count=0;
			iElemType++;
        	continue;
    	}
    	istringstream iss(s);
    	double a[4];  //Cordinates x=a[0],y=a[1],z=a[2]
    	int i=0;
    	do
    	{
        	string subs;
        	iss >> subs;
        	a[i++]=atof(subs.c_str());
    	}while (iss);
    	iss.clear();
		Point Obj;
		Obj.x = a[0];
		Obj.y = a[1];
		Obj.z = a[2];
		if(iElemType == 1) // upper triangles
		{			
			if(count < nodes)
			{
				vec_UpperTrianglePoints.push_back(Obj);
			}
			else if(count >= nodes && count < total)
			{
				std::vector<int>vec_currentFace;
				vec_currentFace.push_back(a[0]);
				vec_currentFace.push_back(a[1]);
				vec_currentFace.push_back(a[2]);
				if(vec_currentFace.size())
				{
					vec_UpperTriangleNodeIndicies.push_back(vec_currentFace);
				}
				Point fv13, fv23, cp;
    			getSubtractionVector(vec_UpperTrianglePoints[a[0]-1],vec_UpperTrianglePoints[a[1]-1],fv13);
    			getSubtractionVector(vec_UpperTrianglePoints[a[2]-1],vec_UpperTrianglePoints[a[1]-1],fv23);

				
				//vector vec_13, vec_23, cp;
				//vec_13 = vector(fv13.x, fv13.y, fv13.z);
				//vec_23 = vector(fv23.x, fv23.y, fv23.z);
				//cp = 0.5*(vec_13 ^ vec_23);
    			Point centroid;
    			centroid.x = (vec_UpperTrianglePoints[a[0]-1].x+vec_UpperTrianglePoints[a[1]-1].x+vec_UpperTrianglePoints[a[2]-1].x)/3.0;
    			centroid.y = (vec_UpperTrianglePoints[a[0]-1].y+vec_UpperTrianglePoints[a[1]-1].y+vec_UpperTrianglePoints[a[2]-1].y)/3.0;
    			centroid.z = (vec_UpperTrianglePoints[a[0]-1].z+vec_UpperTrianglePoints[a[1]-1].z+vec_UpperTrianglePoints[a[2]-1].z)/3.0;
				//float magCp = sqrt(fabs(cp[0]*cp[0])+fabs(cp[1]*cp[1])+fabs(cp[2]*cp[2]));

				crossProduct(fv13, fv23, cp);
				cp.x = 0.5*cp.x;
				cp.y = 0.5*cp.y;
				cp.z = 0.5*cp.z;
				float magCp = sqrt(fabs(cp.x*cp.x)+fabs(cp.y*cp.y)+fabs(cp.z*cp.z));
				printf("Omega : %lf\n",magCp);
				//Info << "Omega : "<<magCp<<endl;
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				//float fDotProduct = cp & vector(centroid.x, centroid.y, centroid.z);
				float fDotProduct = 0.0;				
				dotProduct(cp,centroid,fDotProduct);
				if(fDotProduct < 0)	
				{
					centroid.x = centroid.x*(-1);
					centroid.y = centroid.y*(-1);
					centroid.z = centroid.z*(-1);
				} 
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				dSum += magCp;
				std::pair<Point, float>li_pair;
				li_pair.first = centroid;
				li_pair.second = magCp;
				vec_pair_CentroidDfm_.push_back(li_pair);
				Point p_mod_cp, p_mod_cent;
				p_mod_cp.x = cp.z;
				p_mod_cp.y = cp.x;
				p_mod_cp.z = cp.y;

				p_mod_cent.x = centroid.z;
				p_mod_cent.y = centroid.x;
				p_mod_cent.z = centroid.y;
				vec_Dfm_.push_back(p_mod_cp);
				vec_Cent_.push_back(p_mod_cent);
				vec_omega_.push_back(magCp);
				if (omegaMax_ <  magCp)
        		{
            		omegaMax_ = magCp;
        		}
				//dAve_.boundaryField()[lRayCt] = (cp);
				//lRayCt++;
				//iFieldSize++;
			}
				
    	}
		else if(iElemType == 2) // Quad
		{
			if(count < nodes)
			{
				vec_QuadPoints.push_back(Obj);
			}
			else if(count >= nodes && count < total)
			{
				std::vector<int>vec_currentFace;
				vec_currentFace.push_back(a[0]);
				vec_currentFace.push_back(a[1]);
				vec_currentFace.push_back(a[2]);
				vec_currentFace.push_back(a[3]);
				if(vec_currentFace.size())
				{
					vec_QuadNodeIndicies.push_back(vec_currentFace);
				}
				Point centroid;
				centroid.x = (vec_QuadPoints[a[0]-1].x+vec_QuadPoints[a[1]-1].x+vec_QuadPoints[a[2]-1].x+vec_QuadPoints[a[3]-1].x)/4.0;
     			centroid.y = (vec_QuadPoints[a[0]-1].y+vec_QuadPoints[a[1]-1].y+vec_QuadPoints[a[2]-1].y+vec_QuadPoints[a[3]-1].y)/4.0;
      			centroid.z = (vec_QuadPoints[a[0]-1].z+vec_QuadPoints[a[1]-1].z+vec_QuadPoints[a[2]-1].z+vec_QuadPoints[a[3]-1].z)/4.0;
				Point fv13, fv24, cp;
      			getSubtractionVector(vec_QuadPoints[a[2]-1],vec_QuadPoints[a[0]-1],fv13);
      			getSubtractionVector(vec_QuadPoints[a[3]-1],vec_QuadPoints[a[1]-1],fv24);
				crossProduct(fv13, fv24, cp);
				cp.x = 0.5*cp.x;
				cp.y = 0.5*cp.y;
				cp.z = 0.5*cp.z;
				float magCp = sqrt(fabs(cp.x*cp.x)+fabs(cp.y*cp.y)+fabs(cp.z*cp.z));
				printf("Omega : %lf\n",magCp);
				//Info << "Omega : "<<magCp<<endl;
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				//float fDotProduct = cp & vector(centroid.x, centroid.y, centroid.z);
				float fDotProduct = 0.0;				
				dotProduct(cp,centroid,fDotProduct);
				//vector vec_13, vec_24, cp;
    			//vec_13 = vector(fv13.x,fv13.y,fv13.z);
    			//vec_24 = vector(fv24.x,fv24.y,fv24.z);
      			//cp = 0.5*(vec_13 ^ vec_24);
				//float magCp = sqrt(fabs(cp[0]*cp[0])+fabs(cp[1]*cp[1])+fabs(cp[2]*cp[2]));
				dSum += magCp;
				//float fDotProduct = cp & vector(centroid.x, centroid.y, centroid.z);
				if(fDotProduct < 0)	
				{
					centroid.x = centroid.x*(-1);
					centroid.y = centroid.y*(-1);
					centroid.z = centroid.z*(-1);
				} 
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				std::pair<Point, float>li_pair;
				li_pair.first = centroid;
				li_pair.second = magCp;
				vec_pair_CentroidDfm_.push_back(li_pair);
				Point p_mod_cp, p_mod_cent;
				p_mod_cp.x = cp.z;
				p_mod_cp.y = cp.x;
				p_mod_cp.z = cp.y;

				p_mod_cent.x = centroid.z;
				p_mod_cent.y = centroid.x;
				p_mod_cent.z = centroid.y;
				vec_Dfm_.push_back(p_mod_cp);
				vec_Cent_.push_back(p_mod_cent);
				//Info << "Omega : "<<magCp<<endl;
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				vec_omega_.push_back(magCp);
				if (omegaMax_ <  magCp)
        		{
            		omegaMax_ = magCp;
        		}
				//dAve_.boundaryField()[lRayCt] = (cp);
				//lRayCt++;
				//vec_Dfm_[iFieldSize] = cp;
				//iFieldSize++;
			}	
    	}
		else if(iElemType == 3) // Lower Triangles
		{	
			if(count < nodes)
			{
				vec_LowerTrianglePoints.push_back(Obj);
			}
			else if(count >= nodes && count < total)
			{
				std::vector<int>vec_currentFace;
				vec_currentFace.push_back(a[0]);
				vec_currentFace.push_back(a[1]);
				vec_currentFace.push_back(a[2]);
				if(vec_currentFace.size())
				{
					vec_LowerTriangleNodeIndicies.push_back(vec_currentFace);
				}
				Point fv13, fv23, cp;
    			getSubtractionVector(vec_LowerTrianglePoints[a[0]-1],vec_LowerTrianglePoints[a[1]-1],fv13);
    			getSubtractionVector(vec_LowerTrianglePoints[a[2]-1],vec_LowerTrianglePoints[a[1]-1],fv23);
				crossProduct(fv13, fv23, cp);
				cp.x = 0.5*cp.x;
				cp.y = 0.5*cp.y;
				cp.z = 0.5*cp.z;
				float magCp = sqrt(fabs(cp.x*cp.x)+fabs(cp.y*cp.y)+fabs(cp.z*cp.z));
				printf("Omega : %lf\n",magCp);
				//Info << "Omega : "<<magCp<<endl;
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				//float fDotProduct = cp & vector(centroid.x, centroid.y, centroid.z);
				float fDotProduct = 0.0;				
				
				//vector vec_13, vec_23, cp;
    			//vec_13 = vector(fv13.x,fv13.y,fv13.z);
    			//vec_23 = vector(fv23.x,fv23.y,fv23.z);
    			//cp = 0.5*(vec_13 ^ vec_23);
    			Point centroid;
    			centroid.x = (vec_LowerTrianglePoints[a[0]-1].x+vec_LowerTrianglePoints[a[1]-1].x+vec_LowerTrianglePoints[a[2]-1].x)/3.0;
    			centroid.y = (vec_LowerTrianglePoints[a[0]-1].y+vec_LowerTrianglePoints[a[1]-1].y+vec_LowerTrianglePoints[a[2]-1].y)/3.0;
    			centroid.z = (vec_LowerTrianglePoints[a[0]-1].z+vec_LowerTrianglePoints[a[1]-1].z+vec_LowerTrianglePoints[a[2]-1].z)/3.0;
				dotProduct(cp,centroid,fDotProduct);
    			//float magCp = sqrt(fabs(cp[0]*cp[0])+fabs(cp[1]*cp[1])+fabs(cp[2]*cp[2]));
				dSum += magCp;
				std::pair<Point, float>li_pair;
				li_pair.first = centroid;
				li_pair.second = magCp;
				//float fDotProduct = cp & vector(centroid.x, centroid.y, centroid.z);
				if(fDotProduct < 0)	
				{
					centroid.x = centroid.x*(-1);
					centroid.y = centroid.y*(-1);
					centroid.z = centroid.z*(-1);
				} 
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				vec_pair_CentroidDfm_.push_back(li_pair);
				Point p_mod_cp, p_mod_cent;
				p_mod_cp.x = cp.z;
				p_mod_cp.y = cp.x;
				p_mod_cp.z = cp.y;

				p_mod_cent.x = centroid.z;
				p_mod_cent.y = centroid.x;
				p_mod_cent.z = centroid.y;
				vec_Dfm_.push_back(p_mod_cp);
				vec_Cent_.push_back(p_mod_cent);
				vec_omega_.push_back(magCp);
				//Info << "DotProduct : "<<(cp & vector(centroid.x, centroid.y, centroid.z))<<endl;
				//Info << "Omega : "<<magCp<<endl;
				if (omegaMax_ <  magCp)
        		{
            		omegaMax_ = magCp;
        		}
				//dAve_.boundaryField()[lRayCt] = (cp);
				//lRayCt++;
				//vec_Dfm_[iFieldSize] = cp;
				//iFieldSize++;
			}
				
    	}
    	count++;  
	}
	inbRays = vec_omega_.size();
}
// ************************************************************************* //
void radiativeTransferEquation::dotProduct(Point a, Point b, float &dp)
{
	dp = 0.0;
	dp = (a.x*b.x) +  (a.y*b.y) + (a.z*b.z);
}

// ************************************************************************* //
void radiativeTransferEquation::constructGridUsingCellAndFaces(long l_nbRays, Face *&pFaces, 
															   Cell *&pCells, long &l_nbCells,
															   long &l_nbFaces)
{

	// first assign face parameters
	// Get to know the total number of faces
	long i = 0;
	long j = 0;
	long l_RayCt = 0;
	l_nbFaces = (ni*(ni+1)) + (nj*(nj+1));

	l_nbCells = ni*nj;
	pCells = new Cell[l_nbCells];
	long l_CellCount = 0;

	// First to know which cells are adjecent to particular cell
	for(j = 0; j < nj; ++j)
	{
		for(i = 0; i < ni; ++i)
		{
			/*if(i == 0 || i == ni -1) // this is boundary
			{
				// here west and south cells will be -1
				//pCells[l_CellCount].m_l_CellIndex_West = -1;
				pCells[l_CellCount].m_l_CellIndex_South = -1;
			}
			else
			{*/
				//pCells[l_CellCount].m_l_CellIndex_West = l_CellCount-1;
				if( j == 0)
				{
					pCells[l_CellCount].m_l_CellIndex_South = -1;
				}
				else
				{
					pCells[l_CellCount].m_l_CellIndex_South = l_CellCount - ni;
				}
			/*}*/
			if(i != ni -1)
			{
				pCells[l_CellCount].m_l_CellIndex_East = l_CellCount + 1;
			}
			else
			{
				pCells[l_CellCount].m_l_CellIndex_East = -1;				
			}
			if(i ==0)
			{
				pCells[l_CellCount].m_l_CellIndex_West = -1;
			}
			else
			{
				pCells[l_CellCount].m_l_CellIndex_West = l_CellCount-1;
			}
			if(j != nj - 1)
			{
				pCells[l_CellCount].m_l_CellIndex_North = (ni*(j+1)) + i;
			}
			else
			{
				pCells[l_CellCount].m_l_CellIndex_North = -1;
			}
			l_CellCount++;
		}
	}
	
	pFaces = new Face[l_nbFaces];	
	long l_faceCount = 0; 
	for(i = 0; i < ni; ++i)
	{
		for(j = 0; j < ni+1; ++j)
		{
			if(j ==0 || j == ni)
			{
				pFaces[l_faceCount].m_b_IfBoundary = true;
				pFaces[l_faceCount].m_f_Temp = 0.0;
			}
			else
			{
				pFaces[l_faceCount].m_b_IfBoundary = false;
			}
			l_faceCount++;
		}
	}

	for(j = 0; j < nj; ++j)
	{
		for(i = 0; i < nj+1; ++i)
		{
			if(i ==0 || i == nj)
			{
				pFaces[l_faceCount].m_b_IfBoundary = true;
				pFaces[l_faceCount].m_f_Temp = 0.0;
			}
			else
			{
				pFaces[l_faceCount].m_b_IfBoundary = false;
			}
			l_faceCount++;
		}
	}
	for(i = 0; i < l_nbFaces; ++i)
	{
		pFaces[i].m_l_nbRays = l_nbRays;
		if(pFaces[i].m_d_vec_face_ray_Inten.size())
		{
			pFaces[i].m_d_vec_face_ray_Inten.clear();
		}
		pFaces[i].m_d_vec_face_ray_Inten.resize(l_nbRays);
		for(l_RayCt = 0; l_RayCt < l_nbRays; ++l_RayCt)
		{
			pFaces[i].m_d_vec_face_ray_Inten[l_RayCt] = 0.0;
		}
		pFaces[i].m_f_Area = (1/(float)ni);
		pFaces[i].m_f_Temp = 0.0;
	}
	// In order to distinguish between boundary and internal faces, this loop is necessary
	/*for(j = 0; j < nj + 1; ++j)
	{
		for(i = 0; i < ni + 1; ++i)
		{
			if( i == 0 || i == ni )
			{
				pFaces[l_faceCount].m_b_IfBoundary = true;
			}
			else	
			{
				pFaces[l_faceCount].m_b_IfBoundary = false;
			}
			l_faceCount++;
		}
	}*/
	
	l_CellCount = 0;
	l_faceCount = 0; 
	for(j = 0; j < nj; ++j)
	{
		for(i = 0; i < ni; ++i)
		{
			if(i == 0)
			{
				pCells[l_CellCount].m_l_FaceIndex_West 	= l_faceCount;
				pCells[l_CellCount].m_l_FaceIndex_East 	= l_faceCount+1;
				if(j == 0)
				{
					pCells[l_CellCount].m_l_FaceIndex_South 	= l_nbFaces/2;
					pCells[l_CellCount].m_l_FaceIndex_North 	= (l_nbFaces/2) + 1;
				}
				else
				{
					pCells[l_CellCount].m_l_FaceIndex_South 	= pCells[l_CellCount - ni ].m_l_FaceIndex_North;
					pCells[l_CellCount].m_l_FaceIndex_North 	= pCells[l_CellCount - ni ].m_l_FaceIndex_North + 1;
				}
			}
			else
			{
				pCells[l_CellCount].m_l_FaceIndex_West 		= pCells[l_CellCount - 1].m_l_FaceIndex_East;
				pCells[l_CellCount].m_l_FaceIndex_East 		= pCells[l_CellCount - 1].m_l_FaceIndex_East+1;
				pCells[l_CellCount].m_l_FaceIndex_South 	= pCells[l_CellCount - 1].m_l_FaceIndex_South + ni + 1;
				pCells[l_CellCount].m_l_FaceIndex_North 	= pCells[l_CellCount - 1].m_l_FaceIndex_North + ni + 1;
				
			}	
			pCells[l_CellCount].m_l_nbRays = l_nbRays;
			if(pCells[l_CellCount].m_d_vec_cell_ray_Inten.size())
			{
				pCells[l_CellCount].m_d_vec_cell_ray_Inten.clear();
			}
			pCells[l_CellCount].m_d_vec_cell_ray_Inten.resize(l_nbRays);
			for(l_RayCt = 0; l_RayCt < l_nbRays; ++l_RayCt)
			{
				//pCells[l_CellCount].m_d_vec_cell_ray_Inten[l_RayCt] = dSBConst*pow(fTemp,4)/pi;
				pCells[l_CellCount].m_d_vec_cell_ray_Inten[l_RayCt] = 0.0;
			}
			pCells[l_CellCount].m_f_Vol = ((1/(float)ni)*(1/(float)nj));
			pCells[l_CellCount].m_f_Temp = fTemp;
			++l_CellCount;	
		}
		l_faceCount += ni + 1;
	}
	/*for(int i = 0; i < l_cellSize; ++i)
	{
		std::vector<long>l_vec_FaceIndicies;
		l_vec_FaceIndicies.push_back(i); // first one as south face
		l_vec_FaceIndicies.push_back(i+1); // second one as east face
		l_vec_FaceIndicies.push_back(i+2); // third one as north face
		l_vec_FaceIndicies.push_back(i+3); // fourth one as west face
		pCell[i].
	}*/
	
	
}

// ************************************************************************* //
