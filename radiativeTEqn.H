/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    sphere

Description
    Spherical Geometry 

SourceFiles
    sphere.C

\*---------------------------------------------------------------------------*/

#ifndef radiativeTEqn_H
#define radiativeTEqn_H

#include <stdio.h>
#include <vector>
#include <map>
//#include "Map.H"
typedef int scalar;
//#include "absorptionEmissionModel.H"
//#include "blackBodyEmission.H"
typedef struct point
{
float x;
float y;
float z;
}Point;

class vertex
{
public:
float x;
float y;
float z;
long index;

};


class Cell
{
	
			long m_l_nbCells;
			long m_l_cellIndex;
	public:
			long m_l_FaceIndex_East;
			long m_l_FaceIndex_West;
			long m_l_FaceIndex_North;
			long m_l_FaceIndex_South;

			long m_l_CellIndex_East;
			long m_l_CellIndex_West;
			long m_l_CellIndex_North;
			long m_l_CellIndex_South;
			long m_l_nbRays;
			float m_f_Vol;
			float m_f_Temp;
			std::vector<long>m_l_vec_FaceIndicies;
			std::vector<double> m_d_vec_cell_ray_Inten;	

			// constructors
			//Cell();
			//Cell(long l_nbCells);

			// destructor
			//~Cell();
			
			
};
//class Face;
//class Cell;

class Face
{
	
			//long m_l_nbCells;
			long m_l_FaceIndex;
	public:
			
			long m_l_nbRays;
			bool m_b_IfBoundary;
			float m_f_Temp;
			float m_f_Area;
			//std::vector<long>m_l_vec_FaceIndicies;
			std::vector<double> m_d_vec_face_ray_Inten;	

			// constructors
			//Face();
			//Face(long l_nbCells);

			// destructor
			//~Face();
			
			
};


class radiativeTransferEquation
{

private: 
	
	int inbVerts_;

	int inbFaces_;

	int irecursionLevel_;

	Point v_[12];

	std::vector<Point>vec_TriangleNodes_;

	std::vector<Point>vec_TriangleMidNodes_;

	std::map<Point, long>mapMidPointNode_;

 	// Azimuthal Angles in PI/2 on x-y (from y to x)
	scalar nPhi_;

	// Polar Angle in PI(from Z to X-Y Plane)
	scalar nTheta_;

	// Radius of the Sphere
	float radius_;

	// centre of the Sphere
	Point centrePt_;

	std::map<Point, Point>mapCentroidDFm_;
	//Map<Point>mapCentroidDFm_;

	std::vector<std::pair<Point, Point> >vec_pair_CentroidDfm_;

	std::vector<Point> vec_Dfm_;

	std::vector<Point> vec_Cent_;

	std::vector<float> vec_omega_;

	//- Maximum omega weight
    scalar omegaMax_;

// Private Member Functions

        //- Disallow default bitwise copy construct
	radiativeTransferEquation(const radiativeTransferEquation&);
 	//- Disallow default bitwise assignment
        void operator=(const radiativeTransferEquation&);

	long createnewVertBasedonMiddlePoint(long lNode1, long lNode2);

	void normalize(Point &obj);

	float computeDistanceBetweenTwoPoints(Point a, Point b);

	bool isCollinear(Point obj1, Point obj2, Point obj3);

public:
// Constructors
//- Construct from components
        radiativeTransferEquation
        (
            const scalar nPhi,
            const scalar nTheta,
	    const float radius,
	    const Point centrePt
        );


    //- Destructor
    ~radiativeTransferEquation();

// Member functions

            //- Create Sphere 
            void createSphere
            (
		FILE *&fout
            );

	    void crossProduct(Point a, Point b, Point &cp);

		void dotProduct(Point a, Point b, float &dp);

	    void getSubtractionVector(Point a, Point b, Point &sv);

	    void computeCentroid(std::vector<Point>vecPoints, Point &centroid);

	    std::map<Point, Point>getMapCentroidDFm();

		void setBoundaryCondition(float **&T, float ***&I);

		void extractDataAndGetRayParameters(int &inbRays);

		void constructGridForControlVolume(float *&pfXCV, float *&pfYCV, float **&ppf_Vol);

		void constructGridUsingCellAndFaces(long l_nbRays, Face *&pFaces, Cell *&pCells, long &l_nbCells, long &l_nbFaces);

		void updateIntensity();
		
		void printFaceAndCellInformation(Face *pFaces, Cell *pCells,long l_nbCells, long l_nbFaces );
	    
	    //Map<Point>getMapCentroidDFm();

            //- Return the phi angle
            inline scalar nPhi() const;

            //- Return the theta angle
            inline scalar nTheta() const;

	    //- Return the radius
            inline float radius() const;

	    // Return Centre of the Sphere
	    inline Point pt() const;

	    inline std::vector<std::pair<Point, Point> >getVecPairCentDfm() const;

		


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "sphereI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
