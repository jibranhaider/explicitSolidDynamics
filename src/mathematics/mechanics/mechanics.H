/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::mechanics

Description
    Evaluate quantities widely used in continuum mechanics.

SourceFiles
    mechanics.C

\*---------------------------------------------------------------------------*/

#ifndef mechanics_H
#define mechanics_H

#include "fvc.H"
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mechanics
{
    // Private data

        //- Mesh
         const fvMesh& mesh_;

        // //- Eigen-values
        // vector  eigVal_;

        // //- Eigen-vectors
        // tensor  eigVec_;

private:

    //- Disallow default bitwise copy construct
        mechanics(const mechanics&);

    //- Disallow default bitwise assignment
        void operator=(const mechanics&);

public:

    //- Runtime type information
    TypeName("mechanics");

    // Constructors

        //- Construct from mesh and dictionary
        mechanics
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~mechanics();


    // Member functions

        // Edit

            //- Compute stretch
            volScalarField stretch();

            //- Computes spatial surface normals
            surfaceVectorField spatialNormal();

        // Output

            //- Print centroid of geometry
            void printCentroid () const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //