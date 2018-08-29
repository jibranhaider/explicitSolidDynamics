/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::angularMomentum

Description
    Angular momentum preservation algorithm

SourceFiles
    angularMomentum.C


\*---------------------------------------------------------------------------*/

#ifndef angularMomentum_H
#define angularMomentum_H

#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class angularMomentum
{
    // Private data

    //- mesh
         const fvMesh& mesh_;  

         dimensionedScalar rho_;       

    //- Disallow default bitwise copy construct
        angularMomentum(const angularMomentum&);

    //- Disallow default bitwise assignment
        void operator=(const angularMomentum&);

public:

    // Declare name of the class and its debug switch
        TypeName("angularMomentum");

    // Constructors
        explicit angularMomentum(const fvMesh&, const dimensionedScalar&);

    // Destructor
        virtual ~angularMomentum();

    // Member functions

		// Angular momentum conservation
        void AMconservation
        (
			GeometricField<vector, fvPatchField, volMesh>& ,
            GeometricField<vector, fvPatchField, volMesh>& ,                                      
            const GeometricField<vector, fvPatchField, volMesh>& ,            
            const word& 
        ) const;


        // Print global momentum
        void printGlobalMomentum
        (
            GeometricField<vector, fvPatchField, volMesh>& ,
            GeometricField<vector, fvPatchField, volMesh>&                        
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
