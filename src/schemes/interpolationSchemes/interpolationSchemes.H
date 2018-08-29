/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::interpolationSchemes

Description
    Define various interpolation schemes

SourceFiles
    interpolationSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationSchemes_H
#define interpolationSchemes_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "pointConstraints.H"

#include "volPointInterpolation.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;


class interpolationSchemes
{
    // Private data

    //- mesh
         const fvMesh& mesh_;         

    //- Disallow default bitwise copy construct
        interpolationSchemes(const interpolationSchemes&);

    //- Disallow default bitwise assignment
        void operator=(const interpolationSchemes&);

public:

    // Declare name of the class and its debug switch
        TypeName("interpolationSchemes");

    // Constructors
        explicit interpolationSchemes(const fvMesh&);

    // Destructor
        virtual ~interpolationSchemes();

    // Member functions


		// Surface to cell interpolation 
        volVectorField surfaceToCell  
        (
			const GeometricField<vector, fvsPatchField, surfaceMesh>&		
        ) const;


        template<class Type>
        void pushUntransformedDataNew(List<Type>&) const;


        template<class Type>
        void addSeparatedNew
        (
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;


        // Volume to point interpolation for vectorField
        void volToPoint
        (
            const GeometricField<vector, fvPatchField, volMesh>& U,
            const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
            GeometricField<vector, pointPatchField, pointMesh>& Un                               
        ) const;


        // Point to surface interpolation for vectorField
        void pointToSurface   
        (
            const GeometricField<vector, pointPatchField, pointMesh>& U,
            GeometricField<vector, fvsPatchField, surfaceMesh>& Uf                   
        ) const;


        // Constrained surface to point interpolation for vectorField
        void surfaceToPointConstrained
        (
            const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf,
            GeometricField<vector, pointPatchField, pointMesh>& Un                  
        ) const;
};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
