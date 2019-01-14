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
    Foam::gradientSchemes

Description
    Define various gradient computaion schemes

SourceFiles
    gradientSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gradientSchemes_H
#define gradientSchemes_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;


class gradientSchemes
{
    // Private data

        // Mesh
        const fvMesh& mesh_;
        
        // Owner list
        const labelUList& own_;
        
        // Neighbour list
        const labelUList& nei_;
        
        // Cell centre coordinates
        const volVectorField& X_; 

        // Inverse of distance matrix                           
        volTensorField Ainv_;


private:

        // Disallow default bitwise copy construct
        gradientSchemes(const gradientSchemes&);

        // Disallow default bitwise assignment
        void operator=(const gradientSchemes&);


public:

    // Declare name of the class and its debug switch
    TypeName("gradientSchemes");

    // Constructors
    gradientSchemes(const fvMesh&);


    // Destructor
    virtual ~gradientSchemes();


    // Member functions

        // Inverse distance matrix for gradient
        void distanceMatrix
        (
            GeometricField<tensor, fvPatchField, volMesh>& U
        );

        // Local inverse distance matrix calculation for local gradient
        void distanceMatrixLocal
        (
            GeometricField<tensor, fvPatchField, volMesh>&
        ) const;


        // Least square gradient for a volScalarField
        void gradient
        (
            const GeometricField<scalar, fvPatchField, volMesh>& U,
            GeometricField<vector, fvPatchField, volMesh>& Ugrad
        ) const;

        // Least square gradient for a volVectorField
        void gradient
        (
            const GeometricField<vector, fvPatchField, volMesh>& U,
            GeometricField<tensor, fvPatchField, volMesh>& Ugrad
        ) const;

        // Least square gradient for a volTensorField
        void gradient
        (
            const GeometricField<tensor, fvPatchField, volMesh>& U,
            GeometricField<tensor, fvPatchField, volMesh>& UgradX,
            GeometricField<tensor, fvPatchField, volMesh>& UgradY,
            GeometricField<tensor, fvPatchField, volMesh>& UgradZ
        ) const;

        // Local least square gradient for a volVectorField
        volTensorField localGradient 
        (
            const GeometricField<vector, fvPatchField, volMesh>& U,
            const GeometricField<vector, fvsPatchField, surfaceMesh>& Unei,            
            const GeometricField<tensor, fvPatchField, volMesh>& Ainv                    
        ) const;


        // Reconstruction of a volScalarField to surfaceScalarField
        void reconstruct
        (
            GeometricField<scalar, fvPatchField, volMesh>& U,
            const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
            GeometricField<scalar, fvsPatchField, surfaceMesh>& Uminus,
            GeometricField<scalar, fvsPatchField, surfaceMesh>& Uplus
        );  

        // Reconstruction of a volVectorField to surfaceVectorField
        void reconstruct
        (
            GeometricField<vector, fvPatchField, volMesh>& U,
            const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
            GeometricField<vector, fvsPatchField, surfaceMesh>& Uminus,
            GeometricField<vector, fvsPatchField, surfaceMesh>& Uplus
        );      

        // Reconstruction for a volTensorField to surfaceTensorField
        void reconstruct
        (
            GeometricField<tensor, fvPatchField, volMesh>& U,
            const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
            const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
            const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
            GeometricField<tensor, fvsPatchField, surfaceMesh>& Uminus,
            GeometricField<tensor, fvsPatchField, surfaceMesh>& Uplus
        );  

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //