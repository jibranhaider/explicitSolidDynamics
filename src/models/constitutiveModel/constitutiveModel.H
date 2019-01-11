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
    Foam::constitutiveModel

Description
    Specifies the consitututive model

SourceFiles
    constitutiveModel.C


\*---------------------------------------------------------------------------*/

#ifndef constitutiveModel_H
#define constitutiveModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;


class constitutiveModel
{
    // Private data

        // First Piola Kirchhoff stress tensor 
        volTensorField P_;

        // Pressure 
        volScalarField p_; 

        // Energy of algorithm 
        volScalarField Ealgo_; 

        // Density
        dimensionedScalar rho_; 

        // Young's Modulus of Elasticity                     
        dimensionedScalar E_;

        // Poisson's ratio
        dimensionedScalar nu_;

        // Lame coefficient        
        dimensionedScalar mu_;

        // Lame coefficient
        dimensionedScalar lambda_;
        
        // Bulk modulus 
        dimensionedScalar kappa_;

        // Constitutive model
        word model_;  

        // Pressure wave speed
        dimensionedScalar Up_;

        // Shear wave speed
        dimensionedScalar Us_;

        // Eigen-values
        vector  eigVal_;

        // Eigen-vectors
        tensor  eigVec_;        

private:

    //- Disallow default bitwise copy construct
        constitutiveModel(const constitutiveModel&);

    //- Disallow default bitwise assignment
        void operator=(const constitutiveModel&);

public:

    // Declare name of the class and its debug switch
        TypeName("constitutiveModel");

    // Constructors
        constitutiveModel 
        (
            const volTensorField&, 
            const dictionary& dict
        );


    // Destructor
        virtual ~constitutiveModel();

    // Member functions
         
        const volTensorField& piola() const  
        {
            return P_;  
        }


        const volScalarField& pressure() const  
        {
            return p_;   
        }

        const dimensionedScalar& density() const  
        {
            return rho_;   
        }

        const dimensionedScalar& lambda() const  
        {
            return lambda_;   
        }

        const dimensionedScalar& mu() const  
        {
            return mu_;   
        }

        const dimensionedScalar& kappa() const  
        {
            return kappa_;   
        }

        const dimensionedScalar& Up() const  
        {
            return Up_;   
        }

        const dimensionedScalar& Us() const  
        {
            return Us_;   
        }


        // Non-linear pressure wave speed
        volScalarField Up_nonLinear();  

        // Computes eigen-values and eigen-vectors
        void eigenStructure(const tensor&);

        // Computes spatial surface normals
        surfaceVectorField spatialNormal();


        virtual void correct();

};




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
