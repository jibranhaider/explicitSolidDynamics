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
    Foam::movingDisplacementNodalLinearMomentumPointPatchVectorField

Description
    Prescribes an inlet profile for a circular pipe inlet patch. Flow direction is determined
    automatically from face normals, only velocity magnitude needs to be provided.
    \heading Patch usage
    \table
        Property    		| Description             								| Required  | Default value
        R           		| pipe radius 											| yes     	|
        flowSpeed   		| magnitude of flow velocity      						| yes       |
        deltaByR			| wall function thickness as a fraction of pipe radius  | yes 		|
        approximationType	| name of approximation profile							| no		| exponential
        centrePoint			| point on the centreline of the pipe					| no 		| (0 0 0)
        lambda				| pressure gradient coefficient 						| no 		| 0
    \endtable
    Example of the boundary condition specification:
    \verbatim
    inlet
	{
		// prescribes a parabolic inlet profile
	    type        prescribedPipeInlet;
	    // (optional, default exponential) type of profile used for approximation
	    approximationType exponential;
	    // pipe radius
	    R	$Rpipe;
	    // magnitude of flow velocity into the pipe
	    flowSpeed	$flowVelocity;
	    // (optional, defaults at origin) centrepoint of the pipe used to determine the radius of each face
	    centrepoint	(0 0 0);
	    // inlet BL thickness as a function of radius
	    deltaByR $deltaByR;
	    // (optional, default 0) pressure gradient coefficient (see Polhausen profile literature for explanation)
	    lambda 0.;
	    // dummy value overwritten by the BC, used for initialisation only
	    value		$internalField;
	}
	\endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    prescribedPipeInletPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingDisplacementNodalLinearMomentumPointPatchVectorField_H
#define movingDisplacementNodalLinearMomentumPointPatchVectorField_H

//#include "fvMesh.H"
#include "fvPatchFields.H"
#include "pointPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fvsPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "Switch.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class movingDisplacementNodalLinearMomentumPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingDisplacementNodalLinearMomentumPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{
private:
 
    // Density
    scalar density_;

    // Total displacement
    vector displacement_;

    // End time for displacement
    scalar endTime_;      

public:

    //- Runtime type information
    TypeName("movingDisplacementNodalLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&        
        );

        //- Construct from patch, internal field and dictionary
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&            
        );

        //- Construct by mapping given movingDisplacementNodalLinearMomentumPointPatchVectorField
        //  onto a new patch
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const movingDisplacementNodalLinearMomentumPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&            
        );

        //- Construct as copy
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const movingDisplacementNodalLinearMomentumPointPatchVectorField&      
        );

        virtual autoPtr<pointPatchVectorField> clone() const
        {
            return autoPtr<pointPatchVectorField>
            (
                new movingDisplacementNodalLinearMomentumPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const movingDisplacementNodalLinearMomentumPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new movingDisplacementNodalLinearMomentumPointPatchVectorField(*this,iF)
            );
        }


    // Member functions

        // Mapping functions

            // Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const pointPatchVectorField&,
                const labelList&
            );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //