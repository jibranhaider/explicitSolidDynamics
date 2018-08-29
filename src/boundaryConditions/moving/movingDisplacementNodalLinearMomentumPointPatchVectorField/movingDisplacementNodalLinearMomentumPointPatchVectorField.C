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

\*---------------------------------------------------------------------------*/

#include "movingDisplacementNodalLinearMomentumPointPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    density_(0.0),
    displacement_(vector::zero),
    endTime_(0.0)
{}


movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
    density_( readScalar(dict.lookup("density")) ),  
    displacement_( dict.lookup("displacement") ),  
    endTime_( readScalar(dict.lookup("endTime")) ) 
{
    pointPatchVectorField::operator=(vectorField("value", dict, p.size()));

	updateCoeffs();
}


movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const movingDisplacementNodalLinearMomentumPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,   
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ptf, p, iF, mapper),
    density_(ptf.density_), 
    displacement_(ptf.displacement_),
    endTime_(ptf.endTime_)    
{}


movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const movingDisplacementNodalLinearMomentumPointPatchVectorField& rifvpvf
)
:
    fixedValuePointPatchVectorField(rifvpvf),
    density_(rifvpvf.density_),    
    displacement_(rifvpvf.displacement_),
    endTime_(rifvpvf.endTime_)   
{}


movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const movingDisplacementNodalLinearMomentumPointPatchVectorField& rifvpvf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(rifvpvf, iF),
    density_(rifvpvf.density_),
    displacement_(rifvpvf.displacement_),
    endTime_(rifvpvf.endTime_) 
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingDisplacementNodalLinearMomentumPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


void movingDisplacementNodalLinearMomentumPointPatchVectorField::rmap
(
    const pointPatchVectorField& ptf,    
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}


void movingDisplacementNodalLinearMomentumPointPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField lm_N(this->patchInternalField());
         
    const scalar& t = this->db().time().value();
    const scalar A = 3*10;
    const scalar B = 4*15 / endTime_;
    const scalar C = 5*6 / pow(endTime_,2);

    lm_N = density_ * (displacement_/pow(endTime_,3)) * (A - B*t + C*t*t) * pow(t,2);

    this->operator==(lm_N);   
    fixedValuePointPatchVectorField::updateCoeffs();
}


void movingDisplacementNodalLinearMomentumPointPatchVectorField::write(Ostream& os) const
{
    pointPatchVectorField::write(os);
    os.writeKeyword("density") << density_ << token::END_STATEMENT << nl;
    os.writeKeyword("displacement") << displacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("endTime") << endTime_ << token::END_STATEMENT << nl;     
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    movingDisplacementNodalLinearMomentumPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //