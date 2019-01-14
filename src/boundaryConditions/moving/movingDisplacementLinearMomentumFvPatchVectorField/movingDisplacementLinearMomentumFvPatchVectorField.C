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

\*---------------------------------------------------------------------------*/

#include "movingDisplacementLinearMomentumFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingDisplacementLinearMomentumFvPatchVectorField::
movingDisplacementLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    density_(0.0),
    displacement_(vector::zero),
	endTime_(0.0)
{}


movingDisplacementLinearMomentumFvPatchVectorField::
movingDisplacementLinearMomentumFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    density_( readScalar(dict.lookup("density")) ),  
    displacement_( dict.lookup("displacement") ), 
    endTime_( readScalar(dict.lookup("endTime")) )   
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

	updateCoeffs();
}


movingDisplacementLinearMomentumFvPatchVectorField::
movingDisplacementLinearMomentumFvPatchVectorField
(
    const movingDisplacementLinearMomentumFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    density_(ptf.density_), 
    displacement_(ptf.displacement_),
    endTime_(ptf.endTime_)    
{}


movingDisplacementLinearMomentumFvPatchVectorField::movingDisplacementLinearMomentumFvPatchVectorField
(
    const movingDisplacementLinearMomentumFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    density_(rifvpvf.density_),    
    displacement_(rifvpvf.displacement_),
    endTime_(rifvpvf.endTime_)    
{}


movingDisplacementLinearMomentumFvPatchVectorField::movingDisplacementLinearMomentumFvPatchVectorField
(
    const movingDisplacementLinearMomentumFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    density_(rifvpvf.density_),
    displacement_(rifvpvf.displacement_),
    endTime_(rifvpvf.endTime_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingDisplacementLinearMomentumFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


void movingDisplacementLinearMomentumFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void movingDisplacementLinearMomentumFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ = patch().lookupPatchField<surfaceVectorField, vector>("lm_M");  

    fvsPatchField<vector> lm_C(lm_M_);
         
    const scalar& t = this->db().time().value();
    const scalar A = 3*10;
    const scalar B = 4*15 / endTime_;
    const scalar C = 5*6 / pow(endTime_,2);

    lm_C = density_ * (displacement_/pow(endTime_,3)) * (A - B*t + C*t*t) * pow(t,2);

    this->operator==(lm_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void movingDisplacementLinearMomentumFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("density") << density_ << token::END_STATEMENT << nl;
    os.writeKeyword("displacement") << displacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("endTime") << endTime_ << token::END_STATEMENT << nl;           
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingDisplacementLinearMomentumFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //