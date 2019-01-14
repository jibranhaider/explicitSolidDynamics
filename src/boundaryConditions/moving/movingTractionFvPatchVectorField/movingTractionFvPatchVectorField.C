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

#include "movingTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
	linearMomentum_(vector::zero)
{}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
	linearMomentum_(dict.lookupOrDefault<vector>("linearMomentum",vector::zero))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

	updateCoeffs();
}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const movingTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
	linearMomentum_(ptf.linearMomentum_)
{}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const movingTractionFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    linearMomentum_(rifvpvf.linearMomentum_)
{}


movingTractionFvPatchVectorField::movingTractionFvPatchVectorField
(
    const movingTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    linearMomentum_(rifvpvf.linearMomentum_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


void movingTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void movingTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ = patch().lookupPatchField<surfaceVectorField, vector>("lm_M");
    const fvsPatchField<vector>& t_M_ = patch().lookupPatchField<surfaceVectorField, vector>("t_M");    
    const fvsPatchField<tensor>& nCn_ = patch().lookupPatchField<surfaceTensorField, tensor>("nCn");
    const fvsPatchField<tensor>& iMnCn_ = patch().lookupPatchField<surfaceTensorField, tensor>("iMnCn");

    const fvPatchField<scalar>& Up_ = patch().lookupPatchField<volScalarField, scalar>("Up");
    const fvPatchField<scalar>& Us_ = patch().lookupPatchField<volScalarField, scalar>("Us");    

    fvsPatchField<vector> t_C(lm_M_);
 
    t_C = t_M_ + ( (Up_*nCn_ + Us_*iMnCn_) & (linearMomentum_ - lm_M_) );  

    this->operator==(t_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void movingTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("linearMomentum_") << linearMomentum_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //