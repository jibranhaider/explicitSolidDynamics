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

#include "symmetricTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tractionValue_(vector::zero)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    tractionValue_(vector::zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    tractionValue_ =
        dict.lookupOrDefault<vector>("tractionValue", vector::zero);

    updateCoeffs();
}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tractionValue_(ptf.tractionValue_)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    tractionValue_(rifvpvf.tractionValue_)
{}


symmetricTractionFvPatchVectorField::symmetricTractionFvPatchVectorField
(
    const symmetricTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),
    tractionValue_(rifvpvf.tractionValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void symmetricTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


void symmetricTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void symmetricTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("lm_M");

    const fvsPatchField<vector>& t_M_ =
        patch().lookupPatchField<surfaceVectorField, vector>("t_M");

    const fvsPatchField<tensor>& nCn_ =
        patch().lookupPatchField<surfaceTensorField, tensor>("nCn");

    const fvsPatchField<tensor>& iMnCn_ =
        patch().lookupPatchField<surfaceTensorField, tensor>("iMnCn");

    const fvPatchField<scalar>& Up_ =
        patch().lookupPatchField<volScalarField, scalar>("Up");

    fvsPatchField<vector> t_C(lm_M_);
    t_C = (nCn_ & (t_M_ - Up_*lm_M_)) + (iMnCn_ & tractionValue_);

    this->operator==(t_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void symmetricTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("tractionValue") << tractionValue_ << token::END_STATEMENT
        << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    symmetricTractionFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //