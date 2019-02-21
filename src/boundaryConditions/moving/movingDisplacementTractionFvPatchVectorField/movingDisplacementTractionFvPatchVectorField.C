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

#include "movingDisplacementTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    movingTractionFvPatchVectorField(p, iF),
    density_(0.0),
    displacement_(vector::zero),
    endTime_(0.0)
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    movingTractionFvPatchVectorField(p, iF),
    density_(readScalar(dict.lookup("density"))),
    displacement_(dict.lookup("displacement")),
    endTime_(readScalar(dict.lookup("endTime")))
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const movingDisplacementTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    movingTractionFvPatchVectorField(ptf, p, iF, mapper),
    density_(ptf.density_),
    displacement_(ptf.displacement_),
    endTime_(ptf.endTime_)
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const movingDisplacementTractionFvPatchVectorField& rifvpvf
)
:
    movingTractionFvPatchVectorField(rifvpvf),
    density_(rifvpvf.density_),
    displacement_(rifvpvf.displacement_),
    endTime_(rifvpvf.endTime_)
{}


movingDisplacementTractionFvPatchVectorField::
movingDisplacementTractionFvPatchVectorField
(
    const movingDisplacementTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    movingTractionFvPatchVectorField(rifvpvf, iF),
    density_(rifvpvf.density_),
    displacement_(rifvpvf.displacement_),
    endTime_(rifvpvf.endTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingDisplacementTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    movingTractionFvPatchVectorField::autoMap(m);
}


void movingDisplacementTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    movingTractionFvPatchVectorField::rmap(ptf, addr);
}


void movingDisplacementTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar& t = this->db().time().value();
    const scalar A = 3*10;
    const scalar B = 4*15/endTime_;
    const scalar C = 5*6/pow(endTime_,2);

    linearMomentum() =
        density_*(displacement_/pow(endTime_,3))*(A - B*t + C*t*t)*pow(t,2);

    movingTractionFvPatchVectorField::updateCoeffs();
}


void movingDisplacementTractionFvPatchVectorField::write(Ostream& os) const
{
    movingTractionFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingDisplacementTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //