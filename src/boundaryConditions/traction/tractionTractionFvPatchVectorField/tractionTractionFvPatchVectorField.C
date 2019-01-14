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

#include "tractionTractionFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionTractionFvPatchVectorField::
tractionTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    loadingType_("none"),
    traction_(vector::zero),
    pressure_(0.0)
{}


tractionTractionFvPatchVectorField::
tractionTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    loadingType_(dict.lookupOrDefault<word>("loadingType", "none")),
    traction_(dict.lookupOrDefault("traction", vector::zero)),
    pressure_(dict.lookupOrDefault("pressure", 0.0))  
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (loadingType_ == "traction" && !dict.found("traction") == 1)
    {
        FatalErrorIn ("tractionTractionFvPatchVectorField")   
        << "Keyword 'traction' is undefined in dictionary '" 
        << iF.name() << ".boundaryField." << p.name() << "' for patch type '" << this->type() << "'." << nl 
        << exit(FatalError);
    }

    if (loadingType_ == "pressure" && !dict.found("pressure") == 1)
    {
        FatalErrorIn ("tractionTractionFvPatchVectorField")   
        << "Keyword 'pressure' is undefined in dictionary '" 
        << iF.name() << ".boundaryField." << p.name() << "' for patch type '" << this->type() << "'." << nl 
        << exit(FatalError);
    }

    updateCoeffs();    
}


tractionTractionFvPatchVectorField::
tractionTractionFvPatchVectorField
(
    const tractionTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    loadingType_(ptf.loadingType_),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_)
{}


tractionTractionFvPatchVectorField::
tractionTractionFvPatchVectorField
(
    const tractionTractionFvPatchVectorField& rifvpvf
)
:
    fixedValueFvPatchVectorField(rifvpvf),
    loadingType_(rifvpvf.loadingType_),    
    traction_(rifvpvf.traction_),
    pressure_(rifvpvf.pressure_)  
{}


tractionTractionFvPatchVectorField::
tractionTractionFvPatchVectorField
(
    const tractionTractionFvPatchVectorField& rifvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rifvpvf, iF),   
    loadingType_(rifvpvf.loadingType_),    
    traction_(rifvpvf.traction_),
    pressure_(rifvpvf.pressure_)     
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);   
}


void tractionTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);  
}


void tractionTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<vector>& lm_M_ = patch().lookupPatchField<surfaceVectorField, vector>("lm_M");
    
    fvsPatchField<vector> t_C(lm_M_);
    

    if (loadingType_ == "none")
    {    
        t_C = vector::zero;
    }

    else if (loadingType_ == "traction")
    {    
        t_C = traction_;     
    }

    else if (loadingType_ == "pressure")
    {
        const fvsPatchField<vector>& n_ = patch().lookupPatchField<surfaceVectorField, vector>("n"); 
        t_C = -pressure_ * n_;     
    }

    else
    {
        FatalErrorIn ("tractionTractionFvPatchVectorField::updateCoeffs()")   
        << "Unknown traction type '" << loadingType_ << "'" << nl
        << "Valid types are: 'none', 'traction' and 'pressure'" << nl
        << exit(FatalError);
    }

    this->operator==(t_C);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void tractionTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("loadingType") << loadingType_ << token::END_STATEMENT << nl;
    os.writeKeyword("traction") << traction_ << token::END_STATEMENT << nl;
    os.writeKeyword("pressure") << pressure_ << token::END_STATEMENT << nl;
    writeEntry("value", os);     
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionTractionFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //