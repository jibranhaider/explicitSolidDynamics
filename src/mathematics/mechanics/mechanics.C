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

#include "mechanics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanics, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

mechanics::mechanics
(
    const fvMesh& vm,
    const dictionary& dict
)
:
    mesh_(vm)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

mechanics::~mechanics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

surfaceVectorField mechanics::spatialNormal()
{
    const objectRegistry& db = mesh_.thisDb();
    const volTensorField& F_ = db.lookupObject<volTensorField>("F");

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject("n", mesh_),
            mesh_,
            dimensioned<vector>("n", dimless, pTraits<vector>::one)
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh> n = tsf();

    surfaceTensorField FcInv = inv(fvc::interpolate(F_));
    surfaceVectorField N = mesh_.Sf()/mesh_.magSf();

    n = (FcInv.T() & N)/(mag(FcInv.T() & N));
    tsf.clear();

    return n;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volScalarField mechanics::stretch()
{
    const objectRegistry& db = mesh_.thisDb();
    const volTensorField& F_ = db.lookupObject<volTensorField>("F");
    volTensorField C_ = F_.T() & F_;

    tmp<GeometricField<scalar, fvPatchField, volMesh> > tsf
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject("stretch", mesh_),
            mesh_,
            dimensioned<scalar>("stretch", dimless, pTraits<scalar>::one)
        )
    );
    GeometricField<scalar, fvPatchField, volMesh> stretch = tsf();

    operations op(mesh_);

    forAll(mesh_.cells(), cellID)
    {
        op.eigenStructure(C_[cellID]);
        vector eigVal_ = op.eigenValue();
        stretch[cellID] = min( eigVal_.x(), eigVal_.y());
        stretch[cellID] = ::sqrt(min(stretch[cellID], eigVal_.z()));
    }

    tsf.clear();

    if (Pstream::parRun())
    {
        stretch.correctBoundaryConditions();
    }

    return stretch;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mechanics::printCentroid() const
{
    vector sum = vector::zero;
    scalar vol = gSum(mesh_.V());

    forAll(mesh_.cells(), cellID)
    {
        sum += mesh_.C()[cellID]*mesh_.V()[cellID];
    }

    if( Pstream::parRun() )
    {
        reduce(sum, sumOp<vector>());
    }

    Info << "\nCentroid of geometry = " << sum/vol << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //