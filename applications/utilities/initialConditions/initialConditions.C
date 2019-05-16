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

Application
    initialConditions

Description
    Generates non-standard initial conditions for test cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Read mechanical properties dictionary
    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read run parameters dictionary
    IOdictionary runParameters
    (
        IOobject
        (
            "runParameters",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Test case name
    const word& tutorial = runParameters.lookup("tutorial");

    // Read density
    const dimensionedScalar& rho = mechanicalProperties.lookup("rho");

    // Cell centre coordinates
    const volVectorField& C = mesh.C();

    // Read linear momentum field
    volVectorField lm
    (
        IOobject
        (
            "lm",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
    );

    // Compute linear momentum
    if (tutorial == "twistingColumn" || tutorial == "spinningTop")
    {
        // Read initial angular velocity
        volVectorField omega
        (
            IOobject("omega", mesh),
            mesh,
            dimensionedVector(runParameters.lookup("initialAngularVelocity"))
        );

        if (tutorial == "twistingColumn")
        {
            const scalar& PI = Foam::constant::mathematical::pi;
            const dimensionedScalar height
            (
                "height",
                dimensionSet(0,1,0,0,0,0,0), 6.0
            );

            omega *= Foam::sin(PI*C.component(1)/(2*height));
        }

        lm = rho*omega ^ C;
    }

    else if (tutorial == "taylorImpact")
    {
        // Read initial velocity
        volVectorField v
        (
            IOobject("v", mesh),
            mesh,
            dimensionedVector(runParameters.lookup("initialVelocity"))
        );

        lm = rho*v;
    }

    else
    {
        FatalErrorIn("initialConditions.C")
            << "Valid type entries are 'twistingColumn', 'spinningTop'"
            << "or 'taylorImpact' for tutorial"
            << abort(FatalError);
    }

    lm.write();

    Info<< "\n end\n";

    return 0;
}

// ************************************************************************* //