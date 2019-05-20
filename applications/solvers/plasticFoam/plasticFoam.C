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
    mixedPlasticFoam

Description
    A large strain elasto-plastic solid mechanics solver based on a linear
    momentum/strains mixed formulation. An explicit Total Lagrangian
    formulation utilisiing a monolithic Total Variation Diminishing
    Runge-Kutta time integrator.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "gradientSchemes.H"
#include "interpolationSchemes.H"
#include "angularMomentum.H"
#include "mechanics.H"
#include "operations.H"
#include "plasticityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readControls.H"
    #include "createFields.H"

    while (runTime.loop())
    {
        if (timeStepping == "variable")
        {
            deltaT = (cfl*h)/max(Up_time);
            runTime.setDeltaT(deltaT);
        }

        t += deltaT; tstep++;

        Info<< "\nTime Step =" << tstep << "\ndeltaT = " << deltaT.value()
            << " s" << nl << "Time = " << t.value() << " s" << endl;

        lm.oldTime();
        F.oldTime();
        x.oldTime();
        xF.oldTime();
        xN.oldTime();

        forAll(RKstages, stage)
        {
            #include "gEqns.H"

            if (RKstages[stage] == 0)
            {
                #include "updateVariables.H"
            }
        }

        lm = 0.5*(lm.oldTime() + lm);
        F = 0.5*(F.oldTime() + F);
        x = 0.5*(x.oldTime() + x);
        xF = 0.5*(xF.oldTime() + xF);
        xN = 0.5*(xN.oldTime() + xN);

        #include "updateVariables.H"

        if (runTime.outputTime())
        {
            uN = xN - XN;
            uN.write();
            model.writeOutput();
        }

        Info<< "Percent completed = "
            << (t.value()/runTime.endTime().value())*100 << "%" << endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    return 0;
}

// ************************************************************************* //