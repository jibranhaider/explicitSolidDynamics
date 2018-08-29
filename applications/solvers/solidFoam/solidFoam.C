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

Application
    solidFoam

Description
    A solid mechanics solver based on the cell centred Finite Volume Method 
    to predict large strain deformation. The Total Lagrangian mixed formuation 
	comprises of conservation laws for linear momentum and deformation gradient 
	of the system. It utilises an explicit Total Variation Diminishing two-stage 
    Runge-Kutta time integrator. A discrete angular momentum projection 
	algorithm based on two global Lagrange Multipliers guarantees angular 
	momentum conservation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "interpolationSchemes.H"
#include "angularMomentum.H"
#include "constitutiveModel.H"
#include "gradientSchemes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{ 
	#include "setRootCase.H"	
	#include "createTime.H"
	#include "createMesh.H"	 
	#include "readControls.H"	
	#include "meshData.H" 						
	#include "createFields.H"												


	while (runTime.loop()) 
	{

		if (timeStepping == "variable")
		{
			deltaT = ( cfl*h ) / max(Up_time);
			runTime.setDeltaT(deltaT);
		}		

		t += deltaT; tstep++;	

		Info << "\nTime Step =" << tstep << "\ndeltaT = " << deltaT.value() << " s"
			 << "\nTime = " << t.value() << " s" << endl;


		RKstage = "first";	
		#include "gEqns.H"
													
		RKstage = "second"; 
		#include "updateVariables.H"	
		#include "gEqns.H"

		lm	= 0.5 * ( lm.oldTime() + lm );
		F 	= 0.5 * ( F.oldTime() + F );							
		x 	= 0.5 * ( x.oldTime() + x );
		xF 	= 0.5 * ( xF.oldTime() + xF );			
		xN 	= 0.5 * ( xN.oldTime() + xN );			

		#include "updateVariables.H"	
		uN = xN - xN_0;				
		
		if (runTime.outputTime())
		{
			p.write();
			uN.write();
		}

		Info << "Percentage completed = " 
			 << (t.value()/runTime.endTime().value())*100 << "%" << endl;	
	}

	Info << "\nExecutionTime = " << runTime.elapsedCpuTime() << " s" 
		 << "  ClockTime = " << runTime.elapsedClockTime() << " s" 
	 	 << "\nEnd\n" << endl;
	
	return 0;
}

// ************************************************************************* //



