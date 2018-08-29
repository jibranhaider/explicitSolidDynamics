

rhsLm = fvc::surfaceIntegrate(tC*magSf);

if (angularMomentumConservation == "yes")
{	
	rhsAm = fvc::surfaceIntegrate((xF^tC)*magSf);	
	AM.AMconservation( rhsLm, rhsLm1, rhsAm, RKstage);
}


// Update coordinates
x += deltaT * (lm/rho);
xF += deltaT * (lmC/rho);						
xN += deltaT * (lmN/rho);


// Update linear momentum
lm += deltaT * rhsLm;	


// Update deformation gradient tensor
F += deltaT * fvc::surfaceIntegrate( (lmC/rho) * Sf );	


