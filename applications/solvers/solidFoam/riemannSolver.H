//-------------------------
// ACOUSTIC RIEMANN SOLVER 
//-------------------------	


lmC = 0.5*(lm_M+lm_P) + ( 0.5 * ((nCn/fvc::interpolate(Up)) + (iMnCn/fvc::interpolate(Us))) & (t_P-t_M) );	
tC  = 0.5*(t_M+t_P) + ( 0.5 * ((fvc::interpolate(Up)*nCn) + (fvc::interpolate(Us)*iMnCn)) & (lm_P-lm_M) );


if( Pstream::parRun() ) 
{
	forAll (mesh.cells(), cellID)
	{
		Px[cellID] = vector(P[cellID].xx(), P[cellID].xy(), P[cellID].xz());
		Py[cellID] = vector(P[cellID].yx(), P[cellID].yy(), P[cellID].yz());
		Pz[cellID] = vector(P[cellID].zx(), P[cellID].zy(), P[cellID].zz());									
	}	

	Px.correctBoundaryConditions();
	Py.correctBoundaryConditions();	
	Pz.correctBoundaryConditions();	

	forAll(mesh.boundary(), patchID)
	{
		if (mesh.boundary()[patchID].coupled())				
		{			
			vectorField lm_nei = lm.boundaryField()[patchID].patchNeighbourField();

			tensorField P_nei = P.boundaryField()[patchID].patchNeighbourField();
			vectorField Px_nei = Px.boundaryField()[patchID].patchNeighbourField();
			vectorField Py_nei = Py.boundaryField()[patchID].patchNeighbourField();
			vectorField Pz_nei = Pz.boundaryField()[patchID].patchNeighbourField();						
			
			tensorField lmGrad_nei = lmGrad.boundaryField()[patchID].patchNeighbourField();
			tensorField PxGrad_nei = PxGrad.boundaryField()[patchID].patchNeighbourField();
			tensorField PyGrad_nei = PyGrad.boundaryField()[patchID].patchNeighbourField();
			tensorField PzGrad_nei = PzGrad.boundaryField()[patchID].patchNeighbourField();			
		

			vectorField C_nei = C.boundaryField()[patchID].patchNeighbourField();			
			scalarField Up_nei = Up.boundaryField()[patchID].patchNeighbourField();	
			scalarField Us_nei = Us.boundaryField()[patchID].patchNeighbourField();		
			
			
			forAll (mesh.boundary()[patchID],facei) 
			{
				const label& bCellID = mesh.boundaryMesh()[patchID].faceCells()[facei];				
				const vector& Cf_ = Cf.boundaryField()[patchID][facei];

				const vector& lm_M_ = lm[bCellID] + ( lmGrad[bCellID] & (Cf_ - C[bCellID]) );
				const vector& lm_P_ = lm_nei[facei] + ( lmGrad_nei[facei] & (Cf_ - C_nei[facei]) );					

				const vector& Px_M_ = Px[bCellID] 	+ ( PxGrad[bCellID] & (Cf_ - C[bCellID]) );
				const vector& Px_P_ = Px_nei[facei] + ( PxGrad_nei[facei] & (Cf_ - C_nei[facei]) );		
				const vector& Py_M_ = Py[bCellID] 	+ ( PyGrad[bCellID] & (Cf_ - C[bCellID]) );
				const vector& Py_P_ = Py_nei[facei] + ( PyGrad_nei[facei] & (Cf_ - C_nei[facei]) );																			
				const vector& Pz_M_ = Pz[bCellID] 	+ ( PzGrad[bCellID] & (Cf_ - C[bCellID]) );
				const vector& Pz_P_ = Pz_nei[facei] + ( PzGrad_nei[facei] & (Cf_ - C_nei[facei]) );	


				const tensor& P_M_ = tensor( Px_M_, Py_M_, Pz_M_ );	
				const tensor& P_P_ = tensor( Px_P_, Py_P_, Pz_P_ );
					
			
				const scalar& Up_ = (Up[bCellID] + Up_nei[facei]) / 2.0;
				const scalar& Us_ = (Us[bCellID] + Us_nei[facei]) / 2.0;	
				
				const vector& N_ = N.boundaryField()[patchID][facei];								
				const tensor& nCn_ = nCn.boundaryField()[patchID][facei];
				const tensor& iMnCn_ = iMnCn.boundaryField()[patchID][facei];						
				
	
				lmC.boundaryFieldRef()[patchID][facei] = 0.5*(lm_M_+lm_P_) + 0.5 * ( ((nCn_/Up_) + (iMnCn_/Us_)) & ((P_P_-P_M_) & N_) );	
				tC.boundaryFieldRef()[patchID][facei] = 0.5*((P_M_+P_P_) & N_) + ( 0.5 * ((Up_*nCn_) + (Us_*iMnCn_)) & (lm_P_-lm_M_) );
			}					
		}		
	}
}




lm_b.correctBoundaryConditions();
t_b.correctBoundaryConditions();

forAll(mesh.boundary(), patchID) 
{
	if (mesh.boundary()[patchID].coupled())				
	{
	}

	else
	{
		forAll (mesh.boundary()[patchID],facei) 
		{
			lmC.boundaryFieldRef()[patchID][facei] = lm_b.boundaryField()[patchID][facei];
			tC.boundaryFieldRef()[patchID][facei] = t_b.boundaryField()[patchID][facei];
		}
	}
}


