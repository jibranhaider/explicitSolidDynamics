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

\*---------------------------------------------------------------------------*/

#include "gradientSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  
  
// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
  
gradientSchemes::gradientSchemes
(
	const fvMesh& vm
)
:
	mesh_(vm),

	own_(mesh_.owner()),
	nei_(mesh_.neighbour()),
	X_(mesh_.C()),

    Ainv_                              
    (
        IOobject("Ainv", mesh_),
        mesh_,
        dimensionedTensor("Ainv",dimensionSet(0,2,0,0,0,0,0),tensor::zero)    
    )

{
	gradientSchemes::distanceMatrix(Ainv_);
}

    
// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

  gradientSchemes::~gradientSchemes()
  {}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void gradientSchemes::distanceMatrix
(
    GeometricField<tensor, fvPatchField, volMesh>& U
)
{

	forAll (own_,faceID)
	{
		const label& ownCellID = own_[faceID];
		const label& neiCellID = nei_[faceID];
		const vector& dOwn = X_[neiCellID] - X_[ownCellID];
		const vector& dNei	= X_[ownCellID] - X_[neiCellID];	

		U[ownCellID] += dOwn*dOwn;
		U[neiCellID] += dNei*dNei;
	}

	if( Pstream::parRun() ) 
	{
		forAll(mesh_.boundary(), patchID)
		{
			if (mesh_.boundary()[patchID].coupled())
			{
				vectorField X_nei = X_.boundaryField()[patchID].patchNeighbourField();
				
				forAll (mesh_.boundary()[patchID],facei) 
				{
					const label& bCellID = mesh_.boundaryMesh()[patchID].faceCells()[facei];	
					const vector& d = X_nei[facei] - X_[bCellID];					
					U[bCellID] += d*d;			
				}		
			}		
		}

		U.correctBoundaryConditions();
	}


	U.primitiveFieldRef() = inv(U.internalField());
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrixLocal
(
	GeometricField<tensor, fvPatchField, volMesh>& Ainv
) const
{

	const volVectorField& C = mesh_.C();
	const surfaceVectorField& Cf = mesh_.Cf();
	const scalar& maxInteriorFaceID = mesh_.nInternalFaces()-1;

    const objectRegistry& db = mesh_.thisDb();
	const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

	tmp<GeometricField<tensor, fvPatchField, volMesh> > tvf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject("dummy", mesh_),
            mesh_,
            dimensioned<tensor>("dummy", Ainv.dimensions(), pTraits<tensor>::zero)
        )
    );

	GeometricField<tensor, fvPatchField, volMesh> dCd = tvf();


    forAll (mesh_.faces(), faceID)
    {
    	if (faceID <= maxInteriorFaceID)
    	{
			const label& ownID = own_[faceID];
			const label& neiID = nei_[faceID];

			const vector& dOwn = Cf[faceID] - X_[ownID];	
			const vector& dNei = Cf[faceID] - X_[neiID];	

			dCd[ownID] += dOwn * dOwn;
			dCd[neiID] += dNei * dNei;
    	}

    	else
    	{
			const label& patchID = mesh_.boundaryMesh().whichPatch(faceID);			
			const label& facei = mesh_.boundaryMesh()[patchID].whichFace(faceID);				
			const label& bCellID = mesh_.boundaryMesh()[patchID].faceCells()[facei];
			vector d = Cf.boundaryField()[patchID][facei] - C[bCellID];	 

			dCd[bCellID] += d * d;

			if (lmN_.boundaryField().types()[patchID] == "fixedValue")
			{	
	    		const pointVectorField& xN_0_ = db.lookupObject<pointVectorField> ("xN_0"); 	
				const label& faceID = mesh_.boundary()[patchID].start() + facei;

				forAll (mesh_.faces()[faceID], nodei)
				{
					const label& nodeID = mesh_.faces()[faceID][nodei];				
			
					d = xN_0_[nodeID] - C[bCellID];
					dCd[bCellID] += d * d;

					for (int i=0; i<7; i++)
					{										
						d = ( ( ((i+1)*xN_0_[nodeID]) + ((7-i)*Cf.boundaryField()[patchID][facei]) ) / 8.0 ) - C[bCellID];	
						dCd[bCellID] += d * d;								
					}
				}
			}		

    	}
    }

    Ainv.primitiveFieldRef() = inv(dCd.internalField());	
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
   	const GeometricField<scalar, fvPatchField, volMesh>& U,
   	GeometricField<vector, fvPatchField, volMesh>& Ugrad 
) 	const
{

	Ugrad = dimensionedVector("dummy", Ugrad.dimensions(), vector::zero);  

	forAll (own_, faceID)
	{
		const label& cellID = own_[faceID];
		const label& neiID = nei_[faceID];
		const vector& dcell = X_[neiID] - X_[cellID];	
		const vector& dnei = X_[cellID] - X_[neiID];	

		Ugrad[cellID] += Ainv_[cellID] & (U[neiID]-U[cellID]) * dcell;
		Ugrad[neiID] += Ainv_[neiID] & (U[cellID]-U[neiID]) * dnei;							
	}

	if( Pstream::parRun() ) 
	{
		forAll(mesh_.boundary(), patchID)
		{
			if (mesh_.boundary()[patchID].coupled())
			{
				vectorField X_nei = X_.boundaryField()[patchID].patchNeighbourField();
				scalarField U_nei = U.boundaryField()[patchID].patchNeighbourField();
					
				forAll (mesh_.boundary()[patchID],facei) 
				{
					const label& bCellID = mesh_.boundaryMesh()[patchID].faceCells()[facei];	
					const vector& d = X_nei[facei] - X_[bCellID];				
					Ugrad[bCellID] += Ainv_[bCellID] & (U_nei[facei]-U[bCellID]) * d;		
				}		
			}		
		}

		Ugrad.correctBoundaryConditions();
	}
	
}  


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
   	const GeometricField<vector, fvPatchField, volMesh>& U,
   	GeometricField<tensor, fvPatchField, volMesh>& Ugrad 
) 	const
{

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "dummy",
                mesh_
            ),
            mesh_,
            dimensioned<vector>("dummy", U.dimensions(), pTraits<vector>::zero)
        )
    );
 
	GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
	GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
	GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();


	gradientSchemes::gradient(U.component(0), UgradX);	
	gradientSchemes::gradient(U.component(1), UgradY);
	gradientSchemes::gradient(U.component(2), UgradZ);


	forAll (mesh_.cells(), cellID)
	{
		Ugrad[cellID] = tensor( UgradX[cellID], UgradY[cellID], UgradZ[cellID] );		
	}

	if( Pstream::parRun() ) 
	{	
		Ugrad.correctBoundaryConditions();
	}

	tvf.clear();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
   	const GeometricField<tensor, fvPatchField, volMesh>& U,
   	GeometricField<tensor, fvPatchField, volMesh>& UgradX,
	GeometricField<tensor, fvPatchField, volMesh>& UgradY,
	GeometricField<tensor, fvPatchField, volMesh>& UgradZ
) 	const
{

    const fvMesh& mesh = mesh_;

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "dummy("+U.name()+')',
                U.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
 
	GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
	GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
	GeometricField<vector, fvPatchField, volMesh> Uz = tvf();


	forAll (mesh.cells(), cellID)
	{
		Ux[cellID] = vector(U[cellID].xx(), U[cellID].xy(), U[cellID].xz());
		Uy[cellID] = vector(U[cellID].yx(), U[cellID].yy(), U[cellID].yz());
		Uz[cellID] = vector(U[cellID].zx(), U[cellID].zy(), U[cellID].zz());
	}

	if( Pstream::parRun() ) 
	{			
		Ux.correctBoundaryConditions();
		Uy.correctBoundaryConditions();
		Uz.correctBoundaryConditions();	
	}

	gradientSchemes::gradient(Ux, UgradX);	
	gradientSchemes::gradient(Uy, UgradY);
	gradientSchemes::gradient(Uz, UgradZ);

	tvf.clear();		
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::localGradient    
(
   	const GeometricField<vector, fvPatchField, volMesh>& U,
   	const GeometricField<vector, fvsPatchField, surfaceMesh>& Unei,   	
   	const GeometricField<tensor, fvPatchField, volMesh>& Ainv	 	
) const
{	
	
	const volVectorField& C = mesh_.C();
	const surfaceVectorField& Cf = mesh_.Cf();
	const scalar& maxInteriorFaceID = mesh_.nInternalFaces()-1;

    const objectRegistry& db = mesh_.thisDb();
	const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");


	tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject("dummy",mesh_),
            mesh_,
            dimensioned<vector>("dummy", U.dimensions(), pTraits<vector>::zero)
        )
    );

	GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
	GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
	GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();


	tmp<GeometricField<tensor, fvPatchField, volMesh> > tvf1
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject("dummy", mesh_),
            mesh_,
            dimensioned<tensor>("dummy", U.dimensions()/dimLength, pTraits<tensor>::zero)
        )
    );

	GeometricField<tensor, fvPatchField, volMesh> Ugrad = tvf1();


    forAll (mesh_.faces(), faceID)
    {
    	if (faceID <= maxInteriorFaceID)
    	{
			const label& ownID = mesh_.owner()[faceID];
			const label& neiID = mesh_.neighbour()[faceID];

			const vector& dOwn = Cf[faceID] - C[ownID];	
			const vector& dNei = Cf[faceID] - C[neiID];				

			UgradX[ownID] += Ainv[ownID] & ( (Unei[faceID].x()-U[ownID].x()) * dOwn );
			UgradY[ownID] += Ainv[ownID] & ( (Unei[faceID].y()-U[ownID].y()) * dOwn );
			UgradZ[ownID] += Ainv[ownID] & ( (Unei[faceID].z()-U[ownID].z()) * dOwn );	

			UgradX[neiID] += Ainv[neiID] & ( (Unei[faceID].x()-U[neiID].x()) * dNei );
			UgradY[neiID] += Ainv[neiID] & ( (Unei[faceID].y()-U[neiID].y()) * dNei );
			UgradZ[neiID] += Ainv[neiID] & ( (Unei[faceID].z()-U[neiID].z()) * dNei );	 
    	}

    	else
    	{
			const label& patchID = mesh_.boundaryMesh().whichPatch(faceID); 			
			const label& facei = mesh_.boundaryMesh()[patchID].whichFace(faceID);				
			const label& bCellID = mesh_.boundaryMesh()[patchID].faceCells()[facei];
							
			vector d = Cf.boundaryField()[patchID][facei] - C[bCellID];			

			UgradX[bCellID] += Ainv[bCellID] & ( (Unei.boundaryField()[patchID][facei].x()-U[bCellID].x()) * d );
			UgradY[bCellID] += Ainv[bCellID] & ( (Unei.boundaryField()[patchID][facei].y()-U[bCellID].y()) * d );
			UgradZ[bCellID] += Ainv[bCellID] & ( (Unei.boundaryField()[patchID][facei].z()-U[bCellID].z()) * d );	
    	
			if (lmN_.boundaryField().types()[patchID] == "fixedValue")
			{	
	    		const pointVectorField& xN_0_ = db.lookupObject<pointVectorField> ("xN_0"); 
				const label& faceID = mesh_.boundary()[patchID].start() + facei;

				forAll (mesh_.faces()[faceID], nodei)
				{
					const label& nodeID = mesh_.faces()[faceID][nodei];				
					
					vector d = xN_0_[nodeID] - C[bCellID];							
					UgradX[bCellID] += Ainv[bCellID] & ( (lmN_[nodeID].x()-U[bCellID].x()) * d );
					UgradY[bCellID] += Ainv[bCellID] & ( (lmN_[nodeID].y()-U[bCellID].y()) * d );
					UgradZ[bCellID] += Ainv[bCellID] & ( (lmN_[nodeID].z()-U[bCellID].z()) * d );						
						
					for (int i=0; i<7; i++)
					{										
						d = ( ( ((i+1)*xN_0_[nodeID]) + ((7-i)*Cf.boundaryField()[patchID][facei]) ) / 8.0 ) - C[bCellID];							
								
						UgradX[bCellID] += Ainv[bCellID] & ( (lmN_[nodeID].x()-U[bCellID].x()) * d );
						UgradY[bCellID] += Ainv[bCellID] & ( (lmN_[nodeID].y()-U[bCellID].y()) * d );
						UgradZ[bCellID] += Ainv[bCellID] & ( (lmN_[nodeID].z()-U[bCellID].z()) * d );				
					}									
				} 
	    	}
    	}
    }


    forAll (mesh_.cells(), cellID)
    {
    	Ugrad[cellID] = tensor (UgradX[cellID], UgradY[cellID], UgradZ[cellID]);	
    }


	tvf.clear();
	tvf1.clear();

	return Ugrad;
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
	GeometricField<scalar, fvPatchField, volMesh>& U,
	const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
	GeometricField<scalar, fvsPatchField, surfaceMesh>& Uminus,
	GeometricField<scalar, fvsPatchField, surfaceMesh>& Uplus
)
{
    const fvMesh& mesh = mesh_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
	const volVectorField& C = mesh.C();
	const surfaceVectorField& Cf = mesh.Cf();

	forAll (own, faceID)
	{		
		const label& ownID = own[faceID];
		const label& neiID = nei[faceID];

		Uminus[faceID] = U[ownID] + ( Ugrad[ownID] & (Cf[faceID]-C[ownID]) );
		Uplus[faceID]  = U[neiID] + ( Ugrad[neiID] & (Cf[faceID]-C[neiID]) );
	}


    forAll(mesh.boundary(), patchID)
    {
		forAll(mesh.boundaryMesh()[patchID],face)
		{
			const label& bCellID = mesh.boundaryMesh()[patchID].faceCells()[face];

			U.boundaryFieldRef()[patchID][face] = U[bCellID] + ( Ugrad[bCellID] & ( Cf.boundaryField()[patchID][face]-C[bCellID] ) );
			Uminus.boundaryFieldRef()[patchID][face] = U[bCellID] + ( Ugrad[bCellID] & ( Cf.boundaryField()[patchID][face]-C[bCellID] ) );		
		}    
	}
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
	GeometricField<vector, fvPatchField, volMesh>& U,
	const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
	GeometricField<vector, fvsPatchField, surfaceMesh>& Uminus,
	GeometricField<vector, fvsPatchField, surfaceMesh>& Uplus
)
{
    const fvMesh& mesh = mesh_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
	const volVectorField& C = mesh.C();
	const surfaceVectorField& Cf = mesh.Cf();


	forAll (own, faceID)
	{		
		const label& ownID = own[faceID];
		const label& neiID = nei[faceID];
	
		Uminus[faceID] = U[ownID] + ( Ugrad[ownID] & (Cf[faceID]-C[ownID]) );
		Uplus[faceID]  = U[neiID] + ( Ugrad[neiID] & (Cf[faceID]-C[neiID]) );
	}


    forAll(mesh.boundary(), patchID)
    {
		forAll(mesh.boundaryMesh()[patchID],face)
		{
			const label& bCellID = mesh.boundaryMesh()[patchID].faceCells()[face];
			
			Uminus.boundaryFieldRef()[patchID][face] = U[bCellID] + ( Ugrad[bCellID] & ( Cf.boundaryField()[patchID][face]-C[bCellID] ) );	
		}    
	}

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
	GeometricField<tensor, fvPatchField, volMesh>& U,
	const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
	const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
	const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
	GeometricField<tensor, fvsPatchField, surfaceMesh>& Uminus,
	GeometricField<tensor, fvsPatchField, surfaceMesh>& Uplus
)
{
    const fvMesh& mesh = mesh_;
    const labelUList& own = mesh.owner();
	const volVectorField& C = mesh.C();
	const surfaceVectorField& Cf = mesh.Cf();


    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "dummy("+U.name()+')',
                U.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
 
	GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
	GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
	GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

	forAll (mesh.cells(), cellID)
	{
		Ux[cellID] = vector( U[cellID].xx(), U[cellID].xy(), U[cellID].xz() );
		Uy[cellID] = vector( U[cellID].yx(), U[cellID].yy(), U[cellID].yz() );
		Uz[cellID] = vector( U[cellID].zx(), U[cellID].zy(), U[cellID].zz() );
	}
	

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "dummy("+Uminus.name()+')',
                Uminus.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<vector>("0", Uminus.dimensions(), pTraits<vector>::zero)
        )
    );

	GeometricField<vector, fvsPatchField, surfaceMesh> UminusX = tsf();
	GeometricField<vector, fvsPatchField, surfaceMesh> UminusY = tsf();
	GeometricField<vector, fvsPatchField, surfaceMesh> UminusZ = tsf();

	GeometricField<vector, fvsPatchField, surfaceMesh> UplusX = tsf();
	GeometricField<vector, fvsPatchField, surfaceMesh> UplusY = tsf();
	GeometricField<vector, fvsPatchField, surfaceMesh> UplusZ = tsf();


	forAll (own, faceID)
	{
		UminusX[faceID] = vector( Uminus[faceID].xx(), Uminus[faceID].xy(), Uminus[faceID].xz() );
		UminusY[faceID] = vector( Uminus[faceID].yx(), Uminus[faceID].yy(), Uminus[faceID].yz() );
		UminusZ[faceID] = vector( Uminus[faceID].zx(), Uminus[faceID].zy(), Uminus[faceID].zz() );

		UplusX[faceID] = vector( Uplus[faceID].xx(), Uplus[faceID].xy(), Uplus[faceID].xz() );
		UplusY[faceID] = vector( Uplus[faceID].yx(), Uplus[faceID].yy(), Uplus[faceID].yz() );
		UplusZ[faceID] = vector( Uplus[faceID].zx(), Uplus[faceID].zy(), Uplus[faceID].zz() );
	}

	gradientSchemes::reconstruct(Ux, UxGrad, UminusX, UplusX);
	gradientSchemes::reconstruct(Uy, UyGrad, UminusY, UplusY);
	gradientSchemes::reconstruct(Uz, UzGrad, UminusZ, UplusZ);

	forAll (own, faceID)
	{			
		Uminus[faceID] = tensor(UminusX[faceID], UminusY[faceID], UminusZ[faceID]);
		Uplus[faceID]  = tensor(UplusX[faceID], UplusY[faceID], UplusZ[faceID]);
	}

    forAll(mesh.boundary(), patchID)
    {
		forAll(mesh.boundaryMesh()[patchID],face)
		{
			const label& bCellID = mesh.boundaryMesh()[patchID].faceCells()[face];

			const vector& reconsX = Ux[bCellID] + ( UxGrad[bCellID] & ( Cf.boundaryField()[patchID][face]-C[bCellID] ) );	
			const vector& reconsY = Uy[bCellID] + ( UyGrad[bCellID] & ( Cf.boundaryField()[patchID][face]-C[bCellID] ) );
			const vector& reconsZ = Uz[bCellID] + ( UzGrad[bCellID] & ( Cf.boundaryField()[patchID][face]-C[bCellID] ) );		
		
			U.boundaryFieldRef()[patchID][face] =  tensor(reconsX, reconsY, reconsZ);
			Uminus.boundaryFieldRef()[patchID][face] =  tensor(reconsX, reconsY, reconsZ);
		}    
	}

	tsf.clear();
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
