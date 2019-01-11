/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "interpolationSchemes.H"

#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(interpolationSchemes, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


  
// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //  
interpolationSchemes::interpolationSchemes(const fvMesh& vm)
:
    mesh_(vm)
{}
    
// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
interpolationSchemes::~interpolationSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volVectorField interpolationSchemes::surfaceToCell     
(
    const GeometricField<vector, fvsPatchField, surfaceMesh>& sf        
) const
{   
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    const scalar& maxInteriorFaceID = mesh_.nInternalFaces()-1;

    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");   


    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf_v
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceTocell("+sf.name()+')',
                sf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", sf.dimensions(), pTraits<vector>::zero)
        )
    );

    GeometricField<vector, fvPatchField, volMesh> U = tvf_v();

    tmp<GeometricField<scalar, fvPatchField, volMesh> > tvf_s
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceTocell("+sf.name()+')',
                sf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<scalar>("0", dimless, pTraits<scalar>::zero)
        )
    );

    GeometricField<scalar, fvPatchField, volMesh> w = tvf_s();

    forAll (mesh_.faces(), faceID)
    {
        if (faceID <= maxInteriorFaceID)
        {
            const label& ownID = mesh_.owner()[faceID];
            const label& neiID = mesh_.neighbour()[faceID];

            const vector& dOwn = Cf[faceID] - C[ownID]; 
            const vector& dNei = Cf[faceID] - C[neiID];             

            U[ownID] += sf[faceID] * (1.0/mag(dOwn));
            U[neiID] += sf[faceID] * (1.0/mag(dOwn));

            w[ownID] += 1.0/mag(dOwn);
            w[neiID] += 1.0/mag(dNei);
        }
        else
        {
            const label& patchID = mesh_.boundaryMesh().whichPatch(faceID);         
            const label& facei = mesh_.boundaryMesh()[patchID].whichFace(faceID);               
            const label& bCellID = mesh_.boundaryMesh()[patchID].faceCells()[facei];
            vector d = Cf.boundaryField()[patchID][facei] - C[bCellID];

            U[bCellID] += sf.boundaryField()[patchID][facei] * (1.0/mag(d));
            w[bCellID] += 1.0/mag(d);
 
            if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            {   
                const pointVectorField& xN_0_ = db.lookupObject<pointVectorField> ("xN_0");     
                const label& faceID = mesh_.boundary()[patchID].start() + facei;

                forAll (mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei];             
            
                    d = xN_0_[nodeID] - C[bCellID];                         
                    U[bCellID] += lmN_[nodeID] * (1.0/mag(d));
                    w[bCellID] += (1.0/mag(d));                     
                    
                    for (int i=0; i<7; i++)
                    {                                       
                        d = ( ( ((i+1)*xN_0_[nodeID]) + ((7-i)*Cf.boundaryField()[patchID][facei]) ) / 8.0 ) - C[bCellID];                          
                
                        U[bCellID] +=   lmN_[nodeID] * (1.0/mag(d));
                        w[bCellID] += (1.0/mag(d));                 
                    }
                }
            }
        }
    }

    U.primitiveFieldRef() /= w.internalField();


    tvf_v.clear();
    tvf_s.clear();

    return U;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void interpolationSchemes::pushUntransformedDataNew
(
    List<Type>& pointData
) const
{
    const fvMesh& mesh = mesh_; 
    
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Combine master data with slave data
    forAll(slaves, i)
    {
        const labelList& slavePoints = slaves[i];

        // Copy master data to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elems[i];
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void interpolationSchemes::addSeparatedNew
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{

     if (debug)
     {
         Pout<< "volPointInterpolation::addSeparated" << endl;
     }
 
     typename GeometricField<Type, pointPatchField, pointMesh>::
         Internal& pfi = pf.ref();
 
     typename GeometricField<Type, pointPatchField, pointMesh>::
         Boundary& pfbf = pf.boundaryFieldRef();
 
     forAll(pfbf, patchi)
     {
         if (pfbf[patchi].coupled())
         {
             refCast<coupledPointPatchField<Type>>
                 (pfbf[patchi]).initSwapAddSeparated
                 (
                     Pstream::commsTypes::nonBlocking,
                     pfi
                 );
         }
     }
 
     // Block for any outstanding requests
     Pstream::waitRequests();
 
     forAll(pfbf, patchi)
     {
         if (pfbf[patchi].coupled())
         {
             refCast<coupledPointPatchField<Type>>
                 (pfbf[patchi]).swapAddSeparated
                 (
                     Pstream::commsTypes::nonBlocking,
                     pfi
                 );
         }
    }

}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interpolationSchemes::volToPoint   
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, pointPatchField, pointMesh>& Un                  
) const
{
    
    const fvMesh& mesh = mesh_; 
    const volVectorField& C = mesh.C();    


    if( Pstream::parRun() ) 
    {           
        pointScalarField sum
        (
            IOobject
            (
                "volPointSum",
                mesh.polyMesh::instance(),
                mesh
            ),
            pointMesh::New(mesh),
            dimensionedScalar("zero", dimless, 0.0)
        );      

        forAll (mesh.points(), nodeID)
        {
            Un[nodeID] = vector::zero;
        }   

        forAll (mesh.points(), nodeID)
        {
            forAll (mesh.pointCells()[nodeID], cell)
            {   
                const label& cellID = mesh.pointCells()[nodeID][cell];
                const vector& d = mesh.points()[nodeID] - C[cellID];
                const vector& recons = U[cellID] + ( Ugrad[cellID] & d );   
                const scalar& weight = 1;
                
                Un[nodeID] += recons;
                sum[nodeID] += weight;      
            }
        }       

        pointConstraints::syncUntransformedData(mesh, sum, plusEqOp<scalar>());
        addSeparatedNew(sum);
        pushUntransformedDataNew(sum);
    
        forAll (mesh.points(), nodeID)
        {       
            Un[nodeID] = Un[nodeID] / sum[nodeID];
        }
    
        pointConstraints::syncUntransformedData(mesh, Un, plusEqOp<vector>());
        addSeparatedNew(Un);
        pushUntransformedDataNew(Un);
    }
    
    else
    {
        forAll (mesh.pointCells(), nodeID)
        {
            vector sum = vector::zero;
            scalar weights = 0.0;
            
            forAll (mesh.pointCells()[nodeID], cell)
            {
                const label& cellID = mesh.pointCells()[nodeID][cell];
                const vector& d = mesh.points()[nodeID] - C[cellID];
                const vector& recons = U[cellID] + ( Ugrad[cellID] & d );
                    
                //sum += recons * (1.0/mag(d));
                //weights += (1.0/mag(d));
                
                sum += recons;
                weights += 1.0;
            }
            
            Un[nodeID] = sum / weights; 
        }       
    }   

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interpolationSchemes::pointToSurface
(
    const GeometricField<vector, pointPatchField, pointMesh>& U,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Uf                  
) const
{
    const fvMesh& mesh = mesh_;
    const scalar& maxInteriorFaceID = mesh.nInternalFaces()-1;  

    forAll (mesh.faces(),faceID)
    {
        vector sum = vector::zero;
        scalar weights = 0.0;

        if (faceID <= maxInteriorFaceID)
        {   
            forAll (mesh.faces()[faceID],node)
            {
                const label& nodeID = mesh.faces()[faceID][node];
                sum += U[nodeID];
                weights += 1.0;                             
            }
            
            Uf[faceID] = sum / weights; 
        }
        
        else
        {
            const label& patchID = mesh.boundaryMesh().whichPatch(faceID);
            const label& facei = mesh.boundaryMesh()[patchID].whichFace(faceID);

            forAll (mesh.faces()[faceID],node)
            {
                const label& nodeID = mesh.faces()[faceID][node];       
                sum += U[nodeID];
                weights += 1.0;                             
            }
            
            Uf.boundaryFieldRef()[patchID][facei]= sum / weights;           
        }   
    }
    
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interpolationSchemes::surfaceToPointConstrained
(
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Uf,
    GeometricField<vector, pointPatchField, pointMesh>& Un              
) const
{



}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
