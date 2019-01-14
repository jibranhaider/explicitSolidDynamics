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

#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
defineTypeNameAndDebug(constitutiveModel, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  
// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //  

constitutiveModel::constitutiveModel
(
    const volTensorField& F,
    const dictionary& dict  
)
:
    P_
    (
        IOobject
        (
            "P_",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE            
        ),
        F.mesh(),
        dimensionedTensor("P", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
    ),

    p_
    (
        IOobject
        (
            "p_",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE            
        ),
        F.mesh(),
        dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    Ealgo_
    (
        IOobject 
        (
            "Ealgo_",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        F.mesh(),
        dimensionedScalar("Ealgo", dimensionSet(1,-1,-2,0,0,0,0), 0.0)  
    ),

    rho_ (dict.lookup("rho")),
    E_ (dict.lookup("E")),
    nu_ (dict.lookup("nu")),

    mu_ (E_ / (2.0*(1.0 + nu_))),
    lambda_ ( nu_*E_ / ((1.0 + nu_)*(1.0 - 2.0*nu_)) ),
    kappa_ ( lambda_ + (2.0/3.0)*mu_ ),

    model_ ( dict.lookup("constitutiveModel") ),

    Up_ ( sqrt( (lambda_+2.0*mu_)/rho_ ) ),
    Us_ ( sqrt( mu_/rho_ ) )    
{}
    

// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
constitutiveModel::~constitutiveModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void constitutiveModel::correct()
{
    const fvMesh& mesh_ = P_.mesh();
    const objectRegistry& db = mesh_.thisDb();    
    const volTensorField& H_ = db.lookupObject<volTensorField> ("H");
    const volTensorField& F_ = db.lookupObject<volTensorField> ("F");           
    const volScalarField& J_ = db.lookupObject<volScalarField> ("J");       

    if (model_ == "linearElastic")
    {
        p_ = kappa_ * (tr(F_)-3.0);
        P_ = mu_ * ( F_ + F_.T() - ((2.0/3.0)*tr(F_)*tensor::I) ) + p_*tensor::I;
    }

    else if (model_ == "neoHookean")
    {   
        p_ = kappa_*(J_-1);
        P_ = mu_*pow(J_,(-2.0/3.0))*F_ - ( (mu_/3.0)*pow(J_,(-5.0/3.0))*(F_ && F_)*H_ ) + p_*H_;
    }

    else
    {
        FatalErrorIn("constitutiveModel.C")   
        << "Valid type entries are 'linearElastic' or 'neoHookean' for constitutiveModel"
        << abort(FatalError);       
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volScalarField constitutiveModel::Up_nonLinear()  
{
    const fvMesh& mesh_ = P_.mesh();
    const objectRegistry& db = mesh_.thisDb();
    const volTensorField& F_ = db.lookupObject<volTensorField> ("F");       
    volTensorField C = F_.T() & F_;
     
    tmp<GeometricField<scalar, fvPatchField, volMesh> > tsf
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject("dummy", mesh_),
            mesh_,
            dimensioned<scalar>("dummy", dimless, pTraits<scalar>::one)
        )
    );

    GeometricField<scalar, fvPatchField, volMesh> stretch = tsf();

    forAll (mesh_.cells(), cellID)
    {
        eigenStructure(C[cellID]);
        stretch[cellID] = min( eigVal_.x(), eigVal_.y());
        stretch[cellID] = Foam::sqrt( min( stretch[cellID], eigVal_.z()) );        
    }

    tsf.clear();

    return (Up_/stretch);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void constitutiveModel::eigenStructure(const tensor& ten)
{
    tensor t(ten);

    scalar it_max = 100;
    scalar it_num = 0;
    scalar rot_num = 0;
    scalar size = 3;
    scalar gapq = 0.0;
    
    int i,j;
    int k,l,m,p1,q = 0;

    double term, termp, termq = 0.0;
    double g,h,c,w1 = 0.0;
    double theta,thresh = 0.0;
    double s,t1,tau1 = 0.0;

    vector d( vector(t.xx(),t.yy(),t.zz()) );
    vector bw( vector(t.xx(),t.yy(),t.zz()) );  
    vector zw(vector::zero);
    tensor v1(tensor::I);


    while ( it_num < it_max )
    {
        it_num += 1;

        // The convergence threshold is based on the size of the elements in the strict upper triangle of the matrix.
        thresh = 0.0;

        for ( j = 0; j < size; j++ )
        {
          for ( i = 0; i < j; i++ )
          {
            thresh = thresh + t[i+j*size] * t[i+j*size];
          }
        }

        thresh = Foam::sqrt ( thresh ) / ( 4 * size );

        if ( thresh == 0.0 )
        {
          break;
        }

        for ( p1 = 0; p1 < size; p1++ )
        {
            for ( q = p1 + 1; q < size; q++ )
            {
                gapq = 10.0 * fabs( t[p1+q*size] );
                termp = gapq + fabs( d[p1] );
                termq = gapq + fabs( d[q] );
        
                // Annihilate tiny offdiagonal elements.
                if ( 4 < it_num && termp == fabs( d[p1] ) && termq == fabs( d[q] ) )
                {
                  t[p1+q*size] = 0.0;
                }
            
                //  Otherwise, apply a rotation.
                else if ( thresh <= fabs( t[p1+q*size] ) )
                {
                    h = d[q] - d[p1];
                    term = fabs( h ) + gapq;

                    if ( term == fabs( h ) )
                    {
                        t1 = t[p1+q*size] / h;
                    }
                    else
                    {
                        theta = 0.5 * h / t[p1+q*size];
                        t1 = 1.0 / ( fabs(theta) + Foam::sqrt( 1.0 + theta * theta ) );
                        if ( theta < 0.0 )
                        {
                            t1 = - t1;
                        }
                    }

                    c = 1.0 / Foam::sqrt ( 1.0 + t1 * t1 );
                    s = t1 * c;
                    tau1 = s / ( 1.0 + c );
                    h = t1 * t[p1+q*size];

                    //  Accumulate corrections to diagonal elements.
                    zw[p1] = zw[p1] - h;                 
                    zw[q] = zw[q] + h;
                    d[p1] = d[p1] - h;
                    d[q] = d[q] + h;
                    t[p1+q*size] = 0.0;

                    // Rotate, using information from the upper triangle of A only.
                    for ( j = 0; j < p1; j++ )
                    {
                        g = t[j+p1*size];
                        h = t[j+q*size];
                        t[j+p1*size] = g - s * ( h + g * tau1 );
                        t[j+q*size] = h + s * ( g - h * tau1 );
                    }

                    for ( j = p1 + 1; j < q; j++ )
                    {
                        g = t[p1+j*size];
                        h = t[j+q*size];
                        t[p1+j*size] = g - s * ( h + g * tau1 );
                        t[j+q*size] = h + s * ( g - h * tau1 );
                    }

                    for ( j = q + 1; j < size; j++ )
                    {
                        g = t[p1+j*size];
                        h = t[q+j*size];
                        t[p1+j*size] = g - s * ( h + g * tau1 );
                        t[q+j*size] = h + s * ( g - h * tau1 );
                    }

                    //  Accumulate information in the eigenvector matrix. 
                    for ( j = 0; j < size; j++ )
                    {
                        g = v1[j+p1*size];
                        h = v1[j+q*size];
                        v1[j+p1*size] = g - s * ( h + g * tau1 );
                        v1[j+q*size] = h + s * ( g - h * tau1 );
                    }
                    rot_num = rot_num + 1;
                }
            }
        }
    }


    // Restore upper triangle of input matrix
    for ( i = 0; i < size; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }

    //  Ascending sort the eigenvalues and eigenvectors
    for ( k = 0; k < size - 1; k++ )
    {
        m = k;
        for ( l = k + 1; l < size; l++ )
        {
            if ( d[l] < d[m] )
            {
                m = l;
            }
        }

        if ( m != k )
        {
            t1    = d[m];
            d[m] = d[k];
            d[k] = t1;
            for ( i = 0; i < size; i++ )
            {
                w1        = v1[i+m*size];
                v1[i+m*size] = v1[i+k*size];
                v1[i+k*size] = w1;
            }
        }
    }

    // Corrections for calculating inverse
    tensor sub(tensor::zero);

    for ( i=0; i<3; i++ )
    {
        if ( d[i] < SMALL )
        {
            d[i] = 1;
            sub += d[i] * vector(v1[3*i],v1[3*i+1],v1[3*i+2]) * vector(v1[3*i],v1[3*i+1],v1[3*i+2]);
        } 
    }

    eigVal_ = d;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceVectorField constitutiveModel::spatialNormal()  
{
    const fvMesh& mesh_ = P_.mesh();
    const objectRegistry& db = mesh_.thisDb();
    const volTensorField& F_ = db.lookupObject<volTensorField> ("F"); 

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject("tsf", mesh_),
            mesh_,
            dimensioned<vector>("tsf", dimless, pTraits<vector>::one)
        )
    );

    GeometricField<vector, fvsPatchField, surfaceMesh> n = tsf();

    surfaceTensorField FcInv = inv(fvc::interpolate(F_));    
    surfaceVectorField N = mesh_.Sf() / mesh_.magSf();

    n = (FcInv.T() & N) / ( mag(FcInv.T() & N) );   

    tsf.clear();

    return n;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //