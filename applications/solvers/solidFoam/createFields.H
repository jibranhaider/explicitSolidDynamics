//-----------------------------------------------------------------------------------
// Linear momentum
//-----------------------------------------------------------------------------------
    
    // Cel linear momentum
    volVectorField lm                                        
    (
        IOobject ("lm",runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)             
    );

    // Nodal linear momentum
    pointVectorField lmN
    (
        IOobject ("lmN", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        pMesh
    );


//-----------------------------------------------------------------------------------
// Deformation dependent variables
//-----------------------------------------------------------------------------------

    // Deformation gradient tensor
    volTensorField F                                         
    (
        IOobject ("F", mesh),
        mesh,
        Foam::tensor::I                 
    );


    // Spatial normals
    surfaceVectorField n
    (
        IOobject ("n", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        N  
    );


    // Normal projection tensor
    surfaceTensorField nCn
    (
        IOobject ("nCn", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        n*n   
    );


    // Tangential projection tensor
    surfaceTensorField iMnCn
    (
        IOobject ("iMnCn", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        tensor::I - nCn
    );


    // Inverse of deformation gradient
    volTensorField Finv(inv(F));        


    // Cofactor of deformation
    volTensorField H
    (           
        IOobject ("H", mesh),
        det(F) * Finv.T()
    );


    // Jacobian of deformation                               
    volScalarField J                                         
    (
        IOobject ("J", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        det(F)
    );


//-----------------------------------------------------------------------------------
// Constitutive model
//-----------------------------------------------------------------------------------

    // Constitutive model class
    constitutiveModel model (F, mechanicalProperties);  
    model.correct();


    // Pressure
    volScalarField p                                         
    (
        IOobject ("p", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        model.pressure()  
    );
    p.write();


    // First Piola Kirchhoff stress tensor 
    volTensorField P (model.piola());


//-----------------------------------------------------------------------------------
// Wave speeds
//-----------------------------------------------------------------------------------

    // Longitudinal wave speed
    volScalarField Up
    (           
        IOobject ("Up", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), 
        mesh,
        model.Up() / beta
    );


    // Shear wave speed
    volScalarField Us
    (           
        IOobject ("Us", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), 
        mesh,
        model.Us() * beta
    );  


    // Wave speed for time increment
    volScalarField Up_time (model.Up_nonLinear());



//-----------------------------------------------------------------------------------
// Gradients
//-----------------------------------------------------------------------------------

    // Gradient class
    gradientSchemes grad(mesh);
    grad.distanceMatrixLocal(AvCinv);


    // Gradient of cell linear momentum     
    volTensorField lmGrad                                        
    (
        IOobject ("lmGrad", mesh),
        mesh,
        dimensionedTensor("lmGrad", lm.dimensions()/dimLength, tensor::zero)    
    );  
    grad.gradient(lm, lmGrad);      


    // Components of PK1
    volVectorField Px
    (           
        IOobject ("Px", mesh),  
        mesh,
        dimensionedVector("Px", P.dimensions(), vector::zero)
    );  
    volVectorField Py(Px);
    volVectorField Pz(Px);      


    // Gradient of PK1
    volTensorField PxGrad (fvc::grad(Px));  
    volTensorField PyGrad (fvc::grad(Py));
    volTensorField PzGrad (fvc::grad(Pz));
    grad.gradient(P, PxGrad, PyGrad, PzGrad);



//-----------------------------------------------------------------------------------
// Reconstruction
//-----------------------------------------------------------------------------------


    // Reconstruction of linear momentum
    surfaceVectorField lm_M
    (
        IOobject ("lm_M", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("lm_M", lm.dimensions(), vector::zero)   
    );  
    surfaceVectorField lm_P (lm_M);
    grad.reconstruct( lm, lmGrad, lm_M, lm_P ); 


    // Reconstruction of PK1 stresses
    surfaceTensorField P_M
    (
        IOobject ("P_M", mesh),
        mesh,
        dimensionedTensor("P_M", P.dimensions(), tensor::zero) 
    );          
    surfaceTensorField P_P (P_M);
    grad.reconstruct( P, PxGrad, PyGrad, PzGrad, P_M, P_P) ;


    // Reconstruction of traction
    surfaceVectorField t_M
    (
        IOobject ("t_M",runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        P_M & N 
    );
    surfaceVectorField t_P (P_P & N);   


//-----------------------------------------------------------------------------------
// Riemann solver
//-----------------------------------------------------------------------------------


    // Contact linear momentum
    surfaceVectorField lmC
    (
        IOobject ("lmC", mesh),
        0.5*(lm_M+lm_P) + ( 0.5 * ((nCn/fvc::interpolate(Up)) + (iMnCn/fvc::interpolate(Us))) & (t_P-t_M) )       
    );

    volVectorField lm_b
    (
        IOobject ("lm_b", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );


    // Contact traction
    surfaceVectorField tC
    (
        IOobject ("tC", mesh), 
        0.5*(t_M+t_P) + ( 0.5 * ((fvc::interpolate(Up)*nCn) + (fvc::interpolate(Us)*iMnCn)) & (lm_P-lm_M) ) 
    );

    volVectorField t_b
    (
        IOobject ("t_b", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh 
    );  


    #include "riemannSolver.H"


//-----------------------------------------------------------------------------------
// Constrained approach
//-----------------------------------------------------------------------------------

    // Constrained class
    interpolationSchemes interpolate(mesh); 


    // Cell averaged linear momentum
    volVectorField lmR (interpolate.surfaceToCell(lmC));                                     


    // Local gradient of cell averaged linear momentum
    volTensorField lmRgrad (grad.localGradient(lmR, lmC, AvCinv));  


    // Constrained fluxes
    interpolate.volToPoint(lmR, lmRgrad, lmN);      
    #include "strongBCs.H"
    lmN.correctBoundaryConditions();
    interpolate.pointToSurface(lmN, lmC);       



//-----------------------------------------------------------------------------------
// Angular momentum
//-----------------------------------------------------------------------------------

    // RHS of linear momentum equation
    volVectorField rhsLm                                         
    (
        IOobject ("rhsLm", mesh),
        mesh,
        dimensionedVector("rhsLm",dimensionSet(1,-2,-2,0,0,0,0),vector::zero)    
    );  
    

    // RHS of linear momentum equation (stage 1)
    volVectorField rhsLm1 (rhsLm);


    // RHS of angular momentum equation
    volVectorField rhsAm                                 
    (
        IOobject ("rhsAm", mesh),
        mesh,
        dimensionedVector("rhsAm",dimensionSet(1,-1,-2,0,0,0,0),vector::zero)    
    );


//-----------------------------------------------------------------------------------
// Post processing
//-----------------------------------------------------------------------------------

    // Nodal displacements
    pointVectorField uN
    (
        IOobject ("uN", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pMesh,
        dimensionedVector("uN", dimLength, vector::zero)
    );
    uN.write();



//-----------------------------------------------------------------------------------
// Simulation parameters
//-----------------------------------------------------------------------------------

    scalar tstep = 0;
    dimensionedScalar deltaT ("deltaT", dimTime, runTime.deltaTValue()) ;
    dimensionedScalar t ("t", dimTime, 0.0);

    if (timeStepping == "variable")
    {
        deltaT = ( cfl*h ) / max(Up_time);
        runTime.setDeltaT(deltaT);
    }

    word RKstage("first");


//-----------------------------------------------------------------------------------
// Material parameters
//-----------------------------------------------------------------------------------

    const dimensionedScalar& rho = model.density();
    angularMomentum AM(mesh, rho);

