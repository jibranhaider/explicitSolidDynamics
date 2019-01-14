//-----------------------------------------------------------------------------------
// Deformation dependent properties
//-----------------------------------------------------------------------------------           

    // Spatial normals  
    n = model.spatialNormal();
    nCn = n*n;
    iMnCn = (tensor::I)-(nCn);

    // Strain measures
    Finv = inv(F);       
    H = det(F)*Finv.T();    
    J = det(F); 

    // Wave speeds
    Up_time = model.Up_nonLinear();     


    // Constitutive model
    model.correct();
    p = model.pressure();
    P = model.piola();
    

//-----------------------------------------------------------------------------------
// Linear reconstruction
//-----------------------------------------------------------------------------------           

    // Cell gradients
    grad.gradient(lm, lmGrad);      
    grad.gradient(P, PxGrad, PyGrad, PzGrad);

    // Reconstruction   
    grad.reconstruct( lm, lmGrad, lm_M, lm_P ); 
    grad.reconstruct( P, PxGrad, PyGrad, PzGrad, P_M, P_P) ;
    t_M = P_M & N;
    t_P = P_P & N;


//-----------------------------------------------------------------------------------
// Riemann solver
//-----------------------------------------------------------------------------------

    #include "riemannSolver.H"              


//-----------------------------------------------------------------------------------
// Constrained approach
//-----------------------------------------------------------------------------------

    // Nodal linear momentum                    
    lmR = interpolate.surfaceToCell(lmC);
    lmRgrad = grad.localGradient(lmR, lmC, AvCinv);
    interpolate.volToPoint(lmR, lmRgrad, lmN);  
    #include "strongBCs.H"
    lmN.correctBoundaryConditions();    
    
    // Constrained fluxes
    interpolate.pointToSurface(lmN, lmC);