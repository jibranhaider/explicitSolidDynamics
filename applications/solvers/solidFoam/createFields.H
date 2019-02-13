Info << "Creating mesh normal fields ..." << endl;

 // Material normals
 const surfaceVectorField N = Sf/mesh.magSf();

 // Spatial normals
 surfaceVectorField n
 (
    IOobject("n", mesh),
    N
 );

// Normal projection tensor
surfaceTensorField nCn
(
    IOobject("nCn", mesh),
    n*n
);

// Tangential projection tensor
surfaceTensorField iMnCn
(
    IOobject("iMnCn", mesh),
    tensor::I - nCn
);


Info << "Creating linear momentum fields ..." << endl;

// Cell linear momentum
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

// Nodal linear momentum
pointVectorField lmN
(
    IOobject
    (
        "lmN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh
);


Info << "Creating strain measure fields ..." << endl;

// Deformation gradient tensor
volTensorField F
(
    IOobject("F", mesh),
    mesh,
    Foam::tensor::I
);

// Inverse of deformation gradient
volTensorField Finv(inv(F));

// Cofactor of deformation
volTensorField H
(
    IOobject("H", mesh),
    det(F)*Finv.T()
);

// Jacobian of deformation
volScalarField J
(
    IOobject("J", mesh),
    det(F)
);


Info << "Creating constitutive model ..." << endl;

// Solid model class
constitutiveModel model(F, mechanicalProperties);
model.correct();

// Pressure
volScalarField p
(
    IOobject ("p", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    model.pressure()
);
p.write();

// First Piola Kirchhoff stress tensor
volTensorField P = model.piola();


Info << "Creating fields for wave speeds ..." << endl;

// Longitudinal wave speed
volScalarField Up
(
    IOobject("Up", mesh),
    mesh,
    model.Up()/beta
);

// Shear wave speed
volScalarField Us
(
    IOobject ("Us", mesh),
    mesh,
    model.Us()*beta
);

// Wave speed for time increment
volScalarField Up_time = model.Up_nonLinear();


Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh);
grad.distanceMatrixLocal(AvCinv);

// Gradient of cell linear momentum
volTensorField lmGrad
(
    IOobject("lmGrad", mesh),
    mesh,
    dimensionedTensor("lmGrad", lm.dimensions()/dimLength, tensor::zero)
);
grad.gradient(lm, lmGrad);

// Components of PK1
volVectorField Px
(
    IOobject("Px", mesh),
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


Info << "Creating fields for reconstruction ..." << endl;

// Reconstruction of linear momentum
surfaceVectorField lm_M
(
    IOobject
    (
        "lm_M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("lm_M", lm.dimensions(), vector::zero)
);
surfaceVectorField lm_P = lm_M;
grad.reconstruct(lm, lmGrad, lm_M, lm_P);

// Reconstruction of PK1 stresses
surfaceTensorField P_M
(
    IOobject("P_M", mesh),
    mesh,
    dimensionedTensor("P_M", P.dimensions(), tensor::zero)
);
surfaceTensorField P_P = P_M;
grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P);

// Reconstruction of traction
surfaceVectorField t_M
(
    IOobject
    (
        "t_M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    P_M & N
);
surfaceVectorField t_P = P_P & N;


Info << "Creating fields for riemann solver ..." << endl;

// Contact linear momentum
surfaceVectorField lmC =
(
    0.5*(lm_M + lm_P) + (0.5*((nCn/fvc::interpolate(Up))
  + (iMnCn/fvc::interpolate(Us))) & (t_P - t_M))
);

volVectorField lm_b
(
    IOobject
    (
        "lm_b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Contact traction
surfaceVectorField tC =
(
    0.5*(t_M + t_P) + (0.5 * ((fvc::interpolate(Up)*nCn)
  + (fvc::interpolate(Us)*iMnCn)) & (lm_P - lm_M))
);

volVectorField t_b
(
    IOobject
    (
        "t_b",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Riemann solver
#include "riemannSolver.H"


Info << "Creating fields for the constrained procedure ..." << endl;

// Constrained class
interpolationSchemes interpolate(mesh);

// Cell averaged linear momentum
volVectorField lmR = interpolate.surfaceToCell(lmC);

// Local gradient of cell averaged linear momentum
volTensorField lmRgrad = grad.localGradient(lmR, lmC, AvCinv);

// Constrained fluxes
interpolate.volToPoint(lmR, lmRgrad, lmN);
#include "strongBCs.H"
lmN.correctBoundaryConditions();
interpolate.pointToSurface(lmN, lmC);


Info << "Creating fields for angular momentum ..." << endl;

// RHS of linear momentum equation
volVectorField rhsLm
(
    IOobject("rhsLm", mesh),
    mesh,
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

// RHS of linear momentum equation (stage 1)
volVectorField rhsLm1 = rhsLm;

// RHS of angular momentum equation
volVectorField rhsAm
(
    IOobject("rhsAm", mesh),
    mesh,
    dimensionedVector("rhsAm", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);


Info << "Creating fields for post-processing ..." << endl;

// Nodal displacements
pointVectorField uN
(
    IOobject
    (
        "uN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("uN", dimLength, vector::zero)
);
uN.write();


Info << "Creating variables for time ..." << endl;

// Time
dimensionedScalar t("t", dimTime, 0.0);

// Time step
scalar tstep = 0;

// Time increment
dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());

if (timeStepping == "variable")
{
    deltaT = (cfl*h)/max(Up_time);
    runTime.setDeltaT(deltaT);
}

// Runge-Kutta stage
word RKstage("first");


//-----------------------------------------------------------------------------------
// Material parameters
//-----------------------------------------------------------------------------------

const dimensionedScalar& rho = model.density();
angularMomentum AM(mesh, rho);