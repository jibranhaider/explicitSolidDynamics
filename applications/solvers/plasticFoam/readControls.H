Info << "Reading dictionaries ..." << endl;

// Mechanical properties
IOdictionary mechanicalProperties
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Run parameters dictionary
IOdictionary runParameters
(
    IOobject
    (
        "runParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Control dictionary
IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Finite volume solution dictionary
IOdictionary fvSolution
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info << "Reading data from dictionaries ..." << endl;

// Type of time stepping
const word& timeStepping = controlDict.lookup("timeStepping");

if (timeStepping != "constant" && timeStepping != "variable")
{
    FatalErrorIn("readControls.H")
        << "Valid type entries are 'constant' or 'variable' for timeStepping"
        << abort(FatalError);
}

// Courant-Friedrichs-Lewy Number
const scalar& cfl = readScalar(controlDict.lookup("cfl"));

if (cfl <= 0.0 || cfl > 1.0)
{
    FatalErrorIn("readControls.H")
        << "Valid type entries are '<= 1' or '> 0' for cfl"
        << abort(FatalError);
}

// Stabilisation parameter for near incompressibility
const scalar& beta =
    readScalar(fvSolution.lookup("incompressiblilityCoefficient"));

// Angular momentum conservation
const word& angularMomentumConservation =
    fvSolution.lookup("angularMomentumConservation");

if
(
    angularMomentumConservation != "yes" && angularMomentumConservation != "no"
)
{
    FatalErrorIn("readControls.H")
        << "Valid type entries are 'yes' or 'no' "
        << "for angularMomentumConservation"
        << abort(FatalError);
}