//----------------------------------------------------------------------
// Preliminaries
//----------------------------------------------------------------------

    // Point Mesh (pMesh)
        pointMesh pMesh(mesh);


    // Material face area
        const surfaceScalarField& magSf (mesh.magSf());


    // Material face area normal vector
        const surfaceVectorField& Sf (mesh.Sf());

        
//----------------------------------------------------------------------
// Mesh coordinates
//----------------------------------------------------------------------

    // Material cell center coordinates
        const volVectorField& C = mesh.C();


    // Spatial cell center coordinates
        volVectorField x
        (
            IOobject ("x", mesh),
            C
        );


    // Spatial nodal coordinates 
        pointVectorField xN
        (
            IOobject ("xN", mesh),
            pMesh,
            dimensionedVector ("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
        );
        xN.primitiveFieldRef() = mesh.points();



    // Material nodal coordinates
        pointVectorField xN_0
        (
            IOobject ("xN_0", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            xN
        );      


    // Material face center coordinates
        const surfaceVectorField& Cf = mesh.Cf();               


    // Spatial face center coordinates 
        surfaceVectorField xF = mesh.Cf();



//----------------------------------------------------------------------
// Normals
//----------------------------------------------------------------------

    // Material normals     
        surfaceVectorField N 
        (
            IOobject ("N", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            Sf/mesh.magSf()
        );          




// EDGE TO POINT CONNECTIVITY (edges)
//-----------------------------------
    const edgeList& edges = mesh.edges();
    


// MINIMUM EDGE LENGTH FOR CFL NUMBER (h)
//---------------------------------------
    dimensionedScalar h ("h", dimensionSet(0,1,0,0,0,0,0), GREAT);

    forAll (edges,edgeID)
    {
        forAll (edges[edgeID],point)
        {
            const label& pointID = edges[edgeID][point];
            const label& pointID1 = edges[edgeID][point+1];
            const scalar& edgeLength = mag(mesh.points()[pointID]-mesh.points()[pointID1]);
            
            if (edgeLength < h.value())
            {
                h.value() = edgeLength;
            }                   
            
            break;
        }   
    }

    if( Pstream::parRun() ) 
    {
        reduce(h.value(), minOp<scalar>());
    }




// BOUNDARY PATCHES
//-----------------

    const polyBoundaryMesh& bm = mesh.boundaryMesh();
    const label& symmetricPatchID = bm.findPatchID("symmetric");
    const label& symmetricXpatchID = bm.findPatchID("symmetricX");
    const label& symmetricYpatchID = bm.findPatchID("symmetricY");
    const label& symmetricZpatchID = bm.findPatchID("symmetricZ");  



// A INVERSE FOR CALCULATION OF GRADIENT (Ainv)
//---------------------------------------------

    volTensorField Ainv                              
    (
        IOobject("Ainv", mesh),
        mesh,
        dimensionedTensor("Ainv",dimensionSet(0,2,0,0,0,0,0),tensor::zero)    
    );


    

// Ainv FOR CALCULATION OF GRADIENT FROM CONTACT VALUES (AvCinv)
//--------------------------------------------------------------

    // Initialise AvCinv    
    volTensorField AvCinv(Ainv);


