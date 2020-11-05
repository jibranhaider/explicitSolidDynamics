// MESH

//c = 1;	// 6x10x3
c = 2;		// 12x20x6
//c = 3;	// 18x30x9 (doesn't work)
//c = 4;	// 24x40x12



// OUTPUTS

// Section 1
Point(1) = {6, 0, 0};
Point(2) = {6, 3, 0};
Point(3) = {3, 3, 0};
Point(4) = {3, 0, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Transfinite Line{1,2,3,4} = 3*c+1;

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface{1} = {1,2,3,4};
Recombine Surface{1};


// Section 2
Point(5) = {3, 10, 0};
Point(6) = {0, 10, 0};
Point(7) = {0, 0,  0};

Line(5) = {3,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,4};

Transfinite Line{5} = 7*c+1;
Transfinite Line{6,8} = 3*c+1;
Transfinite Line{7} = 10*c+1;


Line Loop(2) = {5,6,7,8,-3};
Plane Surface(2) = {2};
Transfinite Surface{2} = {5,6,7,4};
Recombine Surface{2};


// Extrusion (Section 1 & 2)
Extrude{0,0,3}{
Surface{1,2};
Layers{3*c};Recombine;
}


// Definition of surfaces for boundary conditions
Physical Surface("loading1") 	= {17};
Physical Surface("loading2") 	= {44};
Physical Surface("free") 		= {21,29,30,1,40,48,52,57,2};

// Definition of a volume
Physical Volume("volume") = {1,2};

