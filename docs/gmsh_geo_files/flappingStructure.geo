//INPUTS


// Cells
cd = 30;
cc = 8;
cz = 10;
p  = 1.03;
h  = 0.2; 

// Dimensions
ri = 0.2 ; 		// inner radius
wi = (25*Pi)/180;	// inner angle
ro = 1;				// outer radius
wo = (15*Pi)/180;	// outer angle




// OUTPUTS

// SECTION #01
section = 1;
Point(1) = {0, 0, 0};
Point(2) = {ri*Sin(wi), ri*Cos(wi), 0};
Point(3) = {ro*Sin(wo), ro*Cos(wo), 0};
Point(4) = {-ro*Sin(wo), ro*Cos(wo), 0};
Point(5) = {-ri*Sin(wi), ri*Cos(wi), 0};

Line(1) = {2,3};
Circle(2) = {3,1,4};
Line(3) = {4,5};
Circle(4) = {5,1,2};
Transfinite Line{1,-3} = cd+1 Using Progression p;
Transfinite Line{2,4} = cc+1;

Line Loop(section) = {1,2,3,4} ;
Plane Surface(section) = {section};
Transfinite Surface{section} = {2,3,4,5};
Recombine Surface{section};


// SECTION #02
section = 2;
Point(6) = {ri*Cos(wi), ri*Sin(wi), 0};
Point(7) = {ro*Cos(wo), ro*Sin(wo), 0};
Point(8) = {ro*Cos(wo), -ro*Sin(wo), 0};
Point(9) = {ri*Cos(wi), -ri*Sin(wi), 0};

Line(5) = {6,7};
Circle(6) = {7,1,8};
Line(7) = {8,9};
Circle(8) = {9,1,6};
Transfinite Line{5,-7} = cd+1 Using Progression p;
Transfinite Line{6,8} = cc+1;

Line Loop(section) = {5,6,7,8} ;
Plane Surface(section) = {section};
Transfinite Surface{section} = {6,7,8,9};
Recombine Surface{section};


// SECTION #03
section = 3;
Point(10) = {ri*Sin(wi), -ri*Cos(wi), 0};
Point(11) = {ro*Sin(wo), -ro*Cos(wo), 0};
Point(12) = {-ro*Sin(wo), -ro*Cos(wo), 0};
Point(13) = {-ri*Sin(wi), -ri*Cos(wi), 0};

Line(9) = {10,11};
Circle(10) = {11,1,12};
Line(11) = {12,13};
Circle(12) = {13,1,10};
Transfinite Line{9,-11} = cd+1 Using Progression p;
Transfinite Line{10,12} = cc+1;

Line Loop(section) = {9,10,11,12} ;
Plane Surface(section) = {section};
Transfinite Surface{section} = {10,11,12,13};
Recombine Surface{section};


// SECTION #04
section = 4;
Point(14) = {-ri*Cos(wi), -ri*Sin(wi), 0};
Point(15) = {-ro*Cos(wo), -ro*Sin(wo), 0};
Point(16) = {-ro*Cos(wo), ro*Sin(wo), 0};
Point(17) = {-ri*Cos(wi), ri*Sin(wi), 0};

Line(13) = {14,15};
Circle(14) = {15,1,16};
Line(15) = {16,17};
Circle(16) = {17,1,14};
Transfinite Line{13,-15} = cd+1 Using Progression p;
Transfinite Line{14,16} = cc+1;

Line Loop(section) = {13,14,15,16} ;
Plane Surface(section) = {section};
Transfinite Surface{section} = {14,15,16,17};
Recombine Surface{section};


// SECTION #05
section = 5;
Circle(17) = {2,1,6};
Circle(18) = {9,1,10};
Circle(19) = {13,1,14};
Circle(20) = {17,1,5};

Line Loop(section) = {4,17,-8,18,-12,19,-16,20} ;
Plane Surface(section) = {section};
Transfinite Line{17,18,19,20} = cc + 1;
Recombine Surface{section};


// Extrusion in the third dimension
Extrude{0,0,h}{
Surface{1,2,3,4,5};
Layers{cz};Recombine;
}

// Definition of surfaces for boundary conditions
Physical Surface("free") 	= {150,3,86,2,64,1,42,4,108,141,133,125,149,81,73,59,51,29,37,103,95,77,55,33,99};
Physical Surface("symmetricZ") 	= {5};


// Definition of a volume
Physical Volume("volume") = {1,2,3,4,5};



