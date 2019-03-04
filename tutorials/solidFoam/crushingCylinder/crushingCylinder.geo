// Number of cells
n = 100;        // 40,000 cells
//n = 150;      // 90,000 cells
//n = 200;      // 160,000 cells
//n = 250;      // 250,000 cells

// cc = Circle circumference
// cy = Cylinder axis

cc = n; cr = 2; cy = n;

// Domain limits
ro = 0.1;           // Radius of circle
ri = 0.099753;      // Square length
ymax = 0.1139;      // Bar height

// Nodes
nc = cc + 1; nr = cr + 1; ny = cy + 1;


// Section
Point(1) = {0, 0, 0};
Point(2) = {0, 0, ri};
Point(3) = {ri, 0, 0};
Point(4) = {ro, 0, 0};
Point(5) = {0, 0, ro};
Circle(1) = {2,1,3};
Line(2) = {3,4};
Circle(3) = {4,1,5};
Line(4) = {5,2};
Transfinite Line{1,3} = nc;
Transfinite Line{2,4} = nr;
Line Loop(1) = {1,2,3,4} ;
Plane Surface(1) = {1};
Transfinite Surface{1} = {2,3,4,5};
Recombine Surface{1};

// Extrusion
Extrude{0,ymax,0}{
Surface{1};
Layers{2*n};Recombine;
}

// Definition of surfaces for boundary conditions
Physical Surface("free") = {13,21};
Physical Surface("loading") = {26};
Physical Surface("symmetric") = {17, 25};
Physical Surface("fixed") = {1};

// Definition of a volume
Physical Volume("volume") = {1};
