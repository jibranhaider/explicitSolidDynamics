// Number of cells
c = 4;      // 384 cells
//c = 8;    // 3072 cells
//c = 16;   // 24576 cells
//c = 32;   // 196608 cells

// Lengths
lx = 1;
ly = 6;
lz = 1;

// Number of cells
cx = c;
cy = 6*c;
cz = c;


Point(1) = {-lx/2, 0, -lz/2};

Extrude {lx,0,0} {
     Point{1}; Layers{cx}; Recombine;
}

Extrude {0,ly,0} {
     Line{1}; Layers{cy}; Recombine;
}

Extrude {0,0,lz} {
     Surface{5}; Layers{cz}; Recombine;
}

// Definition of surfaces for boundary conditions
Physical Surface("topAndSides") = {5,27,26,18,22};
Physical Surface("bottom") = {14};

// Definition of a volume
Physical Volume("volume") = {1};
