// Gmsh project created on Thu Jun  2 17:00:08 2016

/*
Inherit naming from the ratchet geometry

w = DefineNumber[ 4, Name "Parameters/w" ]; // radial width of the geometry 
L = DefineNumber[ 3, Name "Parameters/L" ]; // separation between plane surfaces 
h = DefineNumber[ 1, Name "Parameters/h" ]; // separation between upper surface and particle
a = DefineNumber[ 0.5, Name "Parameters/a" ]; // particle radius

st = DefineNumber[ 0.7, Name "Parameters/st" ]; // meshed space above upper plane 
sb = DefineNumber[ 0.5, Name "Parameters/sb" ]; // meshed space below lower plane

dsr = DefineNumber[ 0.5, Name "Parameters/dsr" ]; // ratchet charge sheet mesh size
dss = DefineNumber[ 0.05, Name "Parameters/dss" ]; // ratchet charge sheet mesh size
*/

// corners
Point(1) = {0, 0, 0, dsr};
Point(2) = {w, 0, 0, dsr};

Point(3) = {0, -sb, 0, dsr};
Point(4) = {w, -sb, 0, dsr};

Point(5) = {0, L, 0, dsr};
Point(6) = {w, L, 0, dsr};

Point(7) = {0, L+st, 0, dsr};
Point(8) = {w, L+st, 0, dsr};

//sphere
Point(9) = {0, L-h, 0, dss};
Point(10) = {0, L-h-a, 0, dss};
Point(11) = {0, L-h-a-a, 0, dss};
Point(12) = {a, L-h-a, 0, dss};

Line(1) = {3, 1};
Line(2) = {1, 11};
Line(3) = {11, 10};
Line(4) = {10, 9};
Line(5) = {9, 5};
Line(6) = {5, 7};
Line(7) = {7, 8};
Line(8) = {8, 6};
Line(9) = {6, 5};
Line(10) = {6, 2};
Line(11) = {2, 4};
Line(12) = {4, 3};
Line(13) = {1, 2};
Circle(14) = {11, 10, 12};
Circle(15) = {12, 10, 9};
Line Loop(16) = {9, 6, 7, 8};
Plane Surface(17) = {16};
Line Loop(18) = {9, -5, -15, -14, -2, 13, -10};
Plane Surface(19) = {18};
Line Loop(20) = {14, 15, -4, -3};
Plane Surface(21) = {20};
Line Loop(22) = {13, 11, 12, 1};
Plane Surface(23) = {22};

// domain labelling
Physical Surface(0) = {19}; // fluid
Physical Surface(1) = {17}; //Si
Physical Surface(2) = {23}; // poly
Physical Surface(3) = {21}; // sph

Physical Line(1) = {13}; // poly
Physical Line(2) = {9}; //si
Physical Line(3) = {15, 14}; //sph
Physical Line(4) = {4, 3, 5, 6, 7, 8, 10, 11, 12, 1, 2}; //edge


