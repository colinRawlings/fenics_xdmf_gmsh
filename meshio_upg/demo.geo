//#gmsh file
Point(6) = {0, 0, 0, 1};
Point(7) = {2, 0, 0, 1};
Point(8) = {0, 2, 0, 1};
Point(9) = {2, 2, 0, 1};
Line(10) = {6, 7};               
Line(20) = {7, 9};           
Line(30) = {9, 8};       
Line(40) = {8, 6};       
Line Loop(1) = {10, 20, 30, 40};

Plane Surface(1) = {1};
Physical Surface(0) = {1};

Physical Line(1) = {40};
Physical Line(2) = {20};
Physical Line(4) = {30};
Physical Line(3) = {10};