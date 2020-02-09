//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Rectangle(1) = {-0.5, -0.5, 0, 1, 1, 0};
//+
Physical Curve("outer_boundary", 1) = {1};
//+
Physical Curve("inner_boundary", 2) = {4, 5, 2, 3};
//+
Physical Surface("inner_surface", 2) = {1};

//+
Curve Loop(2) = {1};
//+
Curve Loop(3) = {4, 5, 2, 3};
//+
Plane Surface(2) = {2, 3};
//+
Physical Surface("outer_surface", 1) = {2};
//+
Characteristic Length {1} = 0.2;
//+
Characteristic Length {5, 2, 3, 4} = 0.05;
