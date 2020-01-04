SetFactory("OpenCASCADE");
//
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 0.5, 0, 1.0};
//+
Point(3) = {0.5, 0, 0, 1.0};
//+
Point(4) = {0, -0.5, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
//+
Point(7) = {1, 0, 0, 1.0};
//+
Line(1) = {5, 4};
//+
Line(2) = {4, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 6};
//+

//+
Circle(5) = {2, 1, 3};
//+
Circle(6) = {3, 1, 4};
//+
Circle(7) = {6, 1, 7};
//+
Circle(8) = {7, 1, 5};
//+
Curve Loop(1) = {5, 6, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 8, 1, -6, -5, 4};
//+
Plane Surface(2) = {2};
//+
//+
Extrude {{0, 1, 0}, {0, 0, 0}, Pi} {
  Surface{2}; Surface{1}; 
}
//+
Physical Surface("inner", 2) = {6, 5};
//+
Physical Surface("outer", 1) = {3, 4};
//+
Physical Volume("outer", 1) = {1};
//+
Physical Volume("inner", 2) = {2};
//+
Physical Point("center", 3) = {1};
//+
Characteristic Length {2, 9, 1, 4, 3} = 0.05;
//+
Characteristic Length {7, 6, 8, 5} = 0.2;
