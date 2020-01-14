SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.05};
//+
Point(2) = {0, 0.5, 0, 0.05};
//+
Point(3) = {0.5, 0.5, 0, 0.05};
//+
Point(4) = {-0.5, 0.5, 0, 0.05};
//+
Point(5) = {-0.5, -0.5, 0, 0.05};
//+
Point(6) = {0.5, -0.5, 0, 0.05};
//+
Point(7) = {0, -0.5, 0, 0.05};


//+
Line(1) = {7, 6};
//+
Line(2) = {6, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {1, 7};
//+
Line(6) = {7, 5};
//+
Line(7) = {5, 4};
//+
Line(8) = {4, 2};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+

//+
Point(8) = {0, 1, 0, 0.05};
//+
Point(9) = {1, 0, 0, 0.05};
//+
Point(10) = {0, -1, 0, 0.05};
//+
Point(11) = {-1, 0, 0, 0.05};
//+
Circle(9) = {8, 1, 9};
//+
Circle(10) = {9, 1, 10};
//+
Circle(11) = {10, 1, 11};
//+
Circle(12) = {11, 1, 8};
//+
Line(13) = {8, 2};
//+
Line(14) = {7, 10};
//+
Curve Loop(3) = {12, 13, -8, -7, -6, 14, 11};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 10, -14, 1, 2, 3, -13};
//+
Plane Surface(4) = {4};
//+
Characteristic Length {8, 9, 10, 11} = 0.2;
//+
Physical Curve("inner", 2) = {3, 2, 1, 6, 7, 8};
//+
Physical Curve("outer", 1) = {9, 10, 11, 12};
//+
Physical Surface("inner", 2) = {2, 1};
//+
Physical Surface("outer", 1) = {3, 4};