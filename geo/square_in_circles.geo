

If(!(Exists(radius_1)))
	radius_1 = DefineNumber[ 1, Name "radius_1" ];
EndIf

If(!(Exists(radius_2)))
    radius_2 = DefineNumber[ 5, Name "radius_2" ];
EndIf

If(!(Exists(square_side)))
    square_side = DefineNumber[ 0.5, Name "square_side" ];
EndIf

If(!(Exists(dx_circle_1_mesh)))
	dx_circle_1_mesh = DefineNumber[ 0.2, Name "dx_circle_1_mesh" ];
EndIf

//+
// SetFactory("Built-in");
//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-0.5*square_side, -0.5*square_side, 0, square_side, square_side, 0.001};
//+
Circle(9) = {0, 0, 0, radius_1, 0, 2*Pi};
//+
Circle(10) = {0, 0, 0, radius_2, 0, 2*Pi};
//+
Curve Loop(2) = {9};
//+
Curve Loop(3) = {2, 3, 4, 5, 6, 7, 8, 1};
//+
Plane Surface(2) = {2, 3};
//+
Curve Loop(4) = {10};
//+
Curve Loop(5) = {9};
//+
Plane Surface(3) = {4, 5};
//+
Characteristic Length {2, 3, 4, 5, 6, 7, 8, 1} = 0.1;
//+
Characteristic Length {9} = dx_circle_1_mesh;
//+
Physical Curve("square_1", 1) = {2, 3, 4, 5, 6, 7, 8, 1};
//+
Physical Curve("circle_1", 2) = {9};
//+
Physical Curve("circle_2", 3) = {10};
//+
Physical Surface("square_1", 1) = {1};
//+
Physical Surface("circle_1", 2) = {2};
//+
Physical Surface("circle_2", 3) = {3};
