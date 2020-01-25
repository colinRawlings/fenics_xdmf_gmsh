//+ Parameter defaults
If(!(Exists(radius)))
	radius = DefineNumber[ 5, Name "radius" ];
EndIf

If(!(Exists(x0)))
	x0 = DefineNumber[ -0.5, Name "x0" ];
EndIf

If(!(Exists(y0)))
	y0 = DefineNumber[ 1, Name "y0" ];
EndIf

If(!(Exists(dx_inner_mesh)))
	dx_inner_mesh = DefineNumber[ 0.1, Name "dx_inner_mesh" ];
EndIf

//+
Point(1) = {x0, y0, 0, 1.0};
//+
Point(2) = {x0+1, y0, 0, 1.0};
//+
Point(3) = {x0+1, y0+1, 0, 1.0};
//+
Point(4) = {x0, y0+1, 0, 1.0};
//+
Point(5) = {radius, 0, 0, 1.0};
//+
Point(6) = {0, radius, 0, 1.0};
//+
Point(7) = {-radius, 0, 0, 1.0};
//+
Point(8) = {0, -radius, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Point(9) = {0, 0, 0, 1.0};
//+
Circle(5) = {8, 9, 5};
//+
Circle(6) = {5, 9, 6};
//+
Circle(7) = {6, 9, 7};
//+
Circle(8) = {7, 9, 8};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(2) = {1, 2};

//

Physical Surface("outer", 1) = {2};
//+
Physical Surface("inner", 2) = {1};
//+
Characteristic Length {1, 2, 3, 4} = dx_inner_mesh;
//+
Physical Curve("outer", 1) = {5, 6, 7, 8};
//+
Physical Curve("inner", 2) = {1, 2, 3, 4};
