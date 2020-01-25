//+ Parameter defaults
If(!(Exists(lc)))
	lc = DefineNumber[ 5, Name "lc" ];
EndIf

If(!(Exists(x0)))
	x0 = DefineNumber[ -1, Name "x0" ];
EndIf

//+
SetFactory("OpenCASCADE");
Rectangle(1) = {x0, 0, 0, 1, 1, 0};
//+
Disk(3) = {0, 0, 0, lc, lc};
//+
Coherence;
//+
Physical Surface("outer", 1) = {2};
//+
Physical Surface("inner", 2) = {1};
//+
Characteristic Length {5, 4, 2, 3} = 0.1;
//+
Physical Curve("outer", 1) = {1};
//+
Physical Curve("inner", 2) = {4, 5, 2, 3};
