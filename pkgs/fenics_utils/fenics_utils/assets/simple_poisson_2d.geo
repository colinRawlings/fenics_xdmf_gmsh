//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Rectangle(2) = {0.2, 0.7, 0, 0.8, 0.2, 0};
//+
Coherence;
//+
Coherence;
//+
Physical Surface("outer", 1) = {3};
//+
Physical Surface("inner", 2) = {2};
//+
Physical Curve("symmetry", 1) = {7, 3};
//+
Physical Curve("symmetry forced", 3) = {9};
//+
Physical Curve("outer", 2) = {2, 1, 8};
//+
Characteristic Length {5, 6, 7, 4} = 0.05;
//+
Characteristic Length {3, 1, 2, 8} = 0.2;
