//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Box(2) = {0.2, 0.7, 0.2, 0.8, 0.2, 0.6};
//+
Coherence;

//+
Physical Surface("outer", 2) = {13, 15, 16, 14, 17};
//+
Physical Surface("symmetry", 1) = {18};
//+
Physical Surface("symmetry forced", 3) = {8};
//+
Physical Volume("inner", 2) = {2};
//+
Physical Volume("outer", 1) = {3};
//+
Characteristic Length {19, 24, 23, 22, 21, 17, 18, 20} = 0.2;
//+
Characteristic Length {9, 10, 12, 11, 15, 13, 14, 16} = 0.05;
