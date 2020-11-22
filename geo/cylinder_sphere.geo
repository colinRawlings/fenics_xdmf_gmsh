
SetFactory("OpenCASCADE");
//+
If(!(Exists(cylinder_radius)))
    cylinder_radius = DefineNumber[ 0.5, Name "cylinder_radius" ];
EndIf
//+
If(!(Exists(cylinder_height)))
    cylinder_height = DefineNumber[ 0.5, Name "cylinder_height" ];
EndIf
//+
If(!(Exists(sphere_radius)))
    sphere_radius = DefineNumber[ 2, Name "sphere_radius" ];
EndIf
//+
If(!(Exists(dx_mesh_outer)))
    dx_mesh_outer = DefineNumber[ 0.5, Name "dx_mesh_outer" ];
EndIf
//+
If(!(Exists(dx_mesh_inner)))
    dx_mesh_inner = DefineNumber[ 0.2, Name "dx_mesh_inner" ];
EndIf
//+
Cylinder(1) = {0, 0, 0, 0, 0, 0.5*cylinder_height, cylinder_radius, 2*Pi};
//+
Cylinder(2) = {0, 0, -0.5*cylinder_height, 0, 0, 0.5*cylinder_height, cylinder_radius, 2*Pi};
//+
Sphere(3) = {0, 0, 0, sphere_radius, -Pi/2, Pi/2, 2*Pi};
//+
Coherence;

//+
Physical Volume("outer", 1) = {3};
//+
Physical Volume("inner", 2) = {1, 2};
//+
Physical Surface("inner", 2) = {3, 2, 4, 5};
//+
Physical Surface("outer", 1) = {1};
//+
Characteristic Length {1, 2} = dx_mesh_outer;
//+
Characteristic Length {3, 4, 5} = dx_mesh_inner;
