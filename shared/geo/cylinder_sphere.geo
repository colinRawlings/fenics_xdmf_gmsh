
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
If(!(Exists(outer_radius)))
    outer_radius = DefineNumber[ 2, Name "outer_radius" ];
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
Cylinder(1) = {0, 0, -0.5*cylinder_height, 0, 0, cylinder_height, cylinder_radius, 2*Pi};
//+
Sphere(2) = {0, 0, 0, sphere_radius, -Pi/2, Pi/2, 2*Pi};
//+
Coherence;
//+
Characteristic Length {3, 4} = dx_mesh_outer;
//+
Characteristic Length {1, 2} = dx_mesh_inner;
//+
Physical Surface("outer", 1) = {4};
//+
Physical Surface("inner", 2) = {3, 1, 2};
//+
Physical Volume("inner", 2) = {1};
//+
Physical Volume("outer", 1) = {2};
