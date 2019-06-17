// Gmsh project created on Fri May  3 10:24:38 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 1, 1.0};
//+
Point(2) = {1, 0, 1, 1.0};
//+
Point(3) = {0, 1, 1, 1.0};
//+
Point(4) = {1, 1, 1, 1.0};
//+
Line(1) = {3, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 3};
//+
Physical Curve("boundary") = {1, 2, 4, 3};
//+
Physical Surface("domain") = {1};
