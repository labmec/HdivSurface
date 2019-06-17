// Gmsh project created on Tue Apr 30 14:17:43 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 0, 1, 1.0};
//+
Point(2) = {0, 1, 1, 1.0};
//+
Point(3) = {0, 0, 1, 1.0};
//+
Point(4) = {1, 1, 1, 1.0};
//+
Line(1) = {2, 4};
//+
Line(2) = {4, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 2};
//+
Physical Curve("boundary") = {1, 2, 3, 4};
//+
Physical Surface("domain") = {1};
