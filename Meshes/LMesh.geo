Merge "LMesh.msh";
Line(9) = {9, 2};
//+
Line(10) = {2, 1};
//+
Line(11) = {1, 3};
//+
Line(12) = {3, 4};
//+
Line(13) = {4, 5};
//+
Line(14) = {5, 7};
//+
Line(15) = {7, 8};
//+
Line(16) = {8, 9};
//+
Curve Loop(1) = {15, 16, 9, 10, 11, 12, 13, 14};
//+
Plane Surface(2) = {1};
