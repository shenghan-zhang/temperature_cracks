h = 0.1;

Point(1) = { 0, 0, 0, h};
Point(2) = {1.0, 0, 0, h};
Point(3) = {1.0, 0.2, 0, h};
Point(4) = { 0, 0.2, 0, h};

Point(5) = { 1.0, 0.3, 0, h};
Point(6) = {0.0, 0.3, 0, h};
Point(7) = {1.0, 0.5, 0, h};
Point(8) = {0.0, 0.5, 0, h};//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {3, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 4};
//+
Line(8) = {5, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 6};
//+
Line Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {6, 7, -4, 5};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {8, 9, 10, -6};
//+
Plane Surface(3) = {3};
//+
Physical Surface("cathode") = {1};
//+
Physical Surface("ceramic") = {2};
//+
Physical Surface("anode") = {3};
