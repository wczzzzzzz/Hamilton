
a = 210;
b = 240;
n = 10;

Point(1) = { -a,  -a, 0.0, a/n};
Point(2) = {  a,  -a, 0.0, a/n};
Point(3) = {  a,   a, 0.0, a/n};
Point(4) = { -a,   a, 0.0, a/n};
Point(5) = {  0,  -a, 0.0, a/n};
Point(6) = { -b,  -b, 0.0};
Point(7) = {  b,  -b, 0.0};
Point(8) = {  b,   b, 0.0};
Point(9) = { -b,   b, 0.0};
Point(10) = {  0,  -b, 0.0};

Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,1};
Line(6) = {6,10};
Line(7) = {10,7};
Line(8) = {7,8};
Line(9) = {8,9};
Line(10) = {9,6};

Curve Loop(1) = {1,2,3,4,5};
Curve Loop(2) = {6,7,8,9,10};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Curve("Γ") = {3,4,5,8,9,10};
Physical Surface("Ω") = {1,2};
Physical Point("Γᵗ") = {5,10};

// Transfinite Curve{1,2} = (n+1)/2;
// Transfinite Curve{3,4,5} = n;
// Transfinite Surface{1};
Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;