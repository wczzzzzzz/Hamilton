
L = 4;
n = 4;
a = 0.9;

Point(1) = {0, 0, 0, a};
Point(2) = {L, 0, 0, a};
Point(3) = {L, L, 0, a};
Point(4) = {0, L, 0, a};
Point(5) = {0, 0, L, a};
Point(6) = {L, 0, L, a};
Point(7) = {L, L, L, a};
Point(8) = {0, L, L, a};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {1, 10, -5, -9};
Curve Loop(4) = {2, 11, -6, -10};
Curve Loop(5) = {3, 12, -7, -11};
Curve Loop(6) = {4, 9, -8, -12}; 

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Physical Surface("Γ¹") = {2};
Physical Surface("Γ²") = {1, 3, 4, 5, 6};
Physical Surface("Γ³") = {};
Physical Surface("Γ⁴") = {};
Physical Surface("Γ⁵") = {};
Physical Surface("Γ⁶") = {};
Physical Volume("Ω") = {1};

Transfinite Curve{1:12} = n+1;
Transfinite Surface{1:6};
Transfinite Volume{1}; 

//Recombine Surface{1:6};
//Recombine Volume{1};

Mesh 3;