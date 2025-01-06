
a = 4.0;
b = 4.0;
n = 21;
m = 101;

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {  a, 0.0, 0.0};
Point(3) = {  a,   b, 0.0};
Point(4) = {0.0,   b, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

//Transfinite Curve{1,2,3,4} = n;
Transfinite Curve{2,4} = m;
Transfinite Curve{1,3} = n;

Transfinite Surface{1};
// Physical Curve("Γᵗ₁") = {2};
// Physical Curve("Γᵗ₂") = {3};
Physical Curve("Γ¹") = {1};
Physical Curve("Γ²") = {2};
Physical Curve("Γ³") = {3};
Physical Curve("Γ⁴") = {4};
Physical Surface("Ω") = {1};

Mesh.Algorithm = 9;
Mesh.MshFileVersion = 2;
Mesh 2;
RecombineMesh;
