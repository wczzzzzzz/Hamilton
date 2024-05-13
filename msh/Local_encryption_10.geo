
a = 4.0;
b = 4.0;
c1 = 0.02;
c2 = 0.5;

Point(1) = {0.0, 0.0, 0.0, c1};
Point(2) = {  a, 0, 0.0, c2};
Point(3) = {  a,  b, 0.0, c1};
Point(4) = {0.0,  b, 0.0, c2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,3};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Curve("Γ₁") = {1};
Physical Curve("Γ₂") = {2};
Physical Curve("Γ₃") = {3};
Physical Curve("Γ₄") = {4};
Physical Surface("Ω") = {1};

Curve{5} In Surface{1};
Field[1] = Distance;
Field[1].CurvesList = [5];
Field[1].Sampling = 100;

Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;
