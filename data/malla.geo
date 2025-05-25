Point(1) = {0, 0, 0};
Point(2) = {4, 0, 0};
Point(3) = {4, 3, 0};
Point(4) = {0, 3, 0};

//Line(1) = {2, 1};
//Line(2) = {3, 2};
//Line(3) = {4, 3};
//Line(4) = {1, 4};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Physical Line("Hola") = {1,2,3,4};
Physical Surface("Hola") = {1};

Transfinite Line{1, 3} = 2;
Transfinite Line{2, 4} = 2;

Transfinite Surface{1};

Mesh 2;

Mesh.MshFileVersion = 2.2;
Save "malla.msh";



