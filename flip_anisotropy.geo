Point(1) = {-20000, 0, 0, 400};
Point(2) = {20000, 0, 0, 400};
Point(3) = {10000, 0, 0, 400};
Point(4) = {-10000, 0, 0, 400};
Point(5) = {-4000, 0, 0, 100};
Point(6) = {4000, 0, 0, 100};
Point(7) = {0, 0, 0, 50};
Point(8) = {-60000, 0, 0, 12500};
Point(9) = {60000, 0, 0, 12500};

Point(11) = {60000, 129000, 0, 25500};
Point(12) = {-60000, 129000, 0, 25500};
Point(13) = {-4000, 129000, 0, 800};
Point(14) = {4000, 129000, 0, 800};

Line(1) = {8, 1};
Line(2) = {1, 4};
Line(3) = {4, 5};
Line(4) = {5, 7};
Line(5) = {7, 6};
Line(6) = {6, 3};
Line(7) = {3, 2};
Line(8) = {2, 9};
Line(9) = {8, 12};
Line(10) = {12, 13};
Line(11) = {13, 14};
Line(12) = {14, 11};
Line(13) = {11, 9};
Line(14) = {5, 13};
Line(15) = {6, 14};
Line Loop(16) = {1, 2, 3, 14, -10, -9};
Plane Surface(17) = {16};
Line Loop(18) = {5, 15, -11, -14, 4};
Plane Surface(19) = {18};
Line Loop(20) = {7, 8, -13, -12, -15, 6};
Plane Surface(21) = {20};




Point(15) = {-60000, -491000, 0, 102000.0};
Point(16) = {60000, -491000, 0, 102000.0};
Line(22) = {8, 15};
Line(23) = {15, 16};
Line(24) = {16, 9};
Line Loop(25) = {23, 24, -8, -7, -6, -5, -4, -3, -2, -1, 22};
Plane Surface(26) = {25};
//taggings!!
Physical Surface("left") = {17};
Physical Surface("middle") = {19};
Physical Surface("right") = {21};
Physical Surface("air_layer") = {26};




Transfinite Line {4, 5, 3, 6} = 1000 Using Progression 1;
Transfinite Line {1, 8} = 100 Using Progression 1;
Transfinite Line {2, 7} = 1200 Using Progression 1;
