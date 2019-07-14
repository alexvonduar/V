function l = line_hmg_from_two_points(p1,p2)
v1 = double([p1(1);p1(2);1]);
v2 = double([p2(1);p2(2);1]);
l=cross(v1,v2);
l=l/sqrt(double(l(1)^2+l(2)^2));