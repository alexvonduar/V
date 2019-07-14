function angle = line_angle2(line)
x1 = line(1,1);
y1 = line(1,2);
x2 = line(1,3);
y2 = line(1,4);
if x2 > x1
    angle = atan2((y2-y1),(x2-x1));
else
    angle = atan2((y1-y2),(x1-x2));
end
