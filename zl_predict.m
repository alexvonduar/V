function zl = zl_predict(ls, dist_max, u0, v0, width, height, params);
hist = [];
ang = [];
for i = 1:size(ls,1)
    ang(i) = line_angle2(ls(i,:));
    if ang(i) < -pi/4
        ang(i) = ang(i) + pi;
    end
    if abs(ang(i)-pi/2) < params.theta_v
        l = line_hmg_from_two_points(ls(i,1:2),ls(i,3:4));
        dist = abs(dot(l, [width/2;height/2;1]));
        if dist < dist_max
            hist(end+1) = ang(i);
        end
    end
end
[N,edges0] = histcounts(hist,(pi/2-pi/8):pi/180:(pi/2+pi/8));
N(find(N <= 5)) = 0;
if sum(N) == 0
    N(round(length(N)/2)) = 1;
end
edges = (edges0(1:end-1)+edges0(2:end))/2;
[max_modes, H] = mnf_modes(N, 1);
if isempty(max_modes)
    max_modes = [1 size(N,2)];
    H = 0;
else
    [~,I] = sort(H, 'descend');
    H = H(I);
    max_modes = max_modes(I,:);
end
zl = [];
for i = 1:size(max_modes,1)
    Ni = zeros(1,size(N,2));
    a = max_modes(i,1);
    b = max_modes(i,2);
    Ni(a:b) = N(a:b);
    m = max(Ni);
    I = find(Ni == m);
    for j = I(1)
        a = edges(j);
        l = abs(v0/sin(a));
        zl(end+1) = u0+l*cos(pi-a);
    end
end
end
