function [modes_homo, modes_offset, modes_left, modes_right, H] = hl_predict(seglines, zenith_homo, u0, v0, width, height, focal, params)
L = params.L_h;
x2 = zenith_homo(1)/zenith_homo(3);
y2 = zenith_homo(2)/zenith_homo(3);
tilt = atan(x2/y2);
if isnan(tilt)
    tilt = 0;
end
offsets = [];
for i = 1:size(seglines,1)
    v = [seglines(i,3)-seglines(i,1) seglines(i,4)-seglines(i,2)];
    scal = dot(v/norm(v),[cos(-tilt), sin(-tilt)]);
    ang = acos(scal);
    if abs(ang) < params.theta_h*pi/180
        offsets(end+1) = dot([(seglines(i,3)+seglines(i,1))/2-u0 (seglines(i,4)+seglines(i,2))/2-v0],[cos(-tilt+pi/2) sin(-tilt+pi/2)]);
    end
end
if isempty(offsets)
    offsets = -height/2:height/2;
end
[maxValue, ~] = max(offsets);
[minValue, ~] = min(offsets);
edges = linspace(minValue, maxValue, L + 1)

%[N,edges] = histcounts(offsets,L);
N = histcounts(offsets, edges);
[max_modes, H] = mnf_modes(N, 1);
if isempty(max_modes)
    max_modes(1,:) = [1 size(N,2)];
    H(1) = -1;
else
    [~,I] = sort(H, 'descend');
    H = H(I);
    max_modes = max_modes(I,:);
end
nmodes = size(max_modes,1);
modes_offset = [];
modes_left = [];
modes_right = [];
modes_homo = [];
for i = 1:nmodes
    Ni = zeros(1,size(N,2));
    a = max_modes(i,1);
    b = max_modes(i,2);
    Ni(a:b) = N(a:b);
    [~,bin] = max(Ni);
    modes_offset(end+1) = edges(1)+bin/L*(edges(L)-edges(1));
    mnf_center = [u0+modes_offset(end)*cos(-tilt+pi/2) v0+modes_offset(end)*sin(-tilt+pi/2)];
    hmnf = line_hmg_from_two_points(mnf_center, mnf_center+[cos(-tilt) sin(-tilt)]);
    modes_left(end+1) = -hmnf(3)/hmnf(2);
    modes_right(end+1) = (-hmnf(1)*width-hmnf(3))/hmnf(2);
    mode_seg = [0, modes_left(1), width, modes_right(1)];
    modes_homo(1:3,end+1) = normalize(mode_seg, width, height, focal);
end
end