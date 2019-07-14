function [samp_homo, samp_left, samp_right] = hl_sample(zenith_homo, modes_homo, modes_offset, modes_left, modes_right, H, u0, v0, width, height, focal, params)
rng(1) % fix random seed
samp_homo = modes_homo;
samp_left = modes_left;
samp_right = modes_right;
S = params.S;
tilt = atan((zenith_homo(1)/zenith_homo(3))/(zenith_homo(2)/zenith_homo(3)));
if isnan(tilt)
    tilt = 0;
end
nsamp = ceil(S/size(modes_homo,2));
x = linspace(-2,2,1e5);
uniformrnd = cumsum(ones(size(x)));
for i = 1:size(modes_homo,2)
    for j = 1:(nsamp-1)
        if H(i) <= 0
            draw = rand()*uniformrnd(end);
            id = min(find(uniformrnd >= draw));
            orand = x(id)*height + modes_offset(i);
        else
            orand = normrnd(0,params.sigma*height) + modes_offset(i);
        end
        mnf_center = [u0+orand*cos(-tilt+pi/2) v0+orand*sin(-tilt+pi/2)];
        hmnf = line_hmg_from_two_points(mnf_center, mnf_center+[cos(-tilt) sin(-tilt)]);
        samp_left(end+1) = -hmnf(3)/hmnf(2);
        samp_right(end+1) = (-hmnf(1)*width-hmnf(3))/hmnf(2);
        mode_seg = [0, samp_left(end), width, samp_right(end)];
        samp_homo(1:3,end+1) = normalize(mode_seg, width, height, focal);
    end
end
end
