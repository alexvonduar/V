function [imR, imref] = orthorectify(im, H, line, s)
max_stretch = 16;
step = round(size(im, 2)/64);
width = size(im,2);
height = size(im,1);
umin = 1e10;
umax = -1e10;
vmin = 1e10;
vmax = -1e10;
imR = [];
imref = [];
for i = 1:step:width
    for j = 1:step:height
        if dot(line, [i; j; 1])*(-1)^s > 0
            pix = [i,i+1;j,j+1;1,1];
            Tpix = H*pix;
            Tpix(:,1) = Tpix(:,1)/Tpix(3,1);
            Tpix(:,2) = Tpix(:,2)/Tpix(3,2);
            u = Tpix(1,1);
            v = Tpix(2,1);
            u2 = Tpix(1,2);
            v2 = Tpix(2,2);
            if abs(u2-u) < max_stretch && abs(v2-v) < max_stretch
%                 mask(j,i) = 1;
                if u < umin umin = u; else if u > umax umax = u; end; end;
                if v < vmin vmin = v; else if v > vmax vmax = v; end; end;
            end
        end
    end
end
heightT = floor(vmax-vmin);
widthT = floor(umax-umin);
if heightT > 0 && widthT > 0
    imref = imref2d([heightT,widthT],[umin,umax],[vmin,vmax]);
    tform = projective2d(H');
    imR = imwarp(im,tform,'OutputView',imref);
end
end
