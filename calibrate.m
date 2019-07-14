%%-------------------------------------------------------------------------
%
% calibrate.m  
% 
% Copyright 2018 Gilles Simon <gsimon@loria.fr> Université de Lorraine
%
% Version 1.0 - September 5, 2018
%
% calibrate.m allows computing the camera focal length and finding a 
% Manhattan frame in a set of VPs. This function is explained in our 
% previous paper, section 4.3:
%
%   "Gilles Simon, Antoine Fond, Marie-Odile Berger.
%   A Simple and Effective Method to Detect Orthogonal Vanishing Points in
%   Uncalibrated Images of Man-Made Environments.
%   Eurographics 2016, May 2016, Lisbon, Portugal." 
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%-------------------------------------------------------------------------

function [focal, manh_vps, confident] = calibrate(zvp, hvp, width, height)
infty = 4;
accuracy = 2;
nt = 2^(5+min(accuracy, log2(width)-5));
step_h = width/nt;
Dt_h = atan(step_h/width);
u0 = width/2;
v0 = height/2;
manh_vps = [];
focal = -1;
confident = -1;
mht = [];
mht_count = 0;
mht_best = 0;
hvp_count = size(hvp,1);
if hvp_count > 1
    for v = 1:hvp_count
        for v2 = (v+1):hvp_count
            x1 = hvp(v,1) - u0;
            y1 = hvp(v,2) - v0;
            x2 = hvp(v2,1) - u0;
            y2 = hvp(v2,2) - v0;
            ff = sqrt(-(x1*x2+y1*y2));
            if ff*ff > 0 && ff > 0.28*width && ff < 3.8*width
                mht_count = mht_count + 1;
                vv1 = [x1/ff;y1/ff;1];
                vv2 = [x2/ff;y2/ff;1];
                vv3 = cross(vv1,vv2);
                mht(mht_count, 1) = ff*vv3(1)/vv3(3)+u0;
                mht(mht_count, 2) = ff*vv3(2)/vv3(3)+v0;
                mht(mht_count, 3) = ff;
                mht(mht_count, 5) = v;
                mht(mht_count, 6) = v2;
                at_vz1 = atan(mht(mht_count, 2)/double(height))/Dt_h;
                at_vz2 = atan(zvp(2)/double(height))/Dt_h;
                mht(mht_count, 4) = min(abs(at_vz1-at_vz2), abs(at_vz1+at_vz2));
            end
        end
    end
end
if mht_count > 0
    [vvpc_min, mht_best] = min(mht(:,4));
    if vvpc_min < 2 || abs(zvp(2)-v0) >= infty*height
        if vvpc_min < 2
            confident = 3;
        else
            [~, mht_best] = max(abs(mht(:,2)-v0));
            confident = 1;
        end
        id1 = mht(mht_best, 5);
        id2 = mht(mht_best, 6);
        x1 = hvp(id1,1);
        x2 = hvp(id2,1);
        y = mht(mht_best, 2);
        if (y > v0 && x1 > x2) || (y < v0 && x1 < x2)
            mht(mht_best, 5) = id2;
            mht(mht_best, 6) = id1;
        end
        focal = mht(mht_best,3);
        manh_vps = [hvp(mht(mht_best, 5),1) hvp(mht(mht_best, 5),2);hvp(mht(mht_best, 6),1) hvp(mht(mht_best, 6),2);zvp(1) zvp(2)];
    end
end
if confident == -1 && hvp_count > 0 && abs(zvp(2) - v0) < infty*height
    id_finite = find(abs(hvp(:,1)-u0) < infty*width);
    if length(id_finite) > 0
        y1 = hvp(id_finite(1),2) - v0;
        y2 = zvp(2) - v0;
        ff = sqrt(-y1*y2);
        if ff*ff > 0 && ff > 0.28*width && ff < 3.8*width
            focal = ff;
            K = [[focal 0 u0];[0 focal v0];[0 0 1]];
            r2 = [hvp(id_finite(1),1)-u0; hvp(id_finite(1),2)-v0; focal];
            r2 = r2/norm(r2);
            r3 = [zvp(1)-u0; zvp(2)-v0; ff];
            r3 = r3/norm(r3);
            vpx = K*cross(r2,r3);
            vpx = vpx/vpx(3);
            manh_vps = [vpx(1)+u0 vpx(2)+v0; hvp(id_finite(1),1) hvp(id_finite(1),2); zvp(1) zvp(2)];
            confident = 2;
        end
    end
end
end