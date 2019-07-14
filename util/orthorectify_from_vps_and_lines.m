%%-------------------------------------------------------------------------
%
% orthorectify_from_vps_and_lines.m  
% 
% Copyright 2018 Gilles Simon <gsimon@loria.fr> Université de Lorraine
%
% Version 1.0 - September 5, 2018
%
% orthorectify_from_vps_and_lines.m allows warping (rectifying) an image so 
% that all the vertical planes present in this image appear as in a frontal 
% view. Please see our webpage https://members.loria.fr/GSimon/software/v/
% to see some examples.
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

function [imR,maskR,transform] = orthorectify_from_vps_and_lines(im, hvps, hvp_groups, zenith, z_group, lines, n_lines_min, K, horizon_line, PLOT)
n_hvp = size(hvps, 1);
width = size(im, 2);
height = size(im, 1);
n_imR = 0;
imR = cell(1);
maskR = cell(1);
transform = cell(1);
vp_association = -1*ones(length(lines));
vp_association(z_group) = 0;
for i = 1:numel(hvp_groups)
    vp_association(hvp_groups{i}) = i;
end
if ~isempty(K)
    Ki = inv(K);
    pp = [K(1,3) K(2,3)];
    pp_zen_line = line_hmg_from_two_points(pp, zenith);
    vp_zen = [zenith(1,1); zenith(1,2); 1];
    y = Ki*vp_zen;
    y = y/norm(y);
    if dot(horizon_line, vp_zen) < 0
        y = -y;
    end
end
for i = 1:n_hvp
    if PLOT
        figure(30000+i); imshow(im); hold on;
    end
    imR{i} = [];
    maskR{i} = [];
    transform{i} = [];
    hvp = [hvps(i,1) hvps(i,2)];
    vp_zen_line = line_hmg_from_two_points(hvp, zenith);
    if vp_zen_line(1) < 0
        vp_zen_line = -vp_zen_line;
    end
    zen_line_normal = vp_zen_line(1:2);
    zen_line_normal = zen_line_normal / norm(zen_line_normal);
    proj_min(1) = 1;
    proj_max(1) = -1;
    proj_min(2) = 1;
    proj_max(2) = -1;
    n_lines_zen(1) = 0;
    n_lines_zen(2) = 0;
    for j = 1:size(lines,1)
        centroid(j,1) = (lines(j,1)+lines(j,3))/2;
        centroid(j,2) = (lines(j,2)+lines(j,4))/2;
        lines_dir(j,1) = lines(j,1)-lines(j,3);
        lines_dir(j,2) = lines(j,2)-lines(j,4);
        lines_dir(j,:) = lines_dir(j,:)/norm(lines_dir(j,:));
        bundle_dir = centroid(j,:)'-zenith(1,1:2)';
        bundle_dir = bundle_dir / norm(bundle_dir);
        proj = dot(zen_line_normal, bundle_dir);
        if vp_association(j) == 0
            for s = 1:2
                if dot(zen_line_normal,bundle_dir)*(-1)^s > 0 
                    n_lines_zen(s) = n_lines_zen(s)+1;
                    if proj > proj_max(s)
                        proj_max(s) = proj;
                        idmax(s) = j;
                    end
                    if proj < proj_min(s)
                        proj_min(s) = proj;
                        idmin(s) = j;
                    end
                end
            end
        end
    end
    for s=1:2
        if n_lines_zen(s) >= n_lines_min
            centroid_min = centroid(idmin(s),1:2);
            centroid_max = centroid(idmax(s),1:2);
            dist_min = line_size([hvp(1), hvp(2), centroid_min(1), centroid_min(2)]);
            dist_max = line_size([hvp(1), hvp(2), centroid_max(1), centroid_max(2)]);
            if dist_min < dist_max
                middle(1) = (hvp(1)+centroid_max(1))/2;
                middle(2) = (hvp(2)+centroid_max(2))/2;
                dist_middle = line_size([hvp, middle]);
                if dist_min < dist_middle && dist_middle < dist_max
                    centroid_min = middle;
                end
            else
                middle(1) = (hvp(1)+centroid_min(1))/2;
                middle(2) = (hvp(2)+centroid_min(2))/2;
                dist_middle = line_size([hvp, middle]);
                if dist_max < dist_middle && dist_middle < dist_min
                    centroid_max = middle;
                end
            end
            lzmin = line_hmg_from_two_points(zenith(1,1:2), centroid_min);
            lzmax = line_hmg_from_two_points(zenith(1,1:2), centroid_max);
            lzmin_V(1) = centroid_min(1)-zenith(1);
            lzmin_V(2) = centroid_min(2)-zenith(2);
            lzmin_V = lzmin_V/norm(lzmin_V);
            lzmax_V(1) = centroid_max(1)-zenith(1);
            lzmax_V(2) = centroid_max(2)-zenith(2);
            lzmax_V = lzmax_V/norm(lzmax_V);
            hvp_amin = 2*pi;
            hvp_amax = -2*pi;
            n_lines_hvp = 0;
            for j = 1:size(lines,1)
                if vp_association(j) == i
                    if dot(vp_zen_line, [lines(j,1); lines(j,2); 1])*(-1)^s > 0 && dot(vp_zen_line, [lines(j,3); lines(j,4); 1])*(-1)^s > 0 && abs(dot(lines_dir(j,:)',lzmin_V)) < cos(pi/8) && abs(dot(lines_dir(j,:)',lzmax_V)) < cos(pi/8)
                        n_lines_hvp = n_lines_hvp + 1;
                        bundle_line = [hvps(i,1) hvps(i,2) centroid(j,1) centroid(j,2)];
                        if PLOT
                            line_plot(bundle_line, 'g:');
                            line_plot(lines(j,:), 'g');
                            drawnow();
                        end
                        a = -line_angle2(bundle_line);
                        if a < hvp_amin
                            hvp_amin = a;
                            idmin2 = j;
                        else
                            if a > hvp_amax
                                hvp_amax = a;
                                idmax2 = j;
                            end
                        end
                    end
                end
            end
            if n_lines_hvp >= n_lines_min
                n_imR = n_imR + 1;
                transform{n_imR}.K = K;
                if length(K) > 0
                    vp = [hvps(i,1); hvps(i,2); 1];
                    x = Ki*vp;
                    x = x/norm(x);
                    if (dot(pp_zen_line, vp) < 0 && dot(horizon_line, vp_zen) < 0) || (dot(pp_zen_line, vp) > 0 && dot(horizon_line, vp_zen) > 0)
                        x = -x;
                    end
                    if (dot(pp_zen_line, vp) < 0 && dot(horizon_line, vp_zen) < 0 && s == 1) || (dot(pp_zen_line, vp) > 0 && dot(horizon_line, vp_zen) < 0  && s == 2) || (dot(pp_zen_line, vp) > 0 && dot(horizon_line, vp_zen) > 0 && s == 1) || (dot(pp_zen_line, vp) < 0 && dot(horizon_line, vp_zen) > 0  && s == 2)
                        x = -x; 
                    end
                    z = cross(x,y);
                    transform{n_imR}.R = [x, y, z];
                    transform{n_imR}.H = K*inv(transform{n_imR}.R)*Ki;
                else
                    lhmin = line_hmg_from_two_points(hvps(i,1:2), centroid(idmin2,1:2));
                    lhmax = line_hmg_from_two_points(hvps(i,1:2), centroid(idmax2,1:2));
                    c(1,1:2) = line_hmg_intersect(lhmin, lzmin);
                    c(2,1:2) = line_hmg_intersect(lhmax, lzmin);
                    c(3,1:2) = line_hmg_intersect(lhmax, lzmax);
                    c(4,1:2) = line_hmg_intersect(lhmin, lzmax);
                    [sc, id] = sort(c(:,1));
                    idul = id(2-(c(id(1),2) < c(id(2),2)));
                    idll = id(2-(c(id(1),2) > c(id(2),2)));
                    up = (idll == mod(idul,4)+1);
                    idc = idul;
                    for j=1:4
                        corners(1:2,j) = c(idc,:)';
                        if up
                            idc = mod(idc,4)+1;
                        else
                            idc = idc-1;
                            if idc == 0 idc = 4; end
                        end
                    end
                    if PLOT
                        line_plot([zenith(1,1) zenith(1,2) centroid_min(1) centroid_min(2)],'r--');
                        line_plot([zenith(1,1) zenith(1,2) centroid_max(1) centroid_max(2)],'r--');
                        line_plot([hvps(i,1) hvps(i,2) centroid(idmin2,1) centroid(idmin2,2)],'r--');
                        line_plot([hvps(i,1) hvps(i,2) centroid(idmax2,1) centroid(idmax2,2)],'r--');
                        lhmin = line_hmg_from_two_points(hvps(i,1:2), centroid(idmin2,1:2));
                        lhmax = line_hmg_from_two_points(hvps(i,1:2), centroid(idmax2,1:2));
                        plot([corners(1,:),corners(1,1)],[corners(2,:),corners(2,1)],'r');
                        drawnow();
                    end
                    BW = poly2mask(corners(1,:), corners(2,:), height, width);
                    M = sum(corners,2)/4;
                    A = sum(BW(:));
                    AR = (sqrt(sum((corners(:,1)-corners(:,2)).^2))+sqrt(sum((corners(:,3)-corners(:,4)).^2)))  /  (sqrt(sum((corners(:,2)-corners(:,3)).^2))+sqrt(sum((corners(:,4)-corners(:,1)).^2)));
                    w = sqrt(A/AR);
                    h = AR*w;
                    pointsR = ([M(1)-w/2,M(1)-w/2,M(1)+w/2,M(1)+w/2;M(2)-h/2,M(2)+h/2,M(2)+h/2,M(2)-h/2]');
                    points = (corners(:,1:4)');
                    tform = estimateGeometricTransform(points,pointsR,'projective');
                    transform{n_imR}.H=tform.T';
                    transform{n_imR}.R = [];
                end
                [imR{n_imR},transform{n_imR}.imref] = orthorectify(im, transform{n_imR}.H, vp_zen_line, s);
                if ~isempty(transform{n_imR}.imref)
                    imwhite = ones(size(im,1),size(im,2));
                    tform = projective2d(transform{n_imR}.H');
                    maskR{n_imR}  = imwarp(imwhite,tform,'OutputView',transform{n_imR}.imref,'Interp','nearest');
                end
            end
        end
    end
end
end