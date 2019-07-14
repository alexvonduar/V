% convert lines or points from cartesian system (image frame) to homogenous system
% according to [Zhai et al., CVPR'2016]
function geometry_homo = normalize(geometry_img, width, height, focal)
center = [width, height] / 2;
if size(geometry_img, 2) == 2
  geometry_normalized = bsxfun(@minus, geometry_img, center) / focal;
else
  geometry_normalized = bsxfun(@minus, geometry_img, [center, center]) / focal;
end
N = size(geometry_normalized, 1);
if size(geometry_normalized, 2) == 2 
  geometry_homo = [geometry_normalized ones(N, 1)]';
%   geometry_homo = bsxfun(@times, geometry_homo, sign(geometry_homo(3,:) + eps));
else  
  x1 = normalize(geometry_img(:,[1,2]), width, height, focal);
  x2 = normalize(geometry_img(:,[3,4]), width, height, focal);
  geometry_homo = cross(x1, x2);
  geometry_homo = bsxfun(@times, geometry_homo, sign(geometry_homo(2,:) + eps));
end
geometry_homo = bsxfun(@rdivide, geometry_homo, sqrt(sum(geometry_homo.^2, 1)));