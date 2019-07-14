% convert lines or points from homogeneous system to cartesian system (image frame)
% according to [Zhai et al., CVPR'2016]
function geometry_img = unnormalize(geometry_homo, width, height, focal, isline)
if nargin == 5 && ~isline  
  p = geometry_homo;
  geometry_img = bsxfun(@rdivide, p(1:2, :), p(3,:))';
  geometry_img = bsxfun(@plus, (geometry_img * focal), [width, height]/2);
else
  h = geometry_homo;
  pl = (0 - width/2) / focal;
  pr = (width - width/2) / focal;
  ly = (height/2 -((h(1,:)*pl + h(3,:))./h(2,:)) * focal)';
  ry = (height/2 -((h(1,:)*pr + h(3,:))./h(2,:)) * focal)';
  geometry_img = [zeros(size(ly)), ly, width*ones(size(ry)), ry];
  geometry_img = line_to_segment(geometry_img);
end