function [zenith_homo, inlierId] = vp_ransac_refinement(lines_homo, opt) % [Zhai et al. 2016]
option = struct();
option.iterNum = 50;
option.thInlrRatio = .02;
option.thDist = sind(opt.theta_con);
[~, inlierId] = ransac_intersection(lines_homo, option);
if opt.debug_fileid > 0
  fprintf(opt.debug_fileid, "ransac inlier %d: ", size(inlierId, 2));
  for i = 1:size(inlierId, 2)
    fprintf(opt.debug_fileid, "%d ", inlierId(1, i) - 1);
  end
  fprintf(opt.debug_fileid, "\n");
end
zenith_homo = lines_normal(lines_homo(:,inlierId), opt);
