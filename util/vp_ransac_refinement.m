function [zenith_homo, inlierId] = vp_ransac_refinement(lines_homo, opt) % [Zhai et al. 2016]
option = struct();
option.iterNum = 50;
option.thInlrRatio = .02;
option.thDist = sind(opt.theta_con);
[~, inlierId] = ransac_intersection(lines_homo, option);
zenith_homo = lines_normal(lines_homo(:,inlierId), opt);
