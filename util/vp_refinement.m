function vps_homo = vp_refinement(lines_homo, vps_homo, horizon_homo, params) % from [Zhai et al.]
rng(1)
niters = params.refine_niters;
inlier_thr = sind(params.theta_con);
ngrps = size(vps_homo,2);
for ig = 1:ngrps
  for it = 1:niters
    good_ids = find(abs(vps_homo(:,ig)'*lines_homo) < inlier_thr)';
    vps_homo(:,ig) = lines_normal(lines_homo(:,good_ids), params, horizon_homo);
  end
end
