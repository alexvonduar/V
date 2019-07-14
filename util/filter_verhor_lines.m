function lines_id = filter_verhor_lines(ls_homo, z_homo,  params) % as in [Zhai et al. 2016]

% filter vertical line segments 
cos_val = abs(ls_homo'*z_homo);
inlier_id = cos_val > sind(params.theta_verline);
lines_id = find(inlier_id);

% filter horizontal line segments 
if ~params.include_infinite_hvps
    ortho_thres = sind(params.theta_horline);
    lhomo = ls_homo(:,lines_id);
    cos_val = abs(lhomo'*[-z_homo(2); z_homo(1); 0]);
    inlier_id = cos_val > ortho_thres;
    lines_id = lines_id(inlier_id);
end