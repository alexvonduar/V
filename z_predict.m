function [zenith_homo, zengroupIds] = z_predict(lines_homo, zenith_line_homo, params, refine)
%% threshold zenith LSs ids

lines_tilt_ortho = atand(-lines_homo(2,:) ./ lines_homo(1,:)); 
zenith_tilt = atand(-zenith_line_homo(2,1) ./ zenith_line_homo(1,1)); 
verticalInd = find(abs(lines_tilt_ortho - zenith_tilt) < params.theta_z);
zenith_homo_pred = lines_normal(lines_homo(:,verticalInd));
   
%% refinement

if refine
    [zenith_homo, ~] = vp_ransac_refinement(lines_homo(:,verticalInd), params);
else
    zenith_homo = zenith_homo_pred;
end

%% LS grouping

[~, groups] = vp_score(zenith_homo, lines_homo, params.score_function);
zengroupIds = groups{1};

function [score, horgroup] = vp_score(vp_homo, lines_homo, score_function)

cos_mat = vp_homo'*lines_homo;
theta_mat = abs(asind(cos_mat));
score_mat = score_function(theta_mat);

tmp = mat2cell(score_mat, ones(size(score_mat,1),1), size(score_mat,2));
horgroup = cellfun(@find, tmp, 'uniformoutput', false);
score = cellfun(@sum, tmp);
