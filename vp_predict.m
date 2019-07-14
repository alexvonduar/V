function [sc, horvps_homo, horgroups] = vp_predict(lines_homo, initialIds, horizon_homo, params)
%% compute intersections between the line segments and the HL candidate

inter_homo = bsxfun(@cross, lines_homo(:,initialIds), horizon_homo);
inter_homo = bsxfun(@rdivide, inter_homo, sqrt(sum(inter_homo.^2,1)));
inter_homo = bsxfun(@times, inter_homo, sign(inter_homo(3,:)+ eps));
inter_pts = bsxfun(@rdivide, inter_homo(1:2,:), inter_homo(3,:));

%% compute the MMMs of the coordinate histogtam

max_modes = [];
p = [];
a = horizon_homo(1);
b = horizon_homo(2);
c = horizon_homo(3);
A_hmg = cross(horizon_homo,[b -a 0]);
A(1) = A_hmg(1)/A_hmg(3);
A(2) = A_hmg(2)/A_hmg(3);
rho = abs(c)/sqrt(a^2+b^2);
rho2 = sqrt(inter_pts(1,:).*inter_pts(1,:)+inter_pts(2,:).*inter_pts(2,:));
if rho > 1
    p = acos(rho./rho2)/pi;
else
    d = sqrt(abs(rho2.*rho2-rho^2));
    I = find(rho2 <= 1);
    if ~isempty(I)
        p(I) = d(I)/pi;
    end
    I = find(rho2 > 1);
    if ~isempty(I)
        d2 = sqrt(rho2.*rho2-1);
        beta = atan(d2);
        p(I) = (beta(I)+d(I)-d2(I))/pi;
    end
end
dt = [b -a]*(inter_pts-A'*ones(1,size(inter_pts,2)));
I = find(dt < 0);
p(I) = -p(I);
[N,edges] = histcounts(p,params.L_vp);
[max_modes, H] = mnf_modes(N,400);
if isempty(max_modes)
    max_modes = [];
    H = 0;
else
    [~,I] = sort(H, 'descend');
    H = H(I);
    max_modes = max_modes(I,:);
end
horgroups = cell(0,0);
scores = [];
horvps_homo = [];
for i = 1:size(max_modes,1)
    Ni = zeros(1,size(N,2));
    a = max_modes(i,1);
    b = max_modes(i,2);
    Ni(a:b) = N(a:b);
    [m,j] = max(Ni);
    p_i = (edges(j)+edges(j+1))/2;
    [~,vpId] = min(abs(p-p_i));
    
    horvps_homo(1:3,end+1) = inter_homo(:,vpId);
    scores(end+1) = m;
    edgesId = intersect(find(p >= edges(a)),find(p <= edges(b+1)));
    horgroups{end+1,1} =  edgesId;
end

if isempty(max_modes)
    scores = [0];
    sc = 0;
else
    %refine the VPs according to [Zhai et al.] and/or compute the scores
    if params.hvp_refinement
        horvps_homo = vp_refinement(lines_homo, horvps_homo, horizon_homo, params);
    end
    [scores, horgroups] = vp_score(horvps_homo, lines_homo, params.score_function);
    
    % sorted by score
    [scores, sortIds] = sort(scores, 'descend');
    horvps_homo = horvps_homo(:,sortIds);
    horgroups = horgroups(sortIds);
    nvps = min(numel(horgroups), 2);
    sc = sum(scores(1:nvps));
end
end

function [score, horgroup] = vp_score(vp_homo, lines_homo, score_function) % from [Zhai et al.]
cos_mat = vp_homo'*lines_homo;
theta_mat = abs(asind(cos_mat));
score_mat = score_function(theta_mat);
tmp = mat2cell(score_mat, ones(size(score_mat,1),1), size(score_mat,2));
horgroup = cellfun(@find, tmp, 'uniformoutput', false);
score = cellfun(@sum, tmp);
end
