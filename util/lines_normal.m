function vp_homo = lines_normal(line_homo, params, b)
% this function solve the homogeneous least squares prob
% A = line_homo', find x s.t. min(||Ax||)
% if b is given, then add a constraint: b'x = 0, which means x is forced to
% be orthogonal to b. [Zhai et al. 2016]

if ~exist('b', 'var')
  
  A = line_homo*line_homo';

  if params.debug_fileid > 0
    fprintf(params.debug_fileid, "lines SVD: \n");
    fprintf(params.debug_fileid, "l: \n");
    for i = 1:size(line_homo,1)
      for j = 1:size(line_homo,2)
        fprintf(params.debug_fileid, "%.1079g ", line_homo(i, j));
      end
      fprintf(params.debug_fileid, "\n");
    end
    fprintf(params.debug_fileid, "lt: \n");
    t = line_homo';
    for i = 1:size(t,1)
      for j = 1:size(t,2)
        fprintf(params.debug_fileid, "%.1079g ", t(i, j));
      end
      fprintf(params.debug_fileid, "\n");
    end
    fprintf(params.debug_fileid, "A: \n");
    for i = 1:size(A,1)
      for j = 1:size(A,2)
        fprintf(params.debug_fileid, "%.1079g ", A(i, j));
      end
      fprintf(params.debug_fileid, "\n");
    end
  end

  [U, ~, ~] = svd(line_homo*line_homo');
  vp_homo = U(:,3);
  
  if params.debug_fileid > 0
    fprintf(params.debug_fileid, "U: \n");
    for i = 1:3
      for j = 1:3
        fprintf(params.debug_fileid, "%.1079g ", U(i, j));
      end
      fprintf(params.debug_fileid, "\n");
    end
  end
  
else
  
  [U, ~, ~] = svd(b);
  p = U(:,2); q = U(:,3);
  Am = line_homo'*[p,q];
  [U, ~, ~] = svd(Am'*Am);
  lambda = U(:,end);
  vp_homo = [p,q]*lambda;
  vp_homo = vp_homo / norm(vp_homo);
  
end

% force to the z-positive semisphere
vp_homo = vp_homo * sign(vp_homo(3)+eps);

% [0 0 0] -> [0 1 0]
if sum(vp_homo) == 0
  vp_homo = [0 1 0]';
end