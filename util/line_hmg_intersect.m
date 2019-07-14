function p = line_hmg_intersect(l1_hmg,l2_hmg)
p_hmg = cross(l1_hmg,l2_hmg);
p(1,1) = p_hmg(1)/p_hmg(3);
p(1,2) = p_hmg(2)/p_hmg(3);
end