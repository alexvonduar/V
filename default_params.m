function params = default_params()
params = struct();

% parameters below are provided in the same order and with the 
% same values as in our [Simon et al., ECCV'2018] paper
params.theta_v = pi/8; 
params.theta_z = 10;
params.L_z = 45;
params.theta_h = 1.5;
params.L_h = 64;
params.sigma = 0.2;
params.S = 300;
params.L_vp = 128;

% put the parameter below to 1 if you need to detect infinite 
% horizontal VPs, e.g. for ortho-rectification of fronto-parallel
% vertical planes
params.include_infinite_hvps = 0; 

% parameters below are used in the code provided by [Zhai et al.,
% CVPR'2016], available here : http://cs.uky.edu/~ted/research/fasthorizon/
% these parameters are used to refine the estimates of the vanishing points
% (zenith and horizontal VPs) as mentioned in our paper and detailed in
% [Zhai et al., CVPR'2016]
% Reference:
% Zhai, M., Workman, S., & Jacobs, N. (2016). Detecting Vanishing Points 
% using Global Image Context in a Non-Manhattan World. In IEEE Conference 
% on Computer Vision and Pattern Recognition (CVPR).
params.theta_con = 1.5; 
params.score_function = @(x) (params.theta_con-x) .* (params.theta_con > x);
params.theta_verline = 15; 
params.theta_horline = 2; 
params.hvp_refinement = true; 
params.refine_niters = 3; 