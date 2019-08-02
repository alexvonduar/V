%%-------------------------------------------------------------------------
%
% >V<  
% 
% Copyright 2018 Gilles Simon <gsimon@loria.fr> Université de Lorraine
%
% Version 1.0 - September 5, 2018
%
% >V< is an implementation of the vanishing point detector described in the
% paper:
%
%   "A-Contrario Horizon-First Vanishing Point Detection Using Second-Order
%   Grouping Laws" by Gilles Simon, Antoine Fond, and Marie-Odile Berger,  
%   European Conference on Computer Vision, Sep 2018, Munich, Germany." 
%
% available on our website, along with a few example results:
% https://members.loria.fr/GSimon/software/v/
%
% If you use the code, please cite the paper. 
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%-------------------------------------------------------------------------

function [hl, hvps, hvp_groups, z, z_group, ls] = V(name, im, width, height, focal, params)
% path to LSD
LSD_BIN_ = 'util/lsd/lsd';

% plot options (0 if not plotted, figure number otherwise)
plots = struct();
plots.lsd = 0; %10
plots.ls_filter = 0; %11
plots.zl = 0; %12
plots.hl_modes = 0; %12
plots.hl_samples = 0; %12

% fix random seed
rng(1)

% principal point is assumed at image center
u0 = width/2;
v0 = height/2;

%% line segment (LS) extraction based on the LSD algorithm
%
% Reference: 
% Rafael Grompone von Gioi, Jérémie Jakubowicz, Jean-Michel Morel, and Gregory Randall, 
% LSD: a Line Segment Detector, Image Processing On Line, 2 (2012), pp. 35–55. 
% https://doi.org/10.5201/ipol.2012.gjmr-lsd
[lsd, ~] =  extract_linesegment(name, im, LSD_BIN_);

% plot the results
if plots.lsd
    figure(plots.lsd); clf(); imshow(im); hold on;
    cmap = colormap(hsv(size(lsd,1)));
    for j = 1:size(lsd,1)
        plot([lsd(j,1);lsd(j,3)],[lsd(j,2);lsd(j,4)], 'Color',[cmap(j,1) cmap(j,2) cmap(j,3)],'LineWidth',2);
    end
    drawnow;
end

if 0 && params.debug_fileid > 0
    [~, I] = sort(lsd(:,1), 'ascend');
    lsd = lsd(I,:);
    fprintf(params.debug_fileid, "%d lsd results ----\n", size(lsd, 1));
    for j = 1:size(lsd, 1)
        fprintf(params.debug_fileid, "[%.149g, %.149g], [%.149g, %.149g]\n", lsd(j,1), lsd(j,2), lsd(j,3), lsd(j,4));
    end
end

%% LS filtering

thres_aligned =  max(width,height)/128.;
length_t = sqrt(width+height)/1.71;
ls = ls_filter(thres_aligned, length_t, lsd);
ls_homo = normalize(ls, width, height, focal);

% plot the results
if plots.ls_filter
    figure(plots.ls_filter); clf(); imshow(im); hold on;
    cmap = colormap(hsv(size(ls,1)));
    for j = 1:size(ls,1)
        plot([ls(j,1);ls(j,3)],[ls(j,2);ls(j,4)], 'Color',[cmap(j,1) cmap(j,2) cmap(j,3)],'LineWidth',2);
    end
    drawnow;
end

if 0 && params.debug_fileid > 0
    fprintf(params.debug_fileid, "%d filtered lsd results ----\n", size(ls, 1));
    for j = 1:size(ls, 1)
        fprintf(params.debug_fileid, "[%.1074g, %.1074g], [%.1074g, %.1074g] homo [%.1074g, %.1074g, %.1074g]\n", ls(j,1), ls(j,2), ls(j,3), ls(j,4), ls_homo(1, j), ls_homo(2, j), ls_homo(3, j));
    end
end

%% ZL and zenith rough predictions

% prediction of the zenith line
dist_max = width/8;
zl = zl_predict(lsd, dist_max, u0, v0, width, height, params);
zl_homo = cell(0,0);
z_homo_cand = cell(0);
z_group_cand = cell(0);
for i = 1:length(zl)
    zl_homo{i} = normalize([zl(i) 0 u0 v0], width, height, focal);
    [z_homo_cand{i}, z_group_cand{i}] = z_predict(ls_homo, zl_homo{i}, params, 0); 
end

% plot the results
if plots.zl
    figure(plots.zl); clf(); imshow(im); hold on;
    cmap = colormap(hsv(length(z_homo_cand)));
    for i = 1:length(z_homo_cand)
        z = unnormalize(z_homo_cand{i}, width, height, focal, 0);
        plot([width/2;z(1)],[height/2;z(2)], 'Color',[cmap(i,1) cmap(i,2) cmap(i,3)],'LineWidth',3);
    end
    drawnow();
end

if params.debug_fileid > 0
    num_z = length(zl);
    fprintf(params.debug_fileid, "num zenith line pred: %d\n", num_z);
    for i = 1:num_z
        fprintf(params.debug_fileid, "zenith line prediction: %d ----\n", i - 1);
        fprintf(params.debug_fileid, "zl: %f zl_homo: [%.1074g, %.1074g, %.1074g]\n", zl(i), zl_homo{i}(1), zl_homo{i}(2), zl_homo{i}(3));
        num_z_homo_cand = size(z_homo_cand{i}, 2);
        fprintf(params.debug_fileid, "num cand %d\n", num_z_homo_cand);
        for j = 1:num_z_homo_cand
            fprintf(params.debug_fileid, "cand %d: [%.13g %.13g %.13g]\n", j - 1, z_homo_cand{i}(1, j), z_homo_cand{i}(2, j), z_homo_cand{i}(3, j));
            num_group_cand = size(z_group_cand{i}, 2);
            fprintf(params.debug_fileid, "group cand %d: ", num_group_cand);
            for k = 1:num_group_cand
                fprintf(params.debug_fileid, "%d ", z_group_cand{i}(j, k) - 1);
            end
            fprintf(params.debug_fileid, "\n");
        end
    end
end

%% choose the best zenith candidate based on the relevance of the predicted HLs

best_z_cand = 1;
best_z_score = 0;
for i = 1:length(zl_homo)
    
    % HL prediction
    [modes_homo, ~, ~, ~, ~] = hl_predict(lsd, z_homo_cand{i}, u0, v0, width, height, focal, params);
    
    % HL scoring (for performance optimization, each zenith candidate is
    % assessed based only on the meaningful HLs (no sampling is performed at
    % that step))
    
    [~, results] = hl_score(modes_homo, ls_homo, z_homo_cand{i}, params);
 
    if params.debug_fileid > 0
      fprintf(params.debug_fileid, "%dth predicted hls -- \n", i);
      n_modes = size(modes_homo, 2);
      for j = 1:n_modes
        fprintf(params.debug_fileid, "modes %d: [%f %f %f]\n", j - 1, modes_homo(1, j), modes_homo(2, j), modes_homo(3, j));
      end
      fprintf(params.debug_fileid, "results: score %f\n", results.score);
      n_hvp = size(results.hvp_homo, 2);
      for j = 1:n_hvp
          fprintf(params.debug_fileid, "hvp %d: [%f %f %f]\n", j - 1, results.hvp_homo(1,j), results.hvp_homo(2,j), results.hvp_homo(3, j));
      end
      n_groups = size(results.hvp_groups, 1);
      for j = 1:n_groups
          fprintf(params.debug_fileid, "group %d: ", j - 1);
          group = results.hvp_groups{j};
          n = size(group, 1);
          for k = 1:n
              fprintf(params.debug_fileid, "%d ", group(k, 1) - 1);
          end
          fprintf(params.debug_fileid, "\n");
      end
    end

    % keep the zenith candidate with highest score
    if results.score > best_z_score
        best_z_cand = i;
        best_z_score = results.score;
    end
end


if params.debug_fileid > 0
    fprintf(params.debug_fileid, "best z: %d score: %f\n", best_z_cand - 1, best_z_score);

    %fprintf(params.debug_fileid, "%d filtered lsd results2 ----\n", size(ls, 1));
    %for j = 1:size(ls, 1)
    %    fprintf(params.debug_fileid, "[%.1074g, %.1074g], [%.1074g, %.1074g] homo [%.1074g, %.1074g, %.1074g]\n", ls(j,1), ls(j,2), ls(j,3), ls(j,4), ls_homo(1, j), ls_homo(2, j), ls_homo(3, j));
    %end
end

% zenith refinement (based on Zhang et al. method)
[z_homo_cand{best_z_cand}, z_group_cand{best_z_cand}] = z_predict(ls_homo, zl_homo{best_z_cand}, params, 1);

if params.debug_fileid > 0
    homo = z_homo_cand{best_z_cand};
    fprintf(params.debug_fileid, "refiend best z: [%f %f %f]\n", homo(1, 1), homo(2, 1), homo(3, 1));
    fprintf(params.debug_fileid, "refined groups: ");
    group = z_group_cand{best_z_cand};
    for i = 1:size(group, 2)
        fprintf(params.debug_fileid, "%d ", group(1, i) - 1);
    end
    fprintf(params.debug_fileid, "\n");
end

% HL prediction
[modes_homo, modes_offset, modes_left, modes_right, H] = hl_predict(lsd, z_homo_cand{best_z_cand}, u0, v0, width, height, focal, params);

if params.debug_fileid > 0
    fprintf(params.debug_fileid, "hl prediction --\n");
    for i = 1:size(modes_homo, 2)
        fprintf(params.debug_fileid, "%d: [%f %f %f] offset %f left %f right %f H %f\n", i - 1, modes_homo(1, i), modes_homo(2, i), modes_homo(3, i), modes_offset(i), modes_left(i), modes_right(i), H(i));
    end
end

% HL sampling
[samp_homo, samp_left, samp_right] = hl_sample(z_homo_cand{best_z_cand}, modes_homo, modes_offset, modes_left, modes_right, H, u0, v0, width, height, focal, params);

if 0 && params.debug_fileid > 0
    fprintf(params.debug_fileid, "hl sampling --\n");
    for i = 1:size(samp_homo, 2)
        fprintf(params.debug_fileid, "%d: [%f %f %f] left %f right %f\n", i - 1, samp_homo(1, i), samp_homo(2, i), samp_homo(3, i), samp_left(i), samp_right(i));
    end
end

% plot the results
if plots.hl_samples
    cmap = colormap(hsv(length(zl_homo)));
    figure(plots.hl_samples);
    if plots.zl ~= plots.hl_samples
        clf(); imshow(im); hold on;
    end
    for j = 1:size(samp_homo,2)
        plot([0;width],[samp_left(j);samp_right(j)], 'Color',[cmap(i,1) cmap(i,2) cmap(i,3)],'LineWidth',1);
    end
end
if plots.hl_modes
    figure(plots.hl_modes);
    if plots.zl ~= plots.hl_modes && plots.hl_samples ~= plots.hl_modes
        clf(); imshow(im); hold on;
    end
    for j = 1:size(modes_homo,2)
        if H(j) > 0
            plot([0;width],[modes_left(j);modes_right(j)],'b--','LineWidth',2);
        end
    end
end

% HL scoring
[hl_homo, results] = hl_score(samp_homo, ls_homo, z_homo_cand{best_z_cand}, params);
hl = unnormalize(hl_homo, width, height, focal, 1);
hvps = unnormalize(results.hvp_homo, width, height, focal, 0);
hvp_groups = results.hvp_groups;
z = unnormalize(z_homo_cand{best_z_cand}, width, height, focal, 0);
z_group = z_group_cand{best_z_cand};

if params.debug_fileid > 0
    fprintf(params.debug_fileid, "z: [%f %f] [%f %f %f]\n", z(1), z(2), z_homo_cand{best_z_cand}(1), z_homo_cand{best_z_cand}(2), z_homo_cand{best_z_cand}(3));
    fprintf(params.debug_fileid, "z group: ");
    for i = 1:length(z_group)
        fprintf(params.debug_fileid, "%d ", z_group(i) - 1);
    end
    fprintf(params.debug_fileid, "\n");
    hvp_count = size(hvps,1);
    for v = 1:hvp_count
        fprintf(params.debug_fileid, "%d vp: [%f %f]\n", v - 1, hvps(v, 1), hvps(v, 2));
    end
    for j = 1:numel(hvp_groups)
        hg = hvp_groups{j};
        fprintf(params.debug_fileid, "%d vp group: ", j - 1);
        for k = 1:length(hg)
            fprintf(params.debug_fileid, "%d ", hg(k) - 1);
        end
        fprintf(params.debug_fileid, "\n");
    end
end
end
