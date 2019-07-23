% function [hl_homo, ortho_ modes_left, modes_right, results] = hl_score(hl_samp, infiniteVPs_homo, seglines, xres, yres, focal, zeniths_homo, zengroups, params)
function [hl_homo, results] = hl_score(hl_samp, ls_homo, z_homo, params)
candidates = repmat(struct(), size(hl_samp,2),1);
nhvps = [];
for i = 1:size(hl_samp,2)
    helpfulIds = filter_verhor_lines(ls_homo, z_homo, params);
    if params.debug_fileid > 0
        fprintf(params.debug_fileid, "%d helpful ids: ", i - 1);
        for j = 1:length(helpfulIds)
            fprintf(params.debug_fileid, "%d ", helpfulIds(j) - 1);
        end
        fprintf(params.debug_fileid, "\n");
    end
    initialIds = 1:numel(helpfulIds);
    candidates(i).horizon_homo = hl_samp(:,i);
    [candidates(i).sc, candidates(i).hvp_homo, hvp_groups] = vp_predict(ls_homo(:,helpfulIds), initialIds, candidates(i).horizon_homo, params);
    candidates(i).hvp_groups = cellfun(@(x) helpfulIds(x), hvp_groups, 'UniformOutput', false);
    nhvps(i) =  size(candidates(i).hvp_homo,2);
end

%% decide the horizon line

horCandidateScores = [candidates.sc];
[~,maxHorCandidateId] = max(horCandidateScores);
hl_homo = candidates(maxHorCandidateId).horizon_homo;

%% output results

results = struct();
results.hvp_groups = candidates(maxHorCandidateId).hvp_groups;
results.hvp_homo = candidates(maxHorCandidateId).hvp_homo;
results.score =  candidates(maxHorCandidateId).sc;
end

