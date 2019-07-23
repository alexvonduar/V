#ifndef _HL_SCORE_HPP_
#define _HL_SCORE_HPP_

#include <opencv2/opencv.hpp>

#include "filter_verhor_lines.hpp"
#include "default_params.hpp"
#include "vp_predict.hpp"

typedef struct
{
    cv::Vec3d horizon_homo;
    double sc;
    std::vector<cv::Vec3d> hvp_homo;
    std::vector<std::vector<int>> hvp_groups;
} Candidate;

/// function [hl_homo, ortho_ modes_left, modes_right, results] = hl_score(hl_samp, infiniteVPs_homo, seglines, xres, yres, focal, zeniths_homo, zengroups, params)
//function [hl_homo, results] = hl_score(hl_samp, ls_homo, z_homo, params)
static inline void hl_score(
    const std::vector<cv::Vec3d> &hl_samp,
    const std::vector<cv::Vec3d> &ls_homo,
    const cv::Vec3d &z_homo,
    const Params &params,
    Candidate &result)
{
    //candidates = repmat(struct(), size(hl_samp,2),1);
    std::vector<Candidate> candidates;
    candidates.reserve(hl_samp.size());
    //nhvps = [];
    std::vector<int> nhvps;
    nhvps.reserve(hl_samp.size());
    //for i = 1:size(hl_samp,2)
    int count = 0;
    for (const auto &samp : hl_samp)
    {
        //helpfulIds = filter_verhor_lines(ls_homo, z_homo, params);
        std::vector<int> helpfulIds;
        filter_verhor_lines(ls_homo, z_homo, params, helpfulIds);
        //--initialIds = 1:numel(helpfulIds);
        if (params.debug_fileid != nullptr)
        {
            fprintf(params.debug_fileid, "%d helpful ids: ", count);
            for (const auto &j : helpfulIds)
            {
                fprintf(params.debug_fileid, "%d ", j);
            }
            fprintf(params.debug_fileid, "\n");
        }
        Candidate candidate;
        //candidates(i).horizon_homo = hl_samp(:,i);
        candidate.horizon_homo = samp;
        //[candidates(i).sc, candidates(i).hvp_homo, hvp_groups] = vp_predict(ls_homo(:,helpfulIds), initialIds, candidates(i).horizon_homo, params);
        std::vector<cv::Vec3d> helpfulLines;
        helpfulLines.reserve(helpfulIds.size());
        std::vector<int> initialIds;
        initialIds.reserve(helpfulIds.size());
        for (int i = 0; i < helpfulIds.size(); ++i)
        {
            auto id = helpfulIds[i];
            helpfulLines.emplace_back(ls_homo[id]);
            initialIds.emplace_back(i);
        }

        std::vector<std::vector<int>> hvp_groups;
        candidate.sc = vp_predict(helpfulLines, initialIds, candidate.horizon_homo, params, candidate.hvp_homo, hvp_groups);
        //candidates(i).hvp_groups = cellfun(@(x) helpfulIds(x), hvp_groups, 'UniformOutput', false);
        for (int i = 0; i < hvp_groups.size(); ++i)
        {
            for (int j = 0; j < hvp_groups[i].size(); ++j)
            {
                auto id = hvp_groups[i][j];
                hvp_groups[i][j] = helpfulIds[id];
            }
        }
        candidate.hvp_groups = hvp_groups;
        //nhvps(i) =  size(candidates(i).hvp_homo,2);
        nhvps.emplace_back(candidate.hvp_homo.size());
        candidates.emplace_back(candidate);
        //end
        ++count;
    }

    // decide the horizon line

    ///horCandidateScores = [candidates.sc];
    ///[~,maxHorCandidateId] = max(horCandidateScores);
    ///hl_homo = candidates(maxHorCandidateId).horizon_homo;
    int best = 0;
    double sc = 0;
    for (int i = 0; i < candidates.size(); ++i)
    {
        if (candidates[i].sc > sc)
        {
            sc = candidates[i].sc;
            best = i;
        }
    }

    // output results

    //results = struct();
    //results.hvp_groups = candidates(maxHorCandidateId).hvp_groups;
    //results.hvp_homo = candidates(maxHorCandidateId).hvp_homo;
    //results.score =  candidates(maxHorCandidateId).sc;
    result = candidates[best];
    //end
}

#endif //_HL_SCORE_HPP_
