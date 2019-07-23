#ifndef _VP_PREDICT_HPP_
#define _VP_PREDICT_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"
#include "util/histcounts.hpp"
#include "util/mnf_modes.hpp"
#include "util/vp_refinement.hpp"

#if 0
/*
function [score, horgroup] = vp_score(vp_homo, lines_homo, score_function) % from [Zhai et al.]
cos_mat = vp_homo'*lines_homo;
theta_mat = abs(asind(cos_mat));
score_mat = score_function(theta_mat);
tmp = mat2cell(score_mat, ones(size(score_mat,1),1), size(score_mat,2));
horgroup = cellfun(@find, tmp, 'uniformoutput', false);
score = cellfun(@sum, tmp);
end
*/
#else
#include "util/vp_score.hpp"
#endif

//function [sc, horvps_homo, horgroups] = vp_predict(lines_homo, initialIds, horizon_homo, params)
static inline double vp_predict(const std::vector<cv::Vec3d> &lines_homo,
                                const std::vector<int> &initialIds,
                                const cv::Vec3d &horizon_homo,
                                const Params &params,
                                std::vector<cv::Vec3d> &horvps_homo,
                                std::vector<std::vector<int>> &horgroups

)
{
    // compute intersections between the line segments and the HL candidate
    std::vector<cv::Vec3d> inter_homo;
    inter_homo.reserve(initialIds.size());
    std::vector<cv::Vec2d> inter_pts;
    inter_pts.reserve(initialIds.size());
    for (const auto &id : initialIds)
    {
        //inter_homo = bsxfun(@cross, lines_homo(:,initialIds), horizon_homo);
        auto inter = lines_homo[id].cross(horizon_homo);
        //inter_homo = bsxfun(@rdivide, inter_homo, sqrt(sum(inter_homo.^2,1)));
        inter /= cv::norm(inter);
        //inter_homo = bsxfun(@times, inter_homo, sign(inter_homo(3,:)+ eps));
        if (inter[2] < 0)
        {
            inter *= -1;
        }
        //inter_pts = bsxfun(@rdivide, inter_homo(1:2,:), inter_homo(3,:));
        inter_pts.emplace_back(cv::Vec2d{inter[0] / inter[2], inter[1] / inter[2]});
        inter_homo.emplace_back(inter);
    }

    // compute the MMMs of the coordinate histogtam

    //max_modes = [];
    //p = [];
    std::vector<double> p;
    p.reserve(inter_pts.size());
    //a = horizon_homo(1);
    auto a = horizon_homo[0];
    //b = horizon_homo(2);
    auto b = horizon_homo[1];
    //c = horizon_homo(3);
    auto c = horizon_homo[2];
    //A_hmg = cross(horizon_homo,[b -a 0]);
    auto A_hmg = horizon_homo.cross(cv::Vec3d{b, -a, 0});
    //A(1) = A_hmg(1)/A_hmg(3);
    //A(2) = A_hmg(2)/A_hmg(3);
    cv::Vec2d A{A_hmg[0] / A_hmg[2], A_hmg[1] / A_hmg[2]};
    //rho = abs(c)/sqrt(a^2+b^2);
    auto rho = std::abs(c) / std::sqrt(a * a + b * b);
    //rho2 = sqrt(inter_pts(1,:).*inter_pts(1,:)+inter_pts(2,:).*inter_pts(2,:));
    //if rho > 1
    //    p = acos(rho./rho2)/pi;
    //else
    //    d = sqrt(abs(rho2.*rho2-rho^2));
    //    I = find(rho2 <= 1);
    //    if ~isempty(I)
    //        p(I) = d(I)/pi;
    //    end
    //    I = find(rho2 > 1);
    //    if ~isempty(I)
    //        d2 = sqrt(rho2.*rho2-1);
    //        beta = atan(d2);
    //        p(I) = (beta(I)+d(I)-d2(I))/pi;
    //    end
    //end
    //dt = [b -a]*(inter_pts-A'*ones(1,size(inter_pts,2)));
    //I = find(dt < 0);
    //p(I) = -p(I);
    cv::Vec2d ba{b, -a};
    for (const auto &pts : inter_pts)
    {
        auto rho2 = cv::norm(pts);
        auto dt = ba.dot(pts - A);
        double _p;
        if (rho > 1)
        {
            _p = std::acos(rho / rho2) / CV_PI;
        }
        else
        {
            auto d = std::sqrt(std::abs(rho2 * rho2 - rho * rho));
            if (rho2 <= 1)
            {
                _p = d / CV_PI;
            }
            else
            {
                auto d2 = std::sqrt(rho2 * rho2 - 1);
                auto beta = std::atan(d2);
                _p = (beta + d - d2) / CV_PI;
            }
        }
        if (dt < 0)
        {
            p.emplace_back(-_p);
        }
        else
        {
            p.emplace_back(_p);
        }
    }

    //[N,edges] = histcounts(p,params.L_vp);
    cv::Mat N;
    std::vector<double> edges;
    histcounts(p, params.L_vp, nullptr, N, edges);

    //[max_modes, H] = mnf_modes(N,400);
    std::vector<std::pair<int, int>> max_modes;
    std::vector<double> H;
    mnf_modes(N, 400, max_modes, H);
    //if isempty(max_modes)
    //    max_modes = [];
    //    H = 0;
    //else
    //    [~,I] = sort(H, 'descend');
    //    H = H(I);
    //    max_modes = max_modes(I,:);
    //end
    if (max_modes.size())
    {
        mnf_modes_sort(max_modes, H);
    }
    //horgroups = cell(0,0);
    ///std::vector<std::vector<int>> horgroups;
    horgroups.reserve(max_modes.size());
    //scores = [];
    std::vector<double> scores;
    scores.reserve(max_modes.size());
    //horvps_homo = [];
    ///std::vector<cv::Vec3d> horvps_homo;
    horvps_homo.reserve(max_modes.size());
    //for i = 1:size(max_modes,1)
    for (const auto mode : max_modes)
    {
        //Ni = zeros(1,size(N,2));
        //a = max_modes(i,1);
        //b = max_modes(i,2);
        //Ni(a:b) = N(a:b);
        auto m = N.at<double>(0, mode.first);
        auto j = mode.first;
        //[m,j] = max(Ni);
        for (int i = j + 1; i <= mode.second; ++i)
        {
            auto _m = N.at<double>(0, i);
            if (_m > m)
            {
                m = _m;
                j = i;
            }
        }
        //p_i = (edges(j)+edges(j+1))/2;
        auto p_i = (edges[j] + edges[j + 1]) / 2;
        //[~,vpId] = min(abs(p-p_i));
        auto vpId = 0;
        auto delta = std::abs(p[0] - p_i);
        for (int i = 1; i < p.size(); ++i)
        {
            auto _delta = std::abs(p[i] - p_i);
            if (_delta < delta)
            {
                delta = _delta;
                vpId = i;
            }
        }

        //horvps_homo(1:3,end+1) = inter_homo(:,vpId);
        horvps_homo.emplace_back(inter_homo[vpId]);
        //scores(end+1) = m;
        scores.emplace_back(m);
        //edgesId = intersect(find(p >= edges(a)),find(p <= edges(b+1)));
        std::vector<int> group;
        group.reserve(p.size());
        for (int i = 0; i < p.size(); ++i)
        {
            if (p[i] >= edges[mode.first] and p[i] <= edges[mode.second + 1])
            {
                group.emplace_back(i);
            }
        }
        //horgroups{end+1,1} =  edgesId;
        horgroups.emplace_back(group);
        //end
    }

    double sc = 0;
    scores.clear();
    horgroups.clear();
    //if isempty(max_modes)
    if (max_modes.size() == 0)
    {
        //scores = [0];
        //sc = 0;
        //-scores.clear();
    }
    else
    {
        //refine the VPs according to [Zhai et al.] and/or compute the scores
        if (params.hvp_refinement)
        {
            //horvps_homo = vp_refinement(lines_homo, horvps_homo, horizon_homo, params);
            //end
            vp_refinement(lines_homo, horvps_homo, horizon_homo, params);
        }
        //[scores, horgroups] = vp_score(horvps_homo, lines_homo, params.score_function);
        for (const auto &horvps : horvps_homo)
        {
            std::vector<int> horgroup;
            auto score = vp_score(horvps, lines_homo, params, horgroup);
            scores.emplace_back(score);
            horgroups.emplace_back(horgroup);
        }

        // sorted by score
        //[scores, sortIds] = sort(scores, 'descend');
        std::vector<std::pair<double, int>> index;
        index.reserve(scores.size());
        for (int i = 0; i < scores.size(); ++i)
        {
            index.emplace_back(std::pair<double, int>{scores[i], i});
        }
        std::sort(index.begin(), index.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
        //horvps_homo = horvps_homo(:,sortIds);
        //horgroups = horgroups(sortIds);
        std::vector<cv::Vec3d> horvps_homo_sorted;
        horvps_homo_sorted.reserve(horvps_homo.size());
        std::vector<std::vector<int>> horgroups_sorted;
        horgroups_sorted.reserve(horgroups.size());
        for (int i = 0; i < index.size(); ++i)
        {
            scores[i] = index[i].first;
            auto id = index[i].second;
            horvps_homo_sorted.emplace_back(horvps_homo[id]);
            horgroups_sorted.emplace_back(horgroups[id]);
        }
        horvps_homo = horvps_homo_sorted;
        horgroups = horgroups_sorted;
        //nvps = min(numel(horgroups), 2);
        auto nvps = horgroups.size() < 2 ? horgroups.size() : 2;
        //sc = sum(scores(1:nvps));
        for (int i = 0; i < nvps; ++i)
        {
            sc += scores[i];
        }
        //end
    }
    return sc;
    //end
}

#endif //_VP_PREDICT_HPP_
