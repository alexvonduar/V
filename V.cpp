#include <string>

#include <opencv2/opencv.hpp>

#include "ls_filter.hpp"
#include "util/normalize.hpp"
#include "default_params.hpp"
#include "zl_predict.hpp"
#include "z_predict.hpp"
#include "hl_predict.hpp"
#include "hl_score.hpp"
#include "hl_sample.hpp"
#include "util/unnormalize.hpp"
#include "lsd.hpp"
#include "draw.hpp"
#include "calibrate.hpp"
#include "orthorectify.hpp"

int V(const Params &params,
      const cv::Mat &img,
      std::vector<cv::Vec2d> &hl,
      std::vector<cv::Vec2d> &hvps,
      std::vector<std::vector<int>> &hvp_groups,
      cv::Vec2d &z,
      std::vector<int> &z_group,
      std::vector<cv::Vec4d> &ls)
{
    //Params params;
    //default_params(params);

    //auto img = cv::imread(fname);
    //if (img.empty())
    //{
    //    return -1;
    //}

    std::vector<cv::Vec4d> lsd;
    // detect line segments
    lsd_detect(img, lsd);

    if (1)
    {
        draw("lsd", img, lsd, cv::Scalar(0, 0, 255), true);
    }

    if (0 and params.debug_fileid != nullptr)
    {
        std::sort(lsd.begin(), lsd.end(), [](const auto &a, const auto &b) { return a[0] < b[0]; });
        auto num_lines = lsd.size();
        fprintf(params.debug_fileid, "%ld lsd results ----\n", num_lines);
        for (int i = 0; i < num_lines; ++i)
        {
            fprintf(params.debug_fileid, "[%.1074g, %.1074g], [%.1074g, %.1074g]\n", lsd[i][0], lsd[i][1], lsd[i][2], lsd[i][3]);
        }
    }

    // filter line segements
    auto width = img.cols;
    auto height = img.rows;

    double focal = width > height ? double(width) / 2. : double(height) / 2.; //fake focal;

    // principal point is assumed at image center
    auto u0 = width / 2;
    auto v0 = height / 2;

    auto thres_aligned = double(std::max(width, height)) / 128.;
    auto length_t = double(std::sqrt(width + height)) / 1.71;
    //std::vector<cv::Vec4d> ls;
    ls_filter(thres_aligned, length_t, lsd, ls);
    std::vector<cv::Vec3d> ls_homo;
    ls_homo.reserve(ls.size());
    for (const auto &l : ls)
    {
        cv::Vec3d l_homo;
        normalize(l, width, height, focal, l_homo);
        ls_homo.emplace_back(l_homo);
    }

    if (1)
    {
        draw("fileted lsd", img, ls, cv::Scalar(0, 0, 255), true);
    }

    if (0 and params.debug_fileid != nullptr)
    {
        auto nfilted = ls.size();
        fprintf(params.debug_fileid, "%ld filtered lsd results ----\n", nfilted);
        for (int i = 0; i < nfilted; ++i)
        {
            fprintf(params.debug_fileid, "[%.1074g, %.1074g], [%.1074g, %.1074g] homo [%.1074g, %.1074g, %.1074g]\n", ls[i][0], ls[i][1], ls[i][2], ls[i][3], ls_homo[i][0], ls_homo[i][1], ls_homo[i][2]);
        }
    }

    // ZL and zenith rough predictions

    // prediction of the zenith line
    auto dist_max = double(width) / 8.;
    // zl = zl_predict(lst, dist_max, u0, v0, width, height, params);
    std::vector<double> zl;
    zl_predict(lsd, dist_max, u0, v0, width, height, params, zl);

    std::vector<cv::Vec3d> zl_homo;
    std::vector<cv::Vec3d> z_homo_cand;
    std::vector<std::vector<int>> z_group_cand;
    for (const auto &_zl : zl)
    {
        cv::Vec3d z_homo;
        normalize(cv::Vec4d{_zl, 0, double(u0), double(v0)}, width, height, focal, z_homo);
        cv::Vec3d homo;
        std::vector<int> group;
        z_predict(ls_homo, z_homo, params, false, homo, group);
        zl_homo.emplace_back(z_homo);
        z_homo_cand.emplace_back(homo);
        z_group_cand.emplace_back(group);
    }

    if (params.debug_fileid != nullptr)
    {
        auto num_z = zl.size();
        fprintf(params.debug_fileid, "num zenith line pred: %ld\n", num_z);
        for (int i = 0; i < num_z; ++i)
        {
            fprintf(params.debug_fileid, "zenith line prediction: %d ----\n", i);
            fprintf(params.debug_fileid, "zl: %f zl_homo: [%.1074g, %.1074g, %.1074g]\n", zl[i], zl_homo[i][0], zl_homo[i][1], zl_homo[i][2]);
            auto num_z_homo_cand = z_homo_cand.size();
            fprintf(params.debug_fileid, "num cand %ld\n", num_z_homo_cand);
            for (int j = 0; j < num_z_homo_cand; ++j)
            {
                fprintf(params.debug_fileid, "cand %d: [%.13g %.13g %.13g]\n", j, z_homo_cand[j][0], z_homo_cand[j][1], z_homo_cand[j][2]);
                auto num_group_cand = z_group_cand[i].size();
                fprintf(params.debug_fileid, "group cand %ld: ", num_group_cand);
                for (int k = 0; k < num_group_cand; ++k)
                {
                    fprintf(params.debug_fileid, "%d ", z_group_cand[i][k]);
                }
                fprintf(params.debug_fileid, "\n");
            }
        }
    }

    // choose the best zenith candidate based on the relevance of the predicted HLs

    int best_z_cand = 0;
    double best_z_score = 0;
    //for i = 1:length(zl_homo)
    for (int i = 0; i < zl_homo.size(); ++i)
    {
        // HL prediction
        //[modes_homo, ~, ~, ~, ~] = hl_predict(lsd, z_homo_cand{i}, u0, v0, width, height, focal, params);
        std::vector<cv::Vec3d> modes_homo;
        std::vector<double> modes_offset;
        std::vector<double> modes_left;
        std::vector<double> modes_right;
        std::vector<double> H;
        hl_predict(lsd, z_homo_cand[i], u0, v0, width, height, focal, params, modes_homo, modes_offset, modes_left, modes_right, H);

        // HL scoring (for performance optimization, each zenith candidate is
        // assessed based only on the meaningful HLs (no sampling is performed at
        // that step))

        //[~, results] = hl_score(modes_homo, ls_homo, z_homo_cand{i}, params);
        Candidate result;
        hl_score(modes_homo, ls_homo, z_homo_cand[i], params, result);

        if (params.debug_fileid > 0)
        {
            fprintf(params.debug_fileid, "%dth predicted hls -- \n", i + 1);
            auto n_modes = modes_homo.size();
            for (int j = 0; j < n_modes; ++j)
            {
                fprintf(params.debug_fileid, "modes %d: [%f %f %f]\n", j, modes_homo[j][0], modes_homo[j][1], modes_homo[j][2]);
            }
            fprintf(params.debug_fileid, "results: score %f\n", result.sc);
            for (int j = 0; j < result.hvp_homo.size(); ++j)
            {
                fprintf(params.debug_fileid, "hvp %d: [%f %f %f]\n", j, result.hvp_homo[j][0], result.hvp_homo[j][1], result.hvp_homo[j][2]);
            }
            auto n_groups = result.hvp_groups.size();
            for (int j = 0; j < n_groups; ++j)
            {
                fprintf(params.debug_fileid, "group %d: ", j);
                auto &group = result.hvp_groups[j];
                auto n = group.size();
                for (int k = 0; k < n; ++k)
                {
                    fprintf(params.debug_fileid, "%d ", group[k]);
                }
                fprintf(params.debug_fileid, "\n");
            }
        }

        // keep the zenith candidate with highest score
        //if results.score > best_z_score
        if (result.sc > best_z_score)
        {
            best_z_cand = i;
            best_z_score = result.sc;
            //end
        }
        //end
    }

    if (params.debug_fileid > 0)
    {
        fprintf(params.debug_fileid, "best z: %d score: %f\n", best_z_cand, best_z_score);

        //auto nfilted = ls.size();
        //fprintf(params.debug_fileid, "%ld filtered lsd results2 ----\n", nfilted);
        //for (int i = 0; i < nfilted; ++i)
        //{
        //    fprintf(params.debug_fileid, "[%.1074g, %.1074g], [%.1074g, %.1074g] homo [%.1074g, %.1074g, %.1074g]\n", ls[i][0], ls[i][1], ls[i][2], ls[i][3], ls_homo[i][0], ls_homo[i][1], ls_homo[i][2]);
        //}
    }

    /// zenith refinement (based on Zhang et al. method)
    //[z_homo_cand{best_z_cand}, z_group_cand{best_z_cand}] = z_predict(ls_homo, zl_homo{best_z_cand}, params, 1);
    z_predict(ls_homo, zl_homo[best_z_cand], params, true, z_homo_cand[best_z_cand], z_group_cand[best_z_cand]);

    if (params.debug_fileid > 0)
    {
        fprintf(params.debug_fileid, "refiend best z: [%f %f %f]\n", z_homo_cand[best_z_cand][0], z_homo_cand[best_z_cand][1], z_homo_cand[best_z_cand][2]);
        fprintf(params.debug_fileid, "refined groups: ");
        auto n = z_group_cand[best_z_cand].size();
        for (int i = 0; i < n; ++i)
        {
            fprintf(params.debug_fileid, "%d ", z_group_cand[best_z_cand][i]);
        }
        fprintf(params.debug_fileid, "\n");
    }

    /// HL prediction
    //[modes_homo, modes_offset, modes_left, modes_right, H] = hl_predict(lsd, z_homo_cand{best_z_cand}, u0, v0, width, height, focal, params);
    std::vector<cv::Vec3d> modes_homo;
    std::vector<double> modes_offset;
    std::vector<double> modes_left;
    std::vector<double> modes_right;
    std::vector<double> H;
    hl_predict(lsd, z_homo_cand[best_z_cand], u0, v0, width, height, focal, params, modes_homo, modes_offset, modes_left, modes_right, H);
    if (params.debug_fileid != nullptr)
    {
        fprintf(params.debug_fileid, "hl prediction --\n");
        for (int i = 0; i < modes_homo.size(); ++i)
        {
            fprintf(params.debug_fileid, "%d: [%f %f %f] offset %f left %f right %f H %f\n", i, modes_homo[i][0], modes_homo[i][1], modes_homo[i][2], modes_offset[i], modes_left[i], modes_right[i], H[i]);
        }
    }

    /// HL sampling
    //[samp_homo, samp_left, samp_right] = hl_sample(z_homo_cand{best_z_cand}, modes_homo, modes_offset, modes_left, modes_right, H, u0, v0, width, height, focal, params);
    std::vector<cv::Vec3d> samp_homo;
    std::vector<double> samp_left;
    std::vector<double> samp_right;
    hl_sample(z_homo_cand[best_z_cand], modes_homo, modes_offset, modes_left, modes_right, H, u0, v0, width, height, focal, params, samp_homo, samp_left, samp_right);

    if (params.debug_fileid > 0)
    {
        fprintf(params.debug_fileid, "hl sampling --\n");
        for (int i = 0; i < samp_homo.size(); ++i)
        {
            fprintf(params.debug_fileid, "%d: [%f %f %f] left %f right %f\n", i, samp_homo[i][0], samp_homo[i][1], samp_homo[i][2], samp_left[i], samp_right[i]);
        }
    }

    /// HL scoring
    //[hl_homo, results] = hl_score(samp_homo, ls_homo, z_homo_cand{best_z_cand}, params);
    Candidate r;
    hl_score(samp_homo, ls_homo, z_homo_cand[best_z_cand], params, r);
    //hl = unnormalize(hl_homo, width, height, focal, 1);
    //std::vector<cv::Vec2d> hl;
    unnormalize(r.horizon_homo, width, height, focal, hl);

    //hvps = unnormalize(results.hvp_homo, width, height, focal, 0);
    //std::vector<cv::Vec2d> hvps;
    hvps.reserve(r.hvp_homo.size());
    for (const auto &hvp_homo : r.hvp_homo)
    {
        cv::Vec2d hvp;
        unnormalize(hvp_homo, width, height, focal, hvp);
        hvps.emplace_back(hvp);
    }

    hvp_groups = r.hvp_groups;
    //z = unnormalize(z_homo_cand{best_z_cand}, width, height, focal, 0);
    //cv::Vec2d z;
    unnormalize(z_homo_cand[best_z_cand], width, height, focal, z);
    //std::vector<int> z_group;
    //z_group = z_group_cand{best_z_cand};
    z_group = z_group_cand[best_z_cand];

    if (params.debug_fileid != nullptr)
    {
        fprintf(params.debug_fileid, "z: [%f %f] [%f %f %f]\n", z[0], z[1], z_homo_cand[best_z_cand][0], z_homo_cand[best_z_cand][1], z_homo_cand[best_z_cand][2]);
        fprintf(params.debug_fileid, "z group: ");
        for (int i = 0; i < z_group.size(); ++i)
        {
            fprintf(params.debug_fileid, "%d ", z_group[i]);
        }
        fprintf(params.debug_fileid, "\n");
        auto hvp_count = hvps.size();
        for (int v = 0; v < hvp_count; ++v)
        {
            fprintf(params.debug_fileid, "%d vp: [%f %f]\n", v, hvps[v][0], hvps[v][1]);
        }
        for (int j = 0; j < hvp_groups.size(); ++j)
        {
            fprintf(params.debug_fileid, "%d vp group: ", j);
            for (int k = 0; k < hvp_groups[j].size(); ++k)
            {
                fprintf(params.debug_fileid, "%d ", hvp_groups[j][k]);
            }
            fprintf(params.debug_fileid, "\n");
        }
    }

    double focal_calbrated;
    std::vector<cv::Vec2d> manh_vps;
    calibrate(z, hvps, width, height, focal, manh_vps);
    assert(focal_calbrated > 0);

    cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
    K.at<double>(0, 0) = focal;
    K.at<double>(1, 1) = focal;
    K.at<double>(0, 2) = (double)width / 2;
    K.at<double>(1, 2) = (double)height / 2;

    auto hl_homo = line_hmg_from_two_points(hl[0], hl[1]);
    std::vector<Transform> transforms;
    for (int i = 0; i < hvps.size(); ++i)
    {
        std::vector<Transform> trans;
        orthorectify(img, hvps[i], hvp_groups[i], z, z_group, ls, 4, K, hl_homo, trans);
        transforms.insert(transforms.end(), trans.begin(), trans.end());
    }

    return 0;
}
