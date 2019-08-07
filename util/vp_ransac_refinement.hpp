#ifndef _VP_RANSAC_REFINEMENT_HPP_
#define _VP_RANSAC_REFINEMENT_HPP_

#include <opencv2/opencv.hpp>
#include "ransac_intersection.hpp"
#include "default_params.hpp"
#include "lines_normal.hpp"

//function [zenith_homo, inlierId] = vp_ransac_refinement(lines_homo, opt) % [Zhai et al. 2016]
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
void vp_ransac_refinement(const std::vector<cv::Vec<T, 3>> &lines_homo, const Params &opt, cv::Vec<T, 3> &zenith_homo, std::vector<int> &inlierId)
{
    //option = struct();
    struct _RANSAC_Parameter<T> option;
    option.iterNum = 50;
    option.thInlrRatio = .02;
    //option.thDist = sind(opt.theta_con);
    option.thDist = std::sin(opt.theta_con * CV_PI / 180.);
    cv::Vec<T, 3> M;
    //--std::vector<int> inlierId;
    //[~, inlierId] = ransac_intersection(lines_homo, option);
    ransac_intersection(lines_homo, option, M, inlierId);
    std::vector<cv::Vec<T, 3>> inlier_lines;
    inlier_lines.reserve(lines_homo.size());
    for (const auto &i : inlierId)
    {
        inlier_lines.emplace_back(lines_homo[i]);
    }
    if (opt.debug_fileid != nullptr)
    {
        fprintf(opt.debug_fileid, "ransac inlier %ld: ", inlierId.size());
        for (const auto &i : inlierId)
        {
            fprintf(opt.debug_fileid, "%d ", i);
        }
        fprintf(opt.debug_fileid, "\n");
    }
    //zenith_homo = lines_normal(lines_homo(:,inlierId));
    lines_normal(inlier_lines, cv::Mat_<T>(), opt, zenith_homo);
}

#endif //_VP_RANSAC_REFINEMENT_HPP_
