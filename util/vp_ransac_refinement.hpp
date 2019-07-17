#ifndef _VP_RANSAC_REFINEMENT_HPP_
#define _VP_RANSAC_REFINEMENT_HPP_

#include <opencv2/opencv.hpp>
#include "ransac_intersection.hpp"
#include "default_params.hpp"
#include "lines_normal.hpp"

//function [zenith_homo, inlierId] = vp_ransac_refinement(lines_homo, opt) % [Zhai et al. 2016]
void vp_ransac_refinement(const std::vector<cv::Vec3d> &lines_homo, const Params &opt, cv::Vec3d &zenith_homo, std::vector<int> &inlierId)
{
    //option = struct();
    RANSAC_Parameter option;
    option.iterNum = 50;
    option.thInlrRatio = .02;
    //option.thDist = sind(opt.theta_con);
    option.thDist = std::sin(opt.theta_con * CV_PI / 180.);
    cv::Vec3d M;
    //--std::vector<int> inlierId;
    //[~, inlierId] = ransac_intersection(lines_homo, option);
    ransac_intersection(lines_homo, option, M, inlierId);
    std::vector<cv::Vec3d> inlier_lines;
    inlier_lines.reserve(lines_homo.size());
    for (const auto &i : inlierId)
    {
        inlier_lines.emplace_back(lines_homo[i]);
    }
    //zenith_homo = lines_normal(lines_homo(:,inlierId));
    lines_normal(inlier_lines, cv::Mat(), zenith_homo);
}

#endif //_VP_RANSAC_REFINEMENT_HPP_
