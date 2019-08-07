#ifndef _FILTER_VERHOR_LINES_HPP_
#define _FILTER_VERHOR_LINES_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"

//function lines_id = filter_verhor_lines(ls_homo, z_homo,  params) % as in [Zhai et al. 2016]
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void filter_verhor_lines(const std::vector<cv::Vec<T, 3>> &ls_homo,
                                       const cv::Vec<T, 3> &z_homo,
                                       const Params &params,
                                       std::vector<int> &lines_id)
{
    // filter vertical line segments
    //cos_val = abs(ls_homo'*z_homo);
    //inlier_id = cos_val > sind(params.theta_verline);
    //lines_id = find(inlier_id);
    lines_id.reserve(ls_homo.size());
    T sind = std::sin(params.theta_verline * CV_PI / 180.);
    for (int i = 0; i < ls_homo.size(); ++i)
    {
        T cos_val = std::abs(ls_homo[i].dot(z_homo));
        if (cos_val > sind)
        {
            lines_id.emplace_back(i);
        }
    }

    // filter horizontal line segments
    //if ~params.include_infinite_hvps
    if (!params.include_infinite_hvps)
    {
        //ortho_thres = sind(params.theta_horline);
        T ortho_thres = std::sin(params.theta_horline * CV_PI / 180.);
        //lhomo = ls_homo(:,lines_id);
        //cos_val = abs(lhomo'*[-z_homo(2); z_homo(1); 0]);
        //inlier_id = cos_val > ortho_thres;
        //lines_id = lines_id(inlier_id);
        std::vector<int> inlier_id;
        inlier_id.reserve(lines_id.size());
        cv::Vec3d z{-z_homo[1], z_homo[0], 0};
        for (const auto &id : lines_id)
        {
            auto cos_val = std::abs(ls_homo[id].dot(z));
            if (cos_val > ortho_thres)
            {
                inlier_id.emplace_back(id);
            }
        }
        lines_id = inlier_id;
    }
    //end
}

#endif //_FILTER_VERHOR_LINES_HPP_
