#ifndef _VP_REFINEMENT_HPP_
#define _VP_REFINEMENT_HPP_

#include <opencv2/opencv.hpp>
#include "default_params.hpp"
#include "util/lines_normal.hpp"

//function vps_homo = vp_refinement(lines_homo, vps_homo, horizon_homo, params) % from [Zhai et al.]
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void vp_refinement(
    const std::vector<cv::Vec<T, 3>> &lines_homo,
    std::vector<cv::Vec<T, 3>> &vps_homo,
    const cv::Vec<T, 3> &horizon_homo,
    const Params &params)
{
    //rng(1)
    auto niters = params.refine_niters;
    //inlier_thr = sind(params.theta_con);
    auto inlier_thr = std::sin(params.theta_con * CV_PI / 180.);
    //ngrps = size(vps_homo,2);
    auto ngrps = vps_homo.size();
    //for ig = 1:ngrps
    for (int ig = 0; ig < ngrps; ++ig)
    {
        //for it = 1:niters
        for (int it = 0; it < niters; ++it)
        {
            //good_ids = find(abs(vps_homo(:,ig)'*lines_homo) < inlier_thr)';
            std::vector<cv::Vec<T, 3>> good_lines;
            good_lines.reserve(lines_homo.size());
            for (const auto &l : lines_homo)
            {
                auto t = std::abs(vps_homo[ig].dot(l));
                if (t < inlier_thr)
                {
                    good_lines.emplace_back(l);
                }
            }
            //vps_homo(:,ig) = lines_normal(lines_homo(:,good_ids), horizon_homo);
            lines_normal(good_lines, horizon_homo, params, vps_homo[ig]);

            //end
        }
        //end
    }
}

#endif //_VP_REFINEMENT_HPP_
