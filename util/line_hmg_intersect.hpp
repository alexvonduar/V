#ifndef _LINE_HMG_INTERSECT_HPP_
#define _LINE_HMG_INTERSECT_HPP_

#include <opencv2/opencv.hpp>

//function p = line_hmg_intersect(l1_hmg,l2_hmg)
static inline void line_hmg_intersect(
    const cv::Vec3d &l1_hmg,
    const cv::Vec3d &l2_hmg,
    cv::Vec2d &p)
{
    //p_hmg = cross(l1_hmg,l2_hmg);
    auto p_hmg = l1_hmg.cross(l2_hmg);
    //p(1,1) = p_hmg(1)/p_hmg(3);
    //p(1,2) = p_hmg(2)/p_hmg(3);
    p[0] = p_hmg[0] / p_hmg[2];
    p[1] = p_hmg[1] / p_hmg[2];
    //end
}

#endif //_LINE_HMG_INTERSECT_HPP_
