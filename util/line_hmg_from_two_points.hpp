#ifndef _LINE_HMG_FROM_TWO_POINTS_HPP_
#define _LINE_HMG_FROM_TWO_POINTS_HPP_

#include <opencv2/opencv.hpp>

//function l = line_hmg_from_two_points(p1,p2)
static inline cv::Vec3d line_hmg_from_two_points(const cv::Vec2d &p1, const cv::Vec2d &p2)
{
    //v1 = double([p1(1);p1(2);1]);
    cv::Vec3d v1{p1[0], p1[1], 1};
    //v2 = double([p2(1);p2(2);1]);
    cv::Vec3d v2{p2[0], p2[1], 1};
    auto l = v1.cross(v2);
    auto norm = std::sqrt(l[0] * l[0] + l[1] * l[1]);
    l[0] /= norm;
    l[1] /= norm;
    l[2] /= norm;
    return l;
}

//function l = line_hmg_from_two_points(p1,p2)
static inline cv::Vec3d line_hmg_from_two_points(const cv::Vec4d &line)
{
    return line_hmg_from_two_points(cv::Vec2d{line[0], line[1]}, cv::Vec2d{line[2], line[3]});
}

#endif //_LINE_HMG_FROM_TWO_POINTS_HPP_
