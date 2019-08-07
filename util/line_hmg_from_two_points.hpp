#ifndef _LINE_HMG_FROM_TWO_POINTS_HPP_
#define _LINE_HMG_FROM_TWO_POINTS_HPP_

#include <opencv2/opencv.hpp>

//function l = line_hmg_from_two_points(p1,p2)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline cv::Vec<T, 3> line_hmg_from_two_points(const cv::Vec<T, 2> &p1, const cv::Vec<T, 2> &p2)
{
    //v1 = double([p1(1);p1(2);1]);
    cv::Vec<T, 3> v1{p1[0], p1[1], 1};
    //v2 = double([p2(1);p2(2);1]);
    cv::Vec<T, 3> v2{p2[0], p2[1], 1};
    auto l = v1.cross(v2);
    auto norm = std::sqrt(l[0] * l[0] + l[1] * l[1]);
    l[0] /= norm;
    l[1] /= norm;
    l[2] /= norm;
    return l;
}

//function l = line_hmg_from_two_points(p1,p2)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline cv::Vec<T, 3> line_hmg_from_two_points(const cv::Vec<T, 4> &line)
{
    return line_hmg_from_two_points(cv::Vec<T, 2>{line[0], line[1]}, cv::Vec<T, 2>{line[2], line[3]});
}

#endif //_LINE_HMG_FROM_TWO_POINTS_HPP_
