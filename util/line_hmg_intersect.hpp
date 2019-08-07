#ifndef _LINE_HMG_INTERSECT_HPP_
#define _LINE_HMG_INTERSECT_HPP_

#include <opencv2/opencv.hpp>

//function p = line_hmg_intersect(l1_hmg,l2_hmg)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline cv::Vec<T, 2> line_hmg_intersect(
    const cv::Vec<T, 3> &l1_hmg,
    const cv::Vec<T, 3> &l2_hmg)
{
    //p_hmg = cross(l1_hmg,l2_hmg);
    auto p_hmg = l1_hmg.cross(l2_hmg);
    //p(1,1) = p_hmg(1)/p_hmg(3);
    //p(1,2) = p_hmg(2)/p_hmg(3);
    cv::Vec<T, 2> p{p_hmg[0] / p_hmg[2], p_hmg[1] / p_hmg[2]};
    return p;
    //-p[0] = p_hmg[0] / p_hmg[2];
    //-p[1] = p_hmg[1] / p_hmg[2];
    //end
}

#endif //_LINE_HMG_INTERSECT_HPP_
