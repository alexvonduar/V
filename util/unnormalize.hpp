#ifndef _UNNORMALIZE_HPP_
#define _UNNORMALIZE_HPP_

#include <opencv2/opencv.hpp>

#include "line_to_segment.hpp"

/// convert lines or points from homogeneous system to cartesian system (image frame)
/// according to [Zhai et al., CVPR'2016]
//function geometry_img = unnormalize(geometry_homo, width, height, focal, isline)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void unnormalize(cv::Vec<T, 3> geometry_homo,
                               const int &width,
                               const int &height,
                               const T &focal,
                               cv::Vec<T, 2> &geometry_img)
{
    //p = geometry_homo;
    //geometry_img = bsxfun(@rdivide, p(1:2, :), p(3,:))';
    //geometry_img = bsxfun(@plus, (geometry_img * focal), [width, height]/2);
    geometry_img = cv::Vec<T, 2>{geometry_homo[0] / geometry_homo[2], geometry_homo[1] / geometry_homo[2]};
    geometry_img *= focal;
    geometry_img += cv::Vec<T, 2>{T(width) / 2, T(height) / 2};
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void unnormalize(cv::Vec<T, 3> geometry_homo,
                               const int &width,
                               const int &height,
                               const T &focal,
                               std::vector<cv::Vec<T, 2>> &geometry_img)
{

    auto h = geometry_homo;
    auto pl = (0 - width / 2) / focal;
    auto pr = (width - width / 2) / focal;
    auto ly = T(height) / 2 - ((h[0] * pl + h[2]) / h[1]) * focal;
    auto ry = T(height) / 2 - ((h[0] * pr + h[2]) / h[1]) * focal;
    //geometry_img = [zeros(size(ly)), ly, width*ones(size(ry)), ry];
    cv::Vec<T, 4> line{0, ly, (T)width, ry};
    //geometry_img = line_to_segment(geometry_img);
    line_to_segment(line, geometry_img);
    //end
}

#endif //_UNNORMALIZE_HPP_
