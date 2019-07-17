#ifndef _UNNORMALIZE_HPP_
#define _UNNORMALIZE_HPP_

#include <opencv2/opencv.hpp>

#include "line_to_segment.hpp"

/// convert lines or points from homogeneous system to cartesian system (image frame)
/// according to [Zhai et al., CVPR'2016]
//function geometry_img = unnormalize(geometry_homo, width, height, focal, isline)
static inline void unnormalize(cv::Vec3d geometry_homo,
                               const int &width,
                               const int &height,
                               const double &focal,
                               cv::Vec2d &geometry_img)
{
    //p = geometry_homo;
    //geometry_img = bsxfun(@rdivide, p(1:2, :), p(3,:))';
    //geometry_img = bsxfun(@plus, (geometry_img * focal), [width, height]/2);
    geometry_img = cv::Vec2d{geometry_homo[0] / geometry_homo[2], geometry_homo[1] / geometry_homo[2]};
    geometry_img *= focal;
    geometry_img += cv::Vec2d{double(width) / 2, double(height) / 2};
}

static inline void unnormalize(cv::Vec3d geometry_homo,
                               const int &width,
                               const int &height,
                               const double &focal,
                               std::vector<cv::Vec2d> &geometry_img)
{

    auto h = geometry_homo;
    auto pl = (0 - width / 2) / focal;
    auto pr = (width - width / 2) / focal;
    auto ly = double(height) / 2 - ((h[0] * pl + h[2]) / h[1]) * focal;
    auto ry = double(height) / 2 - ((h[0] * pr + h[2]) / h[1]) * focal;
    //geometry_img = [zeros(size(ly)), ly, width*ones(size(ry)), ry];
    cv::Vec4d line{0, ly, (double)width, ry};
    //geometry_img = line_to_segment(geometry_img);
    line_to_segment(line, geometry_img);
    //end
}

#endif //_UNNORMALIZE_HPP_
