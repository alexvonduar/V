#ifndef _LINE_ANGLE2_HPP_
#define _LINE_ANGLE2_HPP_

#include <opencv2/opencv.hpp>

//function angle = line_angle2(line)
double line_angle2(const cv::Vec4d &line)
{
    auto x1 = line[0];
    auto y1 = line[1];
    auto x2 = line[2];
    auto y2 = line[3];
    if (x2 > x1)
        return std::atan2((y2 - y1), (x2 - x1));
    else
        return std::atan2((y1 - y2), (x1 - x2));
    //end
}

double line_angle2(const cv::Vec2d &p1, const cv::Vec2d &p2)
{
    auto x1 = p1[0];
    auto y1 = p1[1];
    auto x2 = p2[0];
    auto y2 = p2[1];
    if (x2 > x1)
        return std::atan2((y2 - y1), (x2 - x1));
    else
        return std::atan2((y1 - y2), (x1 - x2));
    //end
}

#endif //_LINE_ANGLE2_HPP_
