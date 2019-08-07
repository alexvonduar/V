#ifndef _LINE_ANGLE2_HPP_
#define _LINE_ANGLE2_HPP_

#include <opencv2/opencv.hpp>

//function angle = line_angle2(line)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
T line_angle2(const cv::Vec<T, 4> &line)
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

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
T line_angle2(const cv::Vec<T, 2> &p1, const cv::Vec<T, 2> &p2)
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
