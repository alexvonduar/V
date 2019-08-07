#ifndef _LINE_SIZE_
#define _LINE_SIZE_

#include <opencv2/opencv.hpp>

//function [lsize] = line_size(seglines)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
T line_size(const cv::Vec<T, 4> &seglines)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto dx = seglines[0] - seglines[2];
    auto dy = seglines[1] - seglines[3];
    auto lsize = std::sqrt(dx * dx + dy * dy);
    return lsize;
    //end
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
T line_size_square(const cv::Vec<T, 4> &seglines)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto dx = seglines[0] - seglines[2];
    auto dy = seglines[1] - seglines[3];
    auto lsize = dx * dx + dy * dy;
    return lsize;
    //end
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
T line_size(const cv::Vec<T, 2> &p1, const cv::Vec<T, 2> &p2)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto d = p1 - p2;
    auto lsize = cv::norm(d);
    return lsize;
    //end
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
T line_size_square(const cv::Vec<T, 2> &p1, const cv::Vec<T, 2> &p2)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto d = p1 - p2;
    auto lsize = d[0] * d[0] + d[1] * d[1];
    return lsize;
    //end
}

#endif //_LINE_SIZE_
