#ifndef _LINE_SIZE_
#define _LINE_SIZE_

#include <opencv2/opencv.hpp>

//function [lsize] = line_size(seglines)
double line_size(const cv::Vec4d& seglines)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto dx = seglines[0] - seglines[2];
    auto dy = seglines[1] - seglines[3];
    auto lsize = std::sqrt(dx * dx + dy * dy);
    return lsize;
    //end
}

double line_size_square(const cv::Vec4d& seglines)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto dx = seglines[0] - seglines[2];
    auto dy = seglines[1] - seglines[3];
    auto lsize = dx * dx + dy * dy;
    return lsize;
    //end
}

double line_size(const cv::Vec2d& p1, const cv::Vec2d& p2)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto d = p1 - p2;
    auto lsize = cv::norm(d);
    return lsize;
    //end
}

double line_size_square(const cv::Vec2d& p1, const cv::Vec2d& p2)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto d = p1 - p2;
    auto lsize = d[0] * d[0] + d[1] * d[1];
    return lsize;
    //end
}

#endif //_LINE_SIZE_
