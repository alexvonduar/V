#ifndef _LINE_SIZE_
#define _LINE_SIZE_

#include <opencv2/opencv.hpp>

//function [lsize] = line_size(seglines)
double line_size(cv::Vec4d seglines)
{
    //lsize = sqrt((seglines(:,1)-seglines(:,3)).^2+(seglines(:,2)-seglines(:,4)).^2);
    auto dx = seglines[0] - seglines[2];
    auto dy = seglines[1] - seglines[3];
    auto lsize = std::sqrt(dx * dx + dy * dy);
    return lsize;
    //end
}

#endif //_LINE_SIZE_
