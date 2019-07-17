#ifndef _ATAND_HPP_
#define _ATAND_HPP_

#include <opencv2/opencv.hpp>

static inline double atan(const double &y, const double &x)
{
    auto r = std::atan2(y, x);
    if (r > (CV_PI / 2))
    {
        r -= CV_PI;
    }
    else if (r < -(CV_PI / 2))
    {
        r += CV_PI;
    }
    return r;
}

static inline double atand(const double &y, const double &x)
{
    auto r = atan(y, x);
    return r * 180 / CV_PI;
}

#endif //_ATAND_HPP_
