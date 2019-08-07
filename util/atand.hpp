#ifndef _ATAND_HPP_
#define _ATAND_HPP_

#include <opencv2/opencv.hpp>

#if 0
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

#else

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value, T>::type>
static inline T atan(const T &y, const T &x)
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

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value, T>::type>
static inline T atand(const T &y, const T &x)
{
    auto r = atan(y, x);
    return r * 180. / CV_PI;
}

#endif

#endif //_ATAND_HPP_
