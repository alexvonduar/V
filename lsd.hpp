#ifndef _LSD_HPP_
#define _LSD_HPP_

#include <opencv2/opencv.hpp>

#if defined(__cplusplus)
extern "C"
{
#endif
#include "util/lsd/lsd.h"
#if defined(__cplusplus)
}
#endif

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void lsd_detect(const cv::Mat &img, std::vector<cv::Vec<T, 4>> &lines)
{
    cv::Mat img_double;
    cv::cvtColor(img, img_double, cv::COLOR_BGR2GRAY);
    img_double.convertTo(img_double, CV_64FC1);
    image_double_s lsd_img_double;
    lsd_img_double.data = (double *)img_double.data;
    lsd_img_double.xsize = img_double.cols;
    lsd_img_double.ysize = img_double.rows;
    ntuple_list out = lsd(&lsd_img_double);
    auto n = out->dim;
    assert(n == 5);
    lines.clear();
    lines.reserve(out->size);
    for (int i = 0; i < out->size; ++i)
    {
        auto x1 = out->values[n * i];
        auto y1 = out->values[n * i + 1];
        auto x2 = out->values[n * i + 2];
        auto y2 = out->values[n * i + 3];
        if ((x2 < x1) or ((x1 == x2) and (y2 < y1)))
        {
            cv::Vec<T, 4> line{x2, y2, x1, y1};
            lines.emplace_back(line);
        }
        else
        {
            cv::Vec<T, 4> line{x1, y1, x2, y2};
            lines.emplace_back(line);
        }
    }

    free_ntuple_list(out);
}

#endif //_LSD_HPP_
