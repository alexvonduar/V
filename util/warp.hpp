#ifndef _WARP_HPP_
#define _WARP_HPP_

#include <opencv2/opencv.hpp>

#include "utils.hpp"

static inline void align_rect(const cv::Mat &H, const cv::Rect &rect, cv::Mat &shift, cv::Size &sz, const cv::Size &max_sz = cv::Size())
{
    std::vector<cv::Point2f> corners{
        cv::Point2f{(float)rect.x, (float)rect.y},
        cv::Point2f{(float)rect.x + rect.width, (float)rect.y},
        cv::Point2f{(float)rect.x + rect.width, (float)rect.y + rect.height},
        cv::Point2f{(float)rect.x, (float)rect.y + rect.height}};
    std::vector<cv::Point2f> warped_corners(4);
    cv::perspectiveTransform(corners, warped_corners, H);
    auto _r = cv::boundingRect(warped_corners);
    sz = _r.size();

    float scale = 1.0;
    if (max_sz.width)
    {
        auto s = (float)max_sz.width / sz.width;
        if (s < scale)
        {
            scale = s;
        }
    }
    if (max_sz.height)
    {
        auto s = (float)max_sz.height / sz.height;
        if (s < scale)
        {
            scale = s;
        }
    }
    if (scale != 1.0)
    {
        for (auto &p : warped_corners)
        {
            p -= cv::Point2f{(float)_r.x, (float)_r.y};
            p *= scale;
        }
        shift = cv::findHomography(corners, warped_corners, cv::noArray(), cv::LMEDS);
        shift.convertTo(shift, CV_32F);
        _r = cv::boundingRect(warped_corners);
        sz = _r.size();
    }
    else
    {
        shift = cv::Mat::eye(3, 3, CV_32F);
        shift.at<float>(0, 2) = -_r.x;
        shift.at<float>(1, 2) = -_r.y;
        shift = shift * H;
    }
}

static inline void warp(const cv::Mat &input, cv::Mat &output, const cv::Mat &H, const cv::Size max_sz = cv::Size())
{
    cv::Mat _H;
    H.convertTo(_H, CV_32F);
    auto sz = input.size();
    auto rect = cv::Rect{0, 0, sz.width, sz.height};
    if (VUtils::is_good_homography(_H, rect))
    {
        cv::Mat shift;
        cv::Size sz;
        align_rect(_H, rect, shift, sz, max_sz);
        output = cv::Mat::zeros(sz, input.type());
        cv::warpPerspective(input, output, shift, sz, cv::INTER_LANCZOS4, cv::BORDER_TRANSPARENT);
    }
    else
    {
        output.release();
    }
}

#endif //_WARP_HPP_
