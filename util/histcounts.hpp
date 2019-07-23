#ifndef _HISTCOUNTS_HPP_
#define _HISTCOUNTS_HPP_

#include <opencv2/opencv.hpp>

// implement MATLAB like histcounts function with
// OpenCV calcHist function

template <typename T>
static inline void histcounts(
    const std::vector<T> &v,
    const int &nbins,
    const float *default_range,
    cv::MatND &hist,
    std::vector<double> &edges)
{
    edges.clear();
    edges.reserve(nbins + 1);
    if (v.size() == 0)
    {
        hist.release();
        return;
    }
    cv::Mat src(1, v.size(), CV_32FC1);
    for (int i = 0; i < v.size(); ++i)
    {
        src.at<float>(0, i) = v[i];
    }

    int channels[] = {0};
    int histSize[] = {nbins};
    float range[] = {0., 0.};
    const float *ranges[] = {range};
    double delta;
    double min, max;
    if (default_range == nullptr)
    {
        auto minmax = std::minmax_element(v.begin(), v.end());
        min = *minmax.first;
        max = *minmax.second;
        delta = max - min;
        delta /= nbins;
    }
    else
    {
        min = default_range[0];
        max = default_range[1];
        delta = max - min;
        delta /= nbins;
    }
    range[0] = min;
    range[1] = max;

    cv::calcHist(&src, 1, channels, cv::Mat(), // do not use mask
                 hist, 1, histSize, ranges,
                 true, // the histogram is uniform
                 false);

    hist = hist.reshape(1, hist.cols);

    hist.convertTo(hist, CV_64FC1);
    auto shist = cv::sum(hist);
    if (shist[0] + 1 == v.size()) {
        hist.at<double>(0, nbins - 1) += 1;
    }
    assert(cv::sum(hist)[0] == v.size());

    for (int i = 0; i < nbins; ++i)
    {
        edges.emplace_back(min + delta * i);
    }
    edges.emplace_back(max);
}

#endif //_HISTCOUNTS_HPP_
