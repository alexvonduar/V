#ifndef _LINE_TO_SEGMENT_HPP_
#define _LINE_TO_SEGMENT_HPP_

#include <opencv2/opencv.hpp>

//function segments = line_to_segment(lines)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void line_to_segment(const cv::Vec<T, 4> &line, std::vector<cv::Vec<T, 2>> &segments)
{
    // generate segments for plotting out of lines
    // lines: [x1, y1, x2, y2]
    // segments: [x11, y11; x12, y12; nan, nan; ...]

    //segments = nan(size(lines,1)*3, 2);
    //segments(1:3:end,:) = lines(:,[1,2]);
    //segments(2:3:end,:) = lines(:,[3,4]);
    segments = std::vector<cv::Vec<T, 2>>{cv::Vec<T, 2>{line[0], line[1]}, cv::Vec<T, 2>{line[2], line[3]}};
}

#endif //_LINE_TO_SEGMENT_HPP_
