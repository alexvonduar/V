#ifndef _COORDINATES_CENTER_HPP_
#define _COORDINATES_CENTER_HPP_

#include <opencv2/opencv.hpp>

//function centered_coordinates = coordinates_center(coordinates, width, height, focal)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void coordinates_center(
    const std::vector<cv::Vec<T, 2>> &coordinates,
    const int &width,
    const int &height,
    const double &focal,
    std::vector<cv::Vec<T, 2>> &centered_coordinates)
{
    centered_coordinates.reserve(coordinates.size());
    //center = [width, height] / 2;
    cv::Vec<T, 2> center{double(width) / 2, double(height) / 2};
    //if size(coordinates, 2) == 2
    //  centered_coordinates = bsxfun(@minus, coordinates, center) / focal;
    //else
    //  centered_coordinates = bsxfun(@minus, coordinates, [center, center]) / focal;
    //end
    for (const auto &v : coordinates)
    {
        centered_coordinates.emplace_back(v - center);
    }
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void coordinates_center(
    const std::vector<cv::Vec<T, 4>> &coordinates,
    const int &width,
    const int &height,
    const T &focal,
    std::vector<cv::Vec<T, 4>> &centered_coordinates)
{
    centered_coordinates.reserve(coordinates.size());
    //center = [width, height] / 2;
    cv::Vec<T, 4> center{double(width) / 2, double(height) / 2, double(width) / 2, double(height) / 2};
    //if size(coordinates, 2) == 2
    //  centered_coordinates = bsxfun(@minus, coordinates, center) / focal;
    //else
    //  centered_coordinates = bsxfun(@minus, coordinates, [center, center]) / focal;
    //end
    for (const auto &v : coordinates)
    {
        centered_coordinates.emplace_back(v - center);
    }
}

#endif //_COORDINATES_CENTER_HPP_
