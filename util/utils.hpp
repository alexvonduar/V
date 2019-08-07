#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <opencv2/opencv.hpp>

namespace VUtils
{

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline bool isContourConvex_(const cv::Point_<T> *p, int n)
{
    cv::Point_<T> prev_pt = p[(n - 2 + n) % n];
    cv::Point_<T> cur_pt = p[n - 1];

    T dx0 = cur_pt.x - prev_pt.x;
    T dy0 = cur_pt.y - prev_pt.y;
    int orientation = 0;

    for (int i = 0; i < n; i++)
    {
        T dxdy0, dydx0;
        T dx, dy;

        prev_pt = cur_pt;
        cur_pt = p[i];

        dx = cur_pt.x - prev_pt.x;
        dy = cur_pt.y - prev_pt.y;
        dxdy0 = dx * dy0;
        dydx0 = dy * dx0;

        // find orientation
        // orient = -dy0 * dx + dx0 * dy;
        // orientation |= (orient > 0) ? 1 : 2;
        orientation |= (dydx0 > dxdy0) ? 1 : ((dydx0 < dxdy0) ? 2 : 3);
        if (orientation == 3)
            return false;

        dx0 = dx;
        dy0 = dy;
    }

    return true;
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline bool isContourConvex_(const cv::Vec<T, 2> *p, int n)
{
    auto prev_pt = p[(n - 2 + n) % n];
    auto cur_pt = p[n - 1];

    T dx0 = cur_pt[0] - prev_pt[0];
    T dy0 = cur_pt[1] - prev_pt[1];
    int orientation = 0;

    for (int i = 0; i < n; i++)
    {
        T dxdy0, dydx0;
        T dx, dy;

        prev_pt = cur_pt;
        cur_pt = p[i];

        dx = cur_pt[0] - prev_pt[0];
        dy = cur_pt[1] - prev_pt[1];
        dxdy0 = dx * dy0;
        dydx0 = dy * dx0;

        // find orientation
        // orient = -dy0 * dx + dx0 * dy;
        // orientation |= (orient > 0) ? 1 : 2;
        orientation |= (dydx0 > dxdy0) ? 1 : ((dydx0 < dxdy0) ? 2 : 3);
        if (orientation == 3)
            return false;

        dx0 = dx;
        dy0 = dy;
    }

    return true;
}

#if 0
static inline bool isContourConvex(cv::InputArray _contour)
{
    cv::Mat contour = _contour.getMat();
    int total = contour.checkVector(2), depth = contour.depth();
    CV_Assert(total >= 0 && (depth == CV_32F || depth == CV_32S));

    if (total == 0)
        return false;

    return depth == CV_32S ? isContourConvex_(contour.ptr<cv::Point>(), total) : isContourConvex_(contour.ptr<cv::Point2f>(), total);
}
#endif

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline bool is_convex_polygon(const std::vector<cv::Point_<T>> &polygon)
{
    assert(polygon.size() >= 3);
    auto p = &polygon[0];
    auto N = polygon.size();
    return isContourConvex_(p, N);
}

template <typename T, std::size_t N, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline bool is_convex_polygon(const std::array<cv::Point_<T>, N> &polygon)
{
    assert(polygon.size() >= 3);
    auto p = &polygon[0];
    return isContourConvex_(p, N);
}

template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline bool is_convex_polygon(const std::vector<cv::Vec<T, 2>> &polygon)
{
    assert(polygon.size() >= 3);
    auto p = &polygon[0];
    auto N = polygon.size();
    return isContourConvex_(p, N);
}

template <typename T, std::size_t N, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline bool is_convex_polygon(const std::array<cv::Vec<T, 2>, N> &polygon)
{
    assert(polygon.size() >= 3);
    auto p = &polygon[0];
    return isContourConvex_(p, N);
}

static inline float identity_similarity(const cv::Mat &hom)
{
    std::vector<cv::Point2f> input{cv::Point2f{0., 0.}, cv::Point2f{640., 0.}, cv::Point2f{640., 640.}, cv::Point2f{0., 640.}};
    std::vector<cv::Point2f> output(input.size());
    cv::perspectiveTransform(input, output, hom);
    if (is_convex_polygon(output))
    {
        cv::Point2f offset = output[0];
        for (auto &pt : output)
        {
            pt -= offset;
        }
        std::vector<cv::Point2f> points;
        auto area = cv::intersectConvexConvex(input, output, points);
        return area;
    }
    else
    {
        return -1.;
    }
}

static inline bool morphological_check(const cv::Mat &hom, const cv::Rect &_rect)
{
    auto rect = _rect;
    if (rect.area() == 0)
    {
        rect = cv::Rect(0, 0, 640, 640);
    }
    std::vector<cv::Point2f> input{
        cv::Point2f{float(rect.x), float(rect.y)},
        cv::Point2f{float(rect.x + rect.width), float(rect.y)},
        cv::Point2f{float(rect.x + rect.width), float(rect.y + rect.height)},
        cv::Point2f{float(rect.x), float(rect.y + rect.height)}};
    std::vector<cv::Point2f> output(input.size());
    cv::perspectiveTransform(input, output, hom);
    /*
    std::cout << "hom: " << hom << std::endl;
    for (int i = 0; i < input.size(); ++i) {
        ZLOGI("[%f, %f] -> [%f, %f]\n", input[i].x, input[i].y, output[i].x, output[i].y);
    }
    */
    /*
    cv::Point2f offset = output[0];

    for (auto& pt : output) {
        pt -= offset;
        //ZLOGI("transformed %f %f\n", pt.x, pt.y);
    }
    */
    bool is_convex = is_convex_polygon(output);
    auto rrect = cv::minAreaRect(output);
    auto aspect_ratio = rrect.size.width > rrect.size.height ? (rrect.size.width / rrect.size.height) : (rrect.size.height / rrect.size.width);
    auto valid_aspect_ratio = aspect_ratio < 20;
    auto valid_scale = false;
    auto valid_offset = false;
    if (is_convex)
    {
        auto original_area = rect.area();
        auto warped_area = rrect.size.area();
        auto scale = (warped_area > original_area) ? (warped_area / original_area) : (original_area / warped_area);
        printf("area: %d %f scale: %f\n", original_area, warped_area, scale);
        valid_scale = scale < 49;
        auto offset_x = std::abs(rrect.center.x);
        auto offset_scale = offset_x / rect.width;
        auto valid_offset_x = offset_scale < 80;
        auto offset_y = std::abs(rrect.center.y);
        offset_scale = offset_y / rect.height;
        printf("offset: %f %f\n", offset_x, offset_y);
        auto valid_offset_y = offset_scale < 80;
        valid_offset = valid_offset_x and valid_offset_y;
        if (!valid_scale or !valid_offset)
        {
            printf("scale: %f offset: [%f %f]\n", scale, offset_x, offset_y);
        }
    }

    auto valid = is_convex and valid_aspect_ratio and valid_scale and valid_offset;
    //if (!valid) {
    printf("is convex: %d aspect ratio: %f[valid %d] scale: %d offset: %d\n", is_convex, aspect_ratio, valid_aspect_ratio, valid_scale, valid_offset);
    //}
    return valid;
}

bool is_good_homography(const cv::Mat &hom, const cv::Rect &rect, const bool &do_morphological_check = true)
{
    if (hom.empty())
    {
        return false;
    }

    cv::Mat sub_homo;
    hom.convertTo(sub_homo, CV_32F);

    cv::Rect roi(0, 0, 2, 2);
    sub_homo = hom(roi);

    auto valid_diag = sub_homo.at<float>(0, 0) > 0 and sub_homo.at<float>(1, 1) > 0;
    auto det = cv::determinant(sub_homo);
    if (valid_diag and det > 0)
    {
        // do morphology test
        if (do_morphological_check)
        {
            if (morphological_check(hom, rect) and morphological_check(hom, cv::Rect{0, 0, 640, 640}))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return true;
        }
    }
    else
    {
        return false;
    }
}

}; // namespace VUtils

#endif //_UTILS_HPP_
