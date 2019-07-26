#ifndef _ORTHORECTIFY_HPP_
#define _ORTHORECTIFY_HPP_

#include <opencv2/opencv.hpp>

#include "line_hmg_from_two_points.hpp"
#include "line_hmg_intersect.hpp"
#include "line_size.hpp"
#include "line_angle2.hpp"

typedef struct
{
    cv::Mat K;
    cv::Mat R;
    cv::Mat H;
} Transform;

void orthorectify(
    const cv::Mat &img,
    const cv::Vec2d &hvp,
    const std::vector<int> &hvp_group,
    const cv::Vec2d &zenith,
    const std::vector<int> &z_group,
    const std::vector<cv::Vec4d> &lines,
    int n_lines_min,
    const cv::Mat &K,
    const cv::Vec3d &horizon_line,
    std::vector<Transform> &transforms)
{
    if (K.empty())
    {
        return;
    }

    auto width = img.cols;
    auto height = img.rows;
    auto n_lines = lines.size();
    std::vector<int> vp_association(n_lines, -1);
    for (const auto &i : z_group)
    {
        vp_association[i] = n_lines;
    }
    for (const auto &i : hvp_group)
    {
        vp_association[i] = i;
    }

    cv::Vec3d hvp_hom{hvp[0], hvp[1], 1};
    auto Ki = K.inv();
    cv::Vec2d pp{K.at<double>(0, 2), K.at<double>(1, 2)};
    auto pp_zen_line = line_hmg_from_two_points(pp, zenith);
    cv::Vec3d vp_zen{zenith[0], zenith[1], 1};
    cv::Mat y = Ki * vp_zen.t();
    y /= cv::norm(y);
    auto t = horizon_line.dot(vp_zen);
    if (t < 0)
    {
        y = -y;
    }

    auto vp_zen_line = line_hmg_from_two_points(hvp, zenith);
    if (vp_zen_line[0] < 0)
    {
        vp_zen_line = -vp_zen_line;
    }
    cv::Vec2d zen_line_normal{vp_zen_line[0], vp_zen_line[1]};
    zen_line_normal /= cv::norm(zen_line_normal);

    std::vector<cv::Vec2d> centroids;
    centroids.reserve(n_lines);
    std::vector<cv::Vec2d> lines_dirs;
    lines_dirs.reserve(n_lines);

    cv::Vec2d proj_min{1, 1};
    cv::Vec2d proj_max{-1, -1};
    cv::Vec2i n_lines_zen{0, 0};
    cv::Vec2i idmax{0, 0};
    cv::Vec2i idmin{0, 0};

    int i = 0;
    for (const auto &line : lines)
    {
        cv::Vec2d p1{line[0], line[1]};
        cv::Vec2d p2{line[2], line[3]};
        auto centroid = p1 + p2;
        centroid /= 2;
        auto lines_dir = p1 - p2;
        lines_dir /= cv::norm(lines_dir);
        auto bundle_dir = centroid - zenith;
        bundle_dir /= cv::norm(bundle_dir);
        auto proj = zen_line_normal.dot(bundle_dir);
        centroids.emplace_back(centroid);
        lines_dirs.emplace_back(lines_dir);
        if (vp_association[i] == n_lines)
        { // z support lines;
            if (proj < 0)
            {
                ++n_lines_zen[0];
                if (proj > proj_max[0])
                {
                    proj_max[0] = proj;
                    idmax[0] = i;
                }
                if (proj < proj_min[0])
                {
                    proj_min[0] = proj;
                    idmin[0] = i;
                }
            }
            else if (proj > 0)
            {
                ++n_lines_zen[1];
                if (proj > proj_max[1])
                {
                    proj_max[1] = proj;
                    idmax[1] = i;
                }
                if (proj < proj_min[1])
                {
                    proj_min[1] = proj;
                    idmin[1] = i;
                }
            }
        }
        ++i;
    }

    for (int s = 0; s < 3; ++s)
    {
        if (n_lines_zen[s] >= n_lines_min)
        {
            auto centroid_min = centroids[idmin[s]];
            auto centroid_max = centroids[idmax[s]];
            auto dist_min = line_size_square(hvp, centroid_min);
            auto dist_max = line_size_square(hvp, centroid_max);
            if (dist_min < dist_max)
            {
                auto middle = hvp + centroid_max;
                middle /= 2;
                auto dist_middle = line_size_square(hvp, middle);
                if (dist_min < dist_middle and dist_middle < dist_max)
                {
                    centroid_min = middle;
                }
            }
            else
            {
                auto middle = hvp + centroid_min;
                middle /= 2;
                auto dist_middle = line_size_square(hvp, middle);
                if (dist_max < dist_middle and dist_middle < dist_min)
                {
                    centroid_max = middle;
                }
            }

            auto lzmin = line_hmg_from_two_points(zenith, centroid_min);
            auto lzmax = line_hmg_from_two_points(zenith, centroid_max);
            auto lzmin_v = centroid_min - zenith;
            auto lzmax_v = centroid_max - zenith;
            lzmin_v /= cv::norm(lzmin_v);
            lzmax_v /= cv::norm(lzmax_v);
            auto hvp_amin = 2 * CV_PI;
            auto hvp_amax = -2 * CV_PI;
            int n_lines_hvp = 0;
            auto idmin2 = 0;
            auto idmax2 = 0;
            for (int j = 0; j < n_lines; ++j)
            {
                if (vp_association[j] != j)
                {
                    continue;
                }
                cv::Vec3d p1_homo{lines[j][0], lines[j][1], 1};
                auto p1_dot = vp_zen_line.dot(p1_homo);
                cv::Vec3d p2_homo{lines[j][2], lines[j][3], 1};
                auto p2_dot = vp_zen_line.dot(p2_homo);
                auto min_dot = lines_dirs[j].dot(lzmin_v);
                auto max_dot = lines_dirs[j].dot(lzmax_v);
                bool r;
                //dot(vp_zen_line, [lines(j,1); lines(j,2); 1])*(-1)^s > 0
                //&& dot(vp_zen_line, [lines(j,3); lines(j,4); 1])*(-1)^s > 0
                //&& abs(dot(lines_dir(j,:)',lzmin_V)) < cos(pi/8)
                //&& abs(dot(lines_dir(j,:)',lzmax_V)) < cos(pi/8)
                auto cos_pi_8 = std::cos(CV_PI / 8.);
                if ((s % 2) == 0)
                {
                    r = p1_dot < 0 and p2_dot < 0;
                }
                else
                {
                    r = p1_dot > 0 and p2_dot > 0;
                }
                r = r and (min_dot < cos_pi_8 and min_dot > -cos_pi_8) and (max_dot < cos_pi_8 and max_dot > cos_pi_8);
                if (r)
                {
                    ++n_lines_hvp;
                    auto a = -line_angle2(hvp, centroids[j]);
                    if (a < hvp_amin)
                    {
                        hvp_amin = a;
                        idmin2 = j;
                    }
                    else if (a > hvp_amax)
                    {
                        hvp_amax = a;
                        idmax2 = j;
                    }
                }
            }

            if (n_lines_hvp > n_lines_min)
            {
                Transform transform;
                transform.K = K;
                if (!K.empty())
                {
                    cv::Mat x = Ki * hvp_hom.t();
                    x /= cv::norm(x);
                    auto d1 = pp_zen_line.dot(hvp_hom); //dot(pp_zen_line, vp)
                    auto d2 = horizon_line.dot(vp_zen); //dot(horizon_line, vp_zen)

                    // TODO: bug here?
                    if ((d1 < 0 && d2 < 0) || (d1 > 0 && d2 > 0))
                    {
                        x = -x;
                    }
                    if ((d1 < 0 && d2 < 0 && s == 0) || (d1 > 0 && d2 < 0 && s == 1) || (d1 > 0 && d2 > 0 && s == 0) || (d1 < 0 && d2 > 0 && s == 1))
                    {
                        x = -x;
                    }
                    auto z = x.cross(y);
                    assert(x.rows == 3 and x.cols == 1);
                    assert(y.rows == 3 and y.cols == 1);
                    assert(z.rows == 3 and z.cols == 1);
                    transform.R = cv::Mat::eye(3, 3, CV_64F);
                    transform.R.at<double>(0, 0) = x.at<double>(0, 0);
                    transform.R.at<double>(0, 1) = y.at<double>(0, 0);
                    transform.R.at<double>(0, 2) = z.at<double>(0, 0);
                    transform.R.at<double>(1, 0) = x.at<double>(1, 0);
                    transform.R.at<double>(1, 1) = y.at<double>(1, 0);
                    transform.R.at<double>(1, 2) = z.at<double>(1, 0);
                    transform.R.at<double>(2, 0) = x.at<double>(2, 0);
                    transform.R.at<double>(2, 1) = y.at<double>(2, 0);
                    transform.R.at<double>(2, 2) = z.at<double>(2, 0);
                    transform.H = K * transform.R.inv() * Ki;
                } /* else {
                    auto lhmin = line_hmg_from_two_points(hvp, centroids[idmin2]);
                    auto lhmax = line_hmg_from_two_points(hvp, centroids[idmax2]);
                    std::vector<cv::Vec2d> c;
                    c.reserve(4);
                    c.emplace_back(line_hmg_intersect(lhmin, lzmin));
                    c.emplace_back(line_hmg_intersect(lhmax, lzmin));
                    c.emplace_back(line_hmg_intersect(lhmax, lzmax));
                    c.emplace_back(line_hmg_intersect(lhmin, lzmax));
                    std::sort(c.begin(), c.end(), [](const auto& a, const auto& b) { return a[0] < b[0]; });
                } */
                transforms.emplace_back(transform);
            }
        }
    }
}

#endif
