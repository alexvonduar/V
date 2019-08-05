#ifndef _ORTHORECTIFY_HPP_
#define _ORTHORECTIFY_HPP_

#include <opencv2/opencv.hpp>

#include "line_hmg_from_two_points.hpp"
#include "line_hmg_intersect.hpp"
#include "line_size.hpp"
#include "line_angle2.hpp"
#include "default_params.hpp"

typedef struct
{
    cv::Mat K;
    cv::Mat R;
    cv::Mat H;
} Transform;

template <typename T, std::size_t N>
static inline auto polygon_area(const std::array<cv::Point_<T>, N> &polygon) -> T
{

    double area = 0.0;
    //int num_points = polygon.size();
    for (int i = 0; i < N; ++i)
    {
        double x1 = polygon[i].x;
        double y1 = polygon[i].y;

        double x2 = polygon[(i + 1) % N].x;
        double y2 = polygon[(i + 1) % N].y;

        area += x1 * y2 - x2 * y1;
    }

    area = area < 0 ? -area : area;

    return area / 2;
}

void orthorectify(
    const cv::Mat &img,
    const std::vector<cv::Vec2d> &hvps,
    const std::vector<std::vector<int>> &hvp_groups,
    const cv::Vec2d &zenith,
    const std::vector<int> &z_group,
    const std::vector<cv::Vec4d> &lines,
    int n_lines_min,
    const cv::Mat &K,
    const cv::Vec3d &horizon_line,
    const Params &params,
    std::vector<Transform> &transforms)
{
    auto width = img.cols;
    auto height = img.rows;
    auto n_lines = lines.size();
    auto n_hvp = hvps.size();
    assert(n_hvp == hvp_groups.size());
    std::vector<int> vp_association(n_lines, -1);
    for (const auto &i : z_group)
    {
        vp_association[i] = n_lines;
    }
    for (int i = 0; i < hvp_groups.size(); ++i)
    {
        for (int j = 0; j < hvp_groups[i].size(); ++j)
        {
            vp_association[hvp_groups[i][j]] = i;
        }
    }

    cv::Mat Ki;
    cv::Vec3d pp_zen_line;
    cv::Vec3d vp_zen{zenith[0], zenith[1], 1};
    cv::Mat y;
    if (!K.empty())
    {
        Ki = K.inv();
        cv::Vec2d pp{K.at<double>(0, 2), K.at<double>(1, 2)};
        pp_zen_line = line_hmg_from_two_points(pp, zenith);

        y = Ki * vp_zen;
        y /= cv::norm(y);
        auto t = horizon_line.dot(vp_zen);
        if (t < 0)
        {
            y = -y;
        }
    }

    for (int i = 0; i < hvps.size(); ++i)
    {
        auto hvp = hvps[i];
        cv::Vec3d hvp_hom{hvp[0], hvp[1], 1};
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

        int c = 0;
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
            if (vp_association[c] == n_lines)
            { // z support lines;
                if (proj < 0)
                {
                    ++n_lines_zen[0];
                    if (proj > proj_max[0])
                    {
                        proj_max[0] = proj;
                        idmax[0] = c;
                    }
                    if (proj < proj_min[0])
                    {
                        proj_min[0] = proj;
                        idmin[0] = c;
                    }
                }
                else if (proj > 0)
                {
                    ++n_lines_zen[1];
                    if (proj > proj_max[1])
                    {
                        proj_max[1] = proj;
                        idmax[1] = c;
                    }
                    if (proj < proj_min[1])
                    {
                        proj_min[1] = proj;
                        idmin[1] = c;
                    }
                }
            }
            ++c;
        }

        for (int s = 0; s < 2; ++s)
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
                    if (vp_association[j] != i)
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
                    r = r and (min_dot < cos_pi_8 and min_dot > -cos_pi_8) and (max_dot < cos_pi_8 and max_dot > -cos_pi_8);
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
                        cv::Mat x = Ki * hvp_hom;
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
                        if (params.debug_fileid != nullptr)
                        {
                            fprintf(params.debug_fileid, "z [%f %f %f]\n", z.at<double>(0, 0), z.at<double>(1, 0), z.at<double>(2, 0));
                            fprintf(params.debug_fileid, "x [%f %f %f]\n", x.at<double>(0, 0), x.at<double>(1, 0), x.at<double>(2, 0));
                            fprintf(params.debug_fileid, "y [%f %f %f]\n", y.at<double>(0, 0), y.at<double>(1, 0), y.at<double>(2, 0));
                        }
                    }
                    else
                    {
                        auto lhmin = line_hmg_from_two_points(hvp, centroids[idmin2]);
                        auto lhmax = line_hmg_from_two_points(hvp, centroids[idmax2]);
                        std::vector<cv::Vec2d> c;
                        c.reserve(4);
                        c.emplace_back(line_hmg_intersect(lhmin, lzmin));
                        c.emplace_back(line_hmg_intersect(lhmax, lzmin));
                        c.emplace_back(line_hmg_intersect(lhmax, lzmax));
                        c.emplace_back(line_hmg_intersect(lhmin, lzmax));
                        std::array<int, 4> ids{0, 1, 2, 3};
                        std::sort(ids.begin(), ids.end(), [&c](const auto &a, const auto &b) { return c[a][0] < c[b][0]; });
                        auto idul = ids[1 - (c[ids[0]][1] < c[ids[1]][1])];
                        auto idll = ids[1 - (c[ids[0]][1] > c[ids[1]][1])];
                        auto up = (idll == ((idul % 4) + 1));
                        auto idc = idul;
                        std::array<cv::Point2d, 4> corners;
                        for (int j = 0; j < 4; ++j)
                        {
                            corners[j] = c[idc];
                            if (up)
                            {
                                ++idc;
                                idc %= 4;
                            }
                            else
                            {
                                if (idc)
                                {
                                    --idc;
                                }
                                else
                                {
                                    idc = 3;
                                }
                            }
                        }
                        auto M = cv::sum(corners);
                        M /= 4;
                        auto A = polygon_area<double, 4>(corners);
                        auto d1 = corners[0] - corners[1];
                        auto d2 = corners[2] - corners[3];
                        auto d3 = corners[1] - corners[2];
                        auto d4 = corners[3] - corners[0];
                        auto s1 = cv::norm(d1) + cv::norm(d2);
                        auto s2 = cv::norm(d3) + cv::norm(d4);
                        auto AR = s1 / s2;
                        auto w = std::sqrt(A / AR);
                        auto h = AR * w;
                        std::array<cv::Point2d, 4> pointsR{
                            cv::Point2d{M[0] - w / 2, M[1] - h / 2},
                            cv::Point2d{M[0] - w / 2, M[1] + h / 2},
                            cv::Point2d{M[0] + w / 2, M[1] + h / 2},
                            cv::Point2d{M[0] + w / 2, M[1] - h / 2}};
                        transform.H = cv::findHomography(corners, pointsR, cv::noArray(), cv::LMEDS);
                        if (params.debug_fileid != nullptr)
                        {
                            fprintf(params.debug_fileid, "M [%f %f]\n", M[0], M[1]);
                            fprintf(params.debug_fileid, "A %f\n", A);
                            fprintf(params.debug_fileid, "AR %f\n", AR);
                            fprintf(params.debug_fileid, "w %f h %f\n", w, h);
                            for (int k = 0; k < pointsR.size(); ++k)
                            {
                                fprintf(params.debug_fileid, "[%f %f] -> [%f %f]\n", corners[k].x, corners[k].y, pointsR[k].x, pointsR[k].y);
                            }
                        }
                    }
                    transforms.emplace_back(transform);
                }
            }
        }
    }
}

#endif
