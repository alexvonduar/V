#ifndef _ZL_PREDICT_HPP_
#define _ZL_PREDICT_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"
#include "util/line_angle2.hpp"
#include "util/line_hmg_from_two_points.hpp"
#include "util/mnf_modes.hpp"
#include "util/histcounts.hpp"

/*
function zl = zl_predict(ls, dist_max, u0, v0, width, height, params);
hist = [];
ang = [];
for i = 1:size(ls,1)
    ang(i) = line_angle2(ls(i,:));
    if ang(i) < -pi/4
        ang(i) = ang(i) + pi;
    end
    if abs(ang(i)-pi/2) < params.theta_v
        l = line_hmg_from_two_points(ls(i,1:2),ls(i,3:4));
        dist = abs(dot(l, [width/2;height/2;1]));
        if dist < dist_max
            hist(end+1) = ang(i);
        end
    end
end
[N,edges0] = histcounts(hist,(pi/2-pi/8):pi/180:(pi/2+pi/8));
N(find(N <= 5)) = 0;
if sum(N) == 0
    N(round(length(N)/2)) = 1;
end
edges = (edges0(1:end-1)+edges0(2:end))/2;
[max_modes, H] = mnf_modes(N, 1);
if isempty(max_modes)
    max_modes = [1 size(N,2)];
    H = 0;
else
    [~,I] = sort(H, 'descend');
    H = H(I);
    max_modes = max_modes(I,:);
end
zl = [];
for i = 1:size(max_modes,1)
    Ni = zeros(1,size(N,2));
    a = max_modes(i,1);
    b = max_modes(i,2);
    Ni(a:b) = N(a:b);
    m = max(Ni);
    I = find(Ni == m);
    for j = I(1)
        a = edges(j);
        l = abs(v0/sin(a));
        zl(end+1) = u0+l*cos(pi-a);
    end
end
end
*/

//function zl = zl_predict(ls, dist_max, u0, v0, width, height, params);
static inline void zl_predict(const std::vector<cv::Vec4d> ls, const double &dist_max, const double &u0, const double &v0, const int &width, const int &height, const Params &params, std::vector<double> &zl)
{
    auto num_lines = ls.size();
    std::vector<double> hist;
    hist.reserve(num_lines);
    std::vector<double> ang;
    ang.reserve(num_lines);
    //hist = [];
    //ang = [];
    //for i = 1:size(ls,1)
    //    ang(i) = line_angle2(ls(i,:));
    //    if ang(i) < -pi/4
    //        ang(i) = ang(i) + pi;
    //    end
    //    if abs(ang(i)-pi/2) < params.theta_v
    //        l = line_hmg_from_two_points(ls(i,1:2),ls(i,3:4));
    //        dist = abs(dot(l, [width/2;height/2;1]));
    //        if dist < dist_max
    //            hist(end+1) = ang(i);
    //        end
    //    end
    //end
    for (const auto &l : ls)
    {
        auto angle = line_angle2(l);
        if (angle < -CV_PI / 4)
        {
            angle += CV_PI;
        }
        if (std::abs(angle - CV_PI / 2) < params.theta_v)
        {
            auto lm = line_hmg_from_two_points(l);
            auto dist = std::abs(lm.dot(cv::Vec3d{double(width) / 2, double(height) / 2, 1.}));
            if (dist < dist_max)
            {
                hist.emplace_back(angle);
            }
        }
        ang.emplace_back(angle);
    }
    //[N,edges0] = histcounts(hist,(pi/2-pi/8):pi/180:(pi/2+pi/8));
    cv::Mat N;
    std::vector<double> edges0;
    float range[] = {CV_PI * 3 / 8, CV_PI * 5 / 8};
    //--cv::calcHist(hist, std::vector<int>{0}, cv::noArray(), N, std::vector<int>{45}, std::vector<float>{CV_PI * 3 / 8, CV_PI * 5 / 8}, false);
    histcounts(hist, 45, range, N, edges0);
    //N(find(N <= 5)) = 0;
    cv::threshold(N, N, 5, 0, cv::THRESH_TOZERO);
    assert(N.rows == 1 and N.cols == 45);
    //if sum(N) == 0
    //    N(round(length(N)/2)) = 1;
    //end
    auto sum = cv::sum(N);
    //assert(sum.channels == 1);
    if (sum[0] == 0)
    {
        // TODO: fix zero case here
        *(N.ptr(0, N.cols / 2)) = 1;
    }

    //edges = (edges0(1:end-1)+edges0(2:end))/2;
    std::vector<double> edges;
    edges.reserve(N.cols);
    for (int i = 0; i < N.cols; ++i)
    {
        edges.emplace_back(CV_PI * 3 / 8 + CV_PI * (2 * i + 1) / 360);
    }

    //[max_modes, H] = mnf_modes(N, 1);
    std::vector<std::pair<int, int>> max_modes;
    std::vector<double> H;
    mnf_modes(N, 1, max_modes, H);

    //if isempty(max_modes)
    //    max_modes = [1 size(N,2)];
    //    H = 0;
    //else
    //    [~,I] = sort(H, 'descend');
    //    H = H(I);
    //    max_modes = max_modes(I,:);
    //end
    if (max_modes.size() == 0)
    {
        max_modes.emplace_back(std::pair<int, int>{0, N.cols});
        H.emplace_back(0);
    }
    else
    {
        std::vector<std::pair<double, int>> sH;
        sH.reserve(H.size());
        for (int i = 0; i < H.size(); ++i)
        {
            sH.emplace_back(std::pair<double, int>{H[i], i});
        }
        std::sort(sH.begin(), sH.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
        std::vector<std::pair<int, int>> s_modes;
        s_modes.reserve(sH.size());
        for (int i = 0; i < sH.size(); ++i)
        {
            H[i] = sH[i].first;
            s_modes.emplace_back(max_modes[sH[i].second]);
        }
        max_modes = s_modes;
    }

    for (const auto &mode : max_modes)
    {
        auto a = mode.first;
        auto b = mode.second;
        auto max = N.at<double>(0, a);
        auto max_id = a;
        for (auto i = a + 1; i < b; ++i)
        {
            auto _max = N.at<double>(0, i);
            if (_max > max)
            {
                max = _max;
                max_id = i;
            }
        }

        auto edge = edges[max_id];
        auto l = std::abs(v0 / std::sin(edge));
        auto z = u0 + l * std::cos(CV_PI - edge);
        zl.emplace_back(z);
        //end
    }
    //end
}

#endif //_ZL_PREDICT_HPP_
