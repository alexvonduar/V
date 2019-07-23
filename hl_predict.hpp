#ifndef _HL_PREDICT_HPP_
#define _HL_PREDICT_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"
#include "util/mnf_modes.hpp"
#include "util/normalize.hpp"
#include "util/line_hmg_from_two_points.hpp"
#include "util/atand.hpp"
#include "util/histcounts.hpp"

//function [modes_homo, modes_offset, modes_left, modes_right, H] = hl_predict(seglines, zenith_homo, u0, v0, width, height, focal, params)
static inline void hl_predict(const std::vector<cv::Vec4d> &seglines,
                              const cv::Vec3d &zenith_homo,
                              const double &u0,
                              const double &v0,
                              const int &width,
                              const int &height,
                              const double &focal,
                              const Params &params,
                              std::vector<cv::Vec3d> &modes_homo,
                              std::vector<double> &modes_offset,
                              std::vector<double> &modes_left,
                              std::vector<double> &modes_right,
                              std::vector<double> &H)
{
    auto L = params.L_h;
    auto x2 = zenith_homo[0] / zenith_homo[2];
    auto y2 = zenith_homo[1] / zenith_homo[2];
    auto tilt = atan(x2, y2);
    if (std::isnan(tilt))
    {
        tilt = 0;
        //end
    }
    //offsets = [];
    std::vector<double> offsets;
    offsets.reserve(std::max(seglines.size(), size_t(height)));
    //for i = 1:size(seglines,1)
    for (const auto &line : seglines)
    {
        //v = [seglines(i,3)-seglines(i,1) seglines(i,4)-seglines(i,2)];
        auto v = cv::Vec2d{line[2] - line[0], line[3] - line[1]};
        //scal = dot(v/norm(v),[cos(-tilt), sin(-tilt)]);
        v /= cv::norm(v);
        auto scal = v.dot(cv::Vec2d{std::cos(-tilt), std::sin(-tilt)});
        //ang = acos(scal);
        auto ang = std::acos(scal);
        //if abs(ang) < params.theta_h*pi/180
        if (std::abs(ang) < (params.theta_h * CV_PI / 180.))
        {
            //offsets(end+1) = dot([(seglines(i,3)+seglines(i,1))/2-u0 (seglines(i,4)+seglines(i,2))/2-v0],[cos(-tilt+pi/2) sin(-tilt+pi/2)]);
            auto v1 = cv::Vec2d{(line[2] + line[0]) / 2 - u0, (line[3] + line[1]) / 2 - v0};
            auto v2 = cv::Vec2d{std::cos(-tilt + CV_PI / 2), std::sin(-tilt + CV_PI / 2)};
            offsets.emplace_back(v1.dot(v2));
            //end
        }
        //end
    }
    //if isempty(offsets)
    //    offsets = -height/2:height/2;
    //end
    if (offsets.size() == 0)
    {
        for (int i = 0; i < height; ++i)
        {
            offsets.emplace_back(-(double(height)) / 2. + i);
        }
    }

    //[N,edges] = histcounts(offsets,L);
    cv::Mat N;
    std::vector<double> edges;
    histcounts(offsets, L, nullptr, N, edges);

    if (params.debug_fileid != nullptr)
    {
        auto nhist = N.cols;
        fprintf(params.debug_fileid, "hist %d: ", nhist);
        for (int i = 0; i < nhist; ++i) {
            fprintf(params.debug_fileid, "%f ", N.at<double>(0, i));
        }
        fprintf(params.debug_fileid, "\n");
        auto noff = offsets.size();
        fprintf(params.debug_fileid, "offsets %ld: ", noff);
        for (int i = 0; i < noff; ++i)
        {
            fprintf(params.debug_fileid, "%f ", offsets[i]);
        }
        fprintf(params.debug_fileid, "\n");
        auto n = edges.size();
        fprintf(params.debug_fileid, "edges %ld: ", n);
        for (int i = 0; i < n; ++i)
        {
            fprintf(params.debug_fileid, "%f ", edges[i]);
        }
        fprintf(params.debug_fileid, "\n");
    }

    //[max_modes, H] = mnf_modes(N, 1);
    std::vector<std::pair<int, int>> max_modes;
    //std::vector<double> H;
    mnf_modes(N, 1, max_modes, H);
    if (max_modes.size() == 0)
    {
        //max_modes(1,:) = [1 size(N,2)];
        max_modes.emplace_back(std::pair<int, int>{0, N.cols});
        //H(1) = -1;
        H = std::vector<double>{-1};
    }
    else
    {
        //[~,I] = sort(H, 'descend');
        //H = H(I);
        //max_modes = max_modes(I,:);
        std::vector<std::pair<double, int>> index;
        index.reserve(H.size());
        for (int i = 0; i < H.size(); ++i)
        {
            index.emplace_back(std::pair<double, int>{H[i], i});
        }
        std::sort(index.begin(), index.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
        std::vector<std::pair<int, int>> modes;
        modes.reserve(H.size());
        for (int i = 0; i < H.size(); ++i)
        {
            H[i] = index[i].first;
            modes.emplace_back(max_modes[index[i].second]);
        }
        max_modes = modes;
    }

    //nmodes = size(max_modes,1);
    auto nmodes = max_modes.size();
    //modes_offset = [];
    //std::vector<double> modes_offset;
    modes_offset.reserve(max_modes.size());
    //modes_left = [];
    //std::vector<double> modes_left;
    modes_left.reserve(max_modes.size());
    //modes_right = [];
    //std::vector<double> modes_right;
    modes_right.reserve(max_modes.size());
    //modes_homo = [];
    //std::vector<cv::Vec3d> modes_homo;
    modes_homo.reserve(max_modes.size());
    //for i = 1:nmodes
    for (const auto &mode : max_modes)
    {
        //Ni = zeros(1,size(N,2));
        //a = max_modes(i,1);
        //b = max_modes(i,2);
        //Ni(a:b) = N(a:b);
        //[~,bin] = max(Ni);
        double max = N.at<double>(0, mode.first);
        int bin = mode.first;
        for (int i = bin + 1; i <= mode.second; ++i)
        {
            if (N.at<double>(0, i) > max)
            {
                max = N.at<double>(0, i);
                bin = i;
            }
        }
        //modes_offset(end+1) = edges(1)+bin/L*(edges(L)-edges(1));
        //auto edge = *minmax.first + (*minmax.second - *minmax.first) * bin / L;
        auto a = edges[0];
        auto b = edges[L - 1];
        auto edge = a + (bin + 1) * (b - a) / L;
        modes_offset.emplace_back(edge);
        //mnf_center = [u0+modes_offset(end)*cos(-tilt+pi/2) v0+modes_offset(end)*sin(-tilt+pi/2)];
        auto mnf_center = cv::Vec2d{u0 + edge * std::cos(-tilt + CV_PI / 2), v0 + edge * std::sin(-tilt + CV_PI / 2)};
        //hmnf = line_hmg_from_two_points(mnf_center, mnf_center+[cos(-tilt) sin(-tilt)]);
        auto tilt_v = cv::Vec2d{std::cos(-tilt), std::sin(-tilt)};
        auto hmnf = line_hmg_from_two_points(mnf_center, mnf_center + tilt_v);
        //modes_left(end+1) = -hmnf(3)/hmnf(2);
        modes_left.emplace_back(-hmnf[2] / hmnf[1]);
        //modes_right(end+1) = (-hmnf(1)*width-hmnf(3))/hmnf(2);
        modes_right.emplace_back((-hmnf[0] * width - hmnf[2]) / hmnf[1]);
        //mode_seg = [0, modes_left(1), width, modes_right(1)];
        cv::Vec4d mode_seg{0., modes_left[0], double(width), modes_right[0]};
        //modes_homo(1:3,end+1) = normalize(mode_seg, width, height, focal);
        cv::Vec3d nline;
        normalize(mode_seg, width, height, focal, nline);
        modes_homo.emplace_back(nline);
        //end
    }
    //end
}

#endif //_HL_PREDICT_HPP_
