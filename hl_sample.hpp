#ifndef _HL_SAMPLE_HPP_
#define _HL_SAMPLE_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"
#include "util/line_hmg_from_two_points.hpp"
#include "util/normalize.hpp"
#include "util/atand.hpp"

//function [samp_homo, samp_left, samp_right] = hl_sample(zenith_homo, modes_homo, modes_offset, modes_left, modes_right, H, u0, v0, width, height, focal, params)

static inline void hl_sample(
    const cv::Vec3d &zenith_homo,
    const std::vector<cv::Vec3d> &modes_homo,
    const std::vector<double> &modes_offset,
    const std::vector<double> &modes_left,
    const std::vector<double> &modes_right,
    const std::vector<double> &H,
    const double &u0,
    const double &v0,
    const int &width,
    const int &height,
    const double &focal,
    const Params &params,
    std::vector<cv::Vec3d> &samp_homo,
    std::vector<double> &samp_left,
    std::vector<double> &samp_right)
{
    //rng(1) % fix random seed
    cv::RNG rng(std::chrono::system_clock::now().time_since_epoch().count());
    samp_homo = modes_homo;
    samp_left = modes_left;
    samp_right = modes_right;
    auto S = params.S;
    //tilt = atan((zenith_homo(1)/zenith_homo(3))/(zenith_homo(2)/zenith_homo(3)));
    auto tilt = atan(zenith_homo[0] / zenith_homo[2], zenith_homo[1] / zenith_homo[2]);
    if (std::isnan(tilt))
    {
        tilt = 0;
    }

    //nsamp = ceil(S/size(modes_homo,2));
    auto num_modes = modes_homo.size();
    auto nsamp = (S + num_modes - 1) / num_modes;
    //x = linspace(-2,2,1e5);
    //uniformrnd = cumsum(ones(size(x)));
    //for i = 1:size(modes_homo,2)
    for (int i = 0; i < num_modes; ++i)
    {
        //for j = 1:(nsamp-1)
        for (int j = 0; j < (nsamp - 1); ++j)
        {
            double orand = 0;
            if (H[i] <= 0)
            {
                //draw = rand()*uniformrnd(end);
                //id = min(find(uniformrnd >= draw));
                //orand = x(id)*height + modes_offset(i);
                auto id = rng.uniform(0, 100000);
                auto orand = double(4.) * id / 100000. - 2.;
            }
            else
            {
                //orand = normrnd(0,params.sigma*height) + modes_offset(i);
                auto orand = rng.gaussian(params.sigma * height) + modes_offset[i];
                //end
            }
            //mnf_center = [u0+orand*cos(-tilt+pi/2) v0+orand*sin(-tilt+pi/2)];
            cv::Vec2d mnf_center{u0 + orand * std::cos(-tilt + CV_PI / 2), v0 + orand * std::sin(-tilt + CV_PI / 2)};
            //hmnf = line_hmg_from_two_points(mnf_center, mnf_center+[cos(-tilt) sin(-tilt)]);
            auto hmnf = line_hmg_from_two_points(mnf_center, mnf_center + cv::Vec2d{std::cos(-tilt), std::sin(-tilt)});
            //samp_left(end+1) = -hmnf(3)/hmnf(2);
            auto left = -hmnf[2] / hmnf[1];
            samp_left.emplace_back(left);
            //samp_right(end+1) = (-hmnf(1)*width-hmnf(3))/hmnf(2);
            auto right = (-hmnf[0] * width - hmnf[2]) / hmnf[1];
            samp_right.emplace_back(right);
            //mode_seg = [0, samp_left(end), width, samp_right(end)];
            cv::Vec4d mode_seg{0, left, (double)width, right};
            //samp_homo(1:3,end+1) = normalize(mode_seg, width, height, focal);
            cv::Vec3d nline;
            normalize(mode_seg, width, height, focal, nline);
            samp_homo.emplace_back(nline);
            //end
        }
        //end
    }
    //end
}

#endif //_HL_SAMPLE_HPP_
