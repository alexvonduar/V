#ifndef _HL_SAMPLE_HPP_
#define _HL_SAMPLE_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"
#include "util/line_hmg_from_two_points.hpp"
#include "util/normalize.hpp"
#include "util/atand.hpp"

//function [samp_homo, samp_left, samp_right] = hl_sample(zenith_homo, modes_homo, modes_offset, modes_left, modes_right, H, u0, v0, width, height, focal, params)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void hl_sample(
    const cv::Vec<T, 3> &zenith_homo,
    const std::vector<cv::Vec<T, 3>> &modes_homo,
    const std::vector<T> &modes_offset,
    const std::vector<T> &modes_left,
    const std::vector<T> &modes_right,
    const std::vector<T> &H,
    const T &u0,
    const T &v0,
    const int &width,
    const int &height,
    const T &focal,
    const Params &params,
    std::vector<cv::Vec<T, 3>> &samp_homo,
    std::vector<T> &samp_left,
    std::vector<T> &samp_right)
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
#if 1
            if (H[i] <= 0)
            {
                //draw = rand()*uniformrnd(end);
                //id = min(find(uniformrnd >= draw));
                //orand = x(id)*height + modes_offset(i);
                auto id = rng.uniform(0, 100000);
                orand = T(4.) * id / 100000. - 2. + modes_offset[i];
                //std::cout << j << ": " << id << " " << orand << std::endl;
            }
            else
            {
                //orand = normrnd(0,params.sigma*height) + modes_offset(i);
                auto id = rng.gaussian(params.sigma * height);
                orand = id + modes_offset[i];
                //std::cout << j << ": " << id << " " << orand << std::endl;
                //end
            }
#else
            orand = -1.9 + (3.8 * j) / (nsamp - 1) + modes_offset[i];
#endif
            //mnf_center = [u0+orand*cos(-tilt+pi/2) v0+orand*sin(-tilt+pi/2)];
            cv::Vec<T, 2> mnf_center{u0 + orand * std::cos(-tilt + CV_PI / 2), v0 + orand * std::sin(-tilt + CV_PI / 2)};
            //hmnf = line_hmg_from_two_points(mnf_center, mnf_center+[cos(-tilt) sin(-tilt)]);
            auto hmnf = line_hmg_from_two_points(mnf_center, cv::Vec<T, 2>{std::cos(-tilt) + mnf_center[0], std::sin(-tilt) + mnf_center[1]});
            if (0 and params.debug_fileid != nullptr)
            {
                fprintf(params.debug_fileid, "%d: %.1079g [%.1079g %.1079g] -> [%.1079g %.1079g %.1079g]\n", j, orand, mnf_center[0], mnf_center[1], hmnf[0], hmnf[1], hmnf[2]);
            }
            //samp_left(end+1) = -hmnf(3)/hmnf(2);
            auto left = -hmnf[2] / hmnf[1];
            samp_left.emplace_back(left);
            //samp_right(end+1) = (-hmnf(1)*width-hmnf(3))/hmnf(2);
            auto right = (-hmnf[0] * width - hmnf[2]) / hmnf[1];
            samp_right.emplace_back(right);
            //mode_seg = [0, samp_left(end), width, samp_right(end)];
            cv::Vec<T, 4> mode_seg{0, left, (T)width, right};
            //samp_homo(1:3,end+1) = normalize(mode_seg, width, height, focal);
            cv::Vec<T, 3> nline;
            normalize(mode_seg, width, height, focal, nline);
            samp_homo.emplace_back(nline);
            //end
        }
        //end
    }
    //end
}

#endif //_HL_SAMPLE_HPP_
