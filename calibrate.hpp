/*-------------------------------------------------------------------------
*
* calibrate.m  
* 
* Copyright 2018 Gilles Simon <gsimon@loria.fr> Université de Lorraine
*
* Version 1.0 - September 5, 2018
*
* calibrate.m allows computing the camera focal length and finding a 
* Manhattan frame in a set of VPs. This function is explained in our 
* previous paper, section 4.3:
*
*   "Gilles Simon, Antoine Fond, Marie-Odile Berger.
*   A Simple and Effective Method to Detect Orthogonal Vanishing Points in
*   Uncalibrated Images of Man-Made Environments.
*   Eurographics 2016, May 2016, Lisbon, Portugal." 
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Affero General Public License for more details.
* 
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------------*/

#ifndef _CALIBRATE_HPP_
#define _CALIBRATE_HPP_

#include <opencv2/opencv.hpp>

typedef struct
{
    double x;
    double y;
    double ff;
    double vvpc;
    int v;
    int v2;
} MHT;

//function [focal, manh_vps, confident] = calibrate(zvp, hvp, width, height)
int calibrate(
    const cv::Vec2d &zvp,
    const std::vector<cv::Vec2d> &hvp,
    const int &width,
    const int &height,
    double &focal,
    std::vector<cv::Vec2d> &manh_vps)
{
    int infty = 4;
    double accuracy = 2;
    //nt = 2^(5+min(accuracy, log2(width)-5));
    auto nt = std::pow(2., 5 + std::min(accuracy, std::log2(width) - 5));
    auto step_h = (double)width / nt;
    auto Dt_h = std::atan(step_h / width);
    auto u0 = (double)width / 2;
    auto v0 = (double)height / 2;
    //manh_vps = [];
    focal = -1;
    int confident = -1;
    //mht = [];
    std::vector<MHT> mht;
    //mht_count = 0;
    int mht_best = 0;
    //hvp_count = size(hvp,1);
    //if hvp_count > 1
    //for v = 1:hvp_count
    for (int v = 0; v < hvp.size(); ++v)
    {
        //for v2 = (v+1):hvp_count
        for (int v2 = v + 1; v2 < hvp.size(); ++v2)
        {
            //x1 = hvp(v,1) - u0;
            auto x1 = hvp[v][0] - u0;
            //y1 = hvp(v,2) - v0;
            auto y1 = hvp[v][1] - v0;
            //x2 = hvp(v2,1) - u0;
            auto x2 = hvp[v2][0] - u0;
            //y2 = hvp(v2,2) - v0;
            auto y2 = hvp[v2][1] - v0;
            //ff = sqrt(-(x1*x2+y1*y2));
            // TODO: fix sqrt here
            auto ff = -(x1 * x2 + y1 * y2);
            //if ff*ff > 0 && ff > 0.28*width && ff < 3.8*width
            if (ff > 0)
            {
                ff = std::sqrt(ff);
                if (ff > 0.28 * width && ff < 3.8 * width)
                {
                    //mht_count = mht_count + 1;
                    //vv1 = [x1/ff;y1/ff;1];
                    cv::Point3d vv1{x1 / ff, y1 / ff, 1.};
                    //vv2 = [x2/ff;y2/ff;1];
                    cv::Point3d vv2{x2 / ff, y2 / ff, 1.};
                    //vv3 = cross(vv1,vv2);
                    // TODO: fix cross product here
                    auto vv3 = vv1.cross(vv2);
                    //mht(mht_count, 1) = ff*vv3(1)/vv3(3)+u0;
                    //mht(mht_count, 2) = ff*vv3(2)/vv3(3)+v0;
                    auto mht2 = ff * vv3.y / vv3.z + v0;
                    //mht(mht_count, 3) = ff;
                    //mht(mht_count, 5) = v;
                    //mht(mht_count, 6) = v2;
                    //at_vz1 = atan(mht(mht_count, 2)/double(height))/Dt_h;
                    auto at_vz1 = std::atan(mht2 / height) / Dt_h;
                    //at_vz2 = atan(zvp(2)/double(height))/Dt_h;
                    auto at_vz2 = std::atan(zvp[1] / height) / Dt_h;
                    //mht(mht_count, 4) = min(abs(at_vz1-at_vz2), abs(at_vz1+at_vz2));
                    mht.emplace_back(MHT{ff * vv3.x / vv3.z + u0,
                                         ff * vv3.y / vv3.z + v0,
                                         ff,
                                         std::min(std::abs(at_vz1 - at_vz2), std::abs(at_vz1 + at_vz2)),
                                         v,
                                         v2});
                }
                //end
            }
            //end
        }
        //end
    }
    //end
    //if mht_count > 0
    if (mht.size())
    {
        //[vvpc_min, mht_best] = min(mht(:,4));
        auto it_mht_best = std::min_element(mht.begin(), mht.end(), [](const auto &a, const auto &b) { return a.vvpc < b.vvpc; });
        mht_best = it_mht_best - mht.begin();
        auto vvpc_min = it_mht_best->vvpc;
        //if vvpc_min < 2 || abs(zvp(2)-v0) >= infty*height
        if (vvpc_min < 2 or std::abs(zvp[1] - v0) >= infty * height)
        {
            if (vvpc_min < 2)
            {
                confident = 3;
            }
            else
            {
                //[~, mht_best] = max(abs(mht(:,2)-v0));
                mht_best = 0;
                auto mht_max = std::abs(mht[0].y - v0);
                for (int i = 1; i < mht.size(); ++i)
                {
                    auto _mht_max = std::abs(mht[i].y - v0);
                    if (_mht_max > mht_max)
                    {
                        mht_best = i;
                        mht_max = _mht_max;
                    }
                }
                confident = 1;
                //end
            }
            //id1 = mht(mht_best, 5);
            auto id1 = mht[mht_best].v;
            //id2 = mht(mht_best, 6);
            auto id2 = mht[mht_best].v2;
            //x1 = hvp(id1,1);
            auto x1 = hvp[id1][0];
            //x2 = hvp(id2,1);
            auto x2 = hvp[id2][0];
            //y = mht(mht_best, 2);
            auto y = mht[mht_best].y;
            //if (y > v0 && x1 > x2) || (y < v0 && x1 < x2)
            if ((y > v0 && x1 > x2) || (y < v0 && x1 < x2))
            {
                //mht(mht_best, 5) = id2;
                mht[mht_best].v = id2;
                //mht(mht_best, 6) = id1;
                mht[mht_best].v2 = id1;
                //end
            }
            //focal = mht(mht_best,3);
            focal = mht[mht_best].ff;
            //manh_vps = [hvp(mht(mht_best, 5),1) hvp(mht(mht_best, 5),2);hvp(mht(mht_best, 6),1) hvp(mht(mht_best, 6),2);zvp(1) zvp(2)];
            manh_vps.emplace_back(hvp[mht[mht_best].v]);
            manh_vps.emplace_back(hvp[mht[mht_best].v2]);
            manh_vps.emplace_back(zvp);
            //end
        }
        //end
    }
    //if confident == -1 && hvp_count > 0 && abs(zvp(2) - v0) < infty*height
    if (confident == -1 and hvp.size() and (std::abs(zvp[1] - v0) < infty * height))
    {
        //id_finite = find(abs(hvp(:,1)-u0) < infty*width);
        auto id_finite = std::find_if(hvp.begin(), hvp.end(), [&u0, &infty, &width](const auto &a) { return (a[1] - u0) < (infty * width); });
        //if length(id_finite) > 0
        if (id_finite != hvp.end())
        {
            //y1 = hvp(id_finite(1),2) - v0;
            auto y1 = (*id_finite)[1] - v0;
            //y2 = zvp(2) - v0;
            auto y2 = zvp[1] - v0;
            //ff = sqrt(-y1*y2);
            auto ff = -y1 * y2;
            //if ff*ff > 0 && ff > 0.28*width && ff < 3.8*width
            if (ff > 0)
            {
                ff = std::sqrt(ff);
                if (ff > 0.28 * width and ff < 3.8 * width)
                {
                    focal = ff;
                    //K = [[focal 0 u0];[0 focal v0];[0 0 1]];
                    cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
                    K.at<double>(0, 0) = focal;
                    K.at<double>(0, 2) = u0;
                    K.at<double>(1, 1) = focal;
                    K.at<double>(1, 2) = v0;
                    //r2 = [hvp(id_finite(1),1)-u0; hvp(id_finite(1),2)-v0; focal];
                    cv::Vec3d r2{(*id_finite)[0] - u0, (*id_finite)[1] - v0, focal};
                    //r2 = r2/norm(r2);
                    r2 /= cv::norm(r2);
                    //r3 = [zvp(1)-u0; zvp(2)-v0; ff];
                    cv::Vec3d r3{zvp[0] - u0, zvp[1] - v0, ff};
                    //r3 = r3/norm(r3);
                    r3 /= cv::norm(r3);
                    //vpx = K*cross(r2,r3);
                    cv::Vec3d cross(r2.cross(r3));
                    cv::Mat vpx = K * cross;
                    //vpx = vpx/vpx(3);
                    auto vp_x = vpx.at<double>(0, 0) / vpx.at<double>(0, 2);
                    auto vp_y = vpx.at<double>(0, 1) / vpx.at<double>(0, 2);
                    //manh_vps = [vpx(1)+u0 vpx(2)+v0; hvp(id_finite(1),1) hvp(id_finite(1),2); zvp(1) zvp(2)];
                    manh_vps.emplace_back(cv::Vec2d{vp_x + u0, vp_y + v0});
                    manh_vps.emplace_back(*id_finite);
                    manh_vps.emplace_back(zvp);
                    confident = 2;
                }
                //end
            }
            //end
        }
        //end
    }
    //end
    return confident;
}

#endif //_CALIBRATE_HPP_
