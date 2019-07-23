// The MIT License (MIT)
//
// Copyright (c) 2016 Menghua Zhai
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef _RANSAC_INTERSECTION_HPP_
#define _RANSAC_INTERSECTION_HPP_

#include <opencv2/opencv.hpp>

typedef struct
{
    int iterNum;        //50;
    double thInlrRatio; // .02;
    double thDist;      // = sind(opt.theta_con);
} RANSAC_Parameter;

//function [M, index] = randIndex(lines_homo, maxIndex,len)
void randIndex(const std::vector<cv::Vec3d> &lines_homo, const int &maxIndex, const int &len, cv::Vec3d &M, std::vector<int> &index)
{
    //INDEX = RANDINDEX(MAXINDEX,LEN)
    //   randomly, non-repeatedly select LEN integers from 1:MAXINDEX
    cv::RNG rng(std::chrono::system_clock::now().time_since_epoch().count());
    //std::vector<int> index;
    index.resize(len);

    //index = randperm(maxIndex, len);
    rng.fill(index, cv::RNG::UNIFORM, 0, maxIndex);

    // normal of great circle consisting of 2 randomly sampled points.
    //M = cross(lines_homo(:, index(1)), lines_homo(:, index(2)));
    M = lines_homo[index[0]].cross(lines_homo[index[1]]);
    //M = M / norm(M);
    M /= cv::norm(M);
    //M = M * sign(M(2) + eps); // stay in y-positive plane
    if (M[1] < 0)
    {
        M *= -1;
    }
}

//function d = calcdist(M, X)
double calcdist(const cv::Vec3d &M, const cv::Vec3d &X)
{
    //d = abs(M'*X);
    return std::abs(M.dot(X));
}

//function [M, inlierIdx] = ransac_intersection(lines_homo, ransacCoef)
void ransac_intersection(const std::vector<cv::Vec3d> &lines_homo, const RANSAC_Parameter &ransacCoef, cv::Vec3d &M, std::vector<int> &inlierIdx)
{
    //rng(1);

    int minPtNum = 2; // sample 2 LSs to find the VP
    auto iterNum = ransacCoef.iterNum;
    auto thInlrRatio = ransacCoef.thInlrRatio;
    auto thDist = ransacCoef.thDist;
    //ptNum = size(lines_homo,2);
    auto ptNum = lines_homo.size();
    //thInlr = round(thInlrRatio*ptNum);
    auto thInlr = std::round(thInlrRatio * ptNum);

    //if ptNum < minPtNum
    if (ptNum < minPtNum)
    {
        //M = zeros(3,2);
        M.zeros();
        inlierIdx.clear();
        //N = 0;
        return;
        //end
    }

    //inlrNum = zeros(1,iterNum);
    std::vector<int> inlrNum; //(iterNum, 0);
    inlrNum.reserve(iterNum);
    //fLib = cell(1,iterNum);
    std::vector<cv::Vec3d> fLib;
    fLib.reserve(iterNum);

    //for p = 1:iterNum
    for (int p = 0; p < iterNum; )//++p)
    {
        // 1. fit using random points
        //M = randIndex(lines_homo, ptNum, minPtNum);
        std::vector<int> index;
        randIndex(lines_homo, ptNum, minPtNum, M, index);

        // 2. count the inliers, if more than thInlr, refit; else iterate
        //dist = calcdist(M, lines_homo);

        //inlier1 = find(dist < thDist);
        std::vector<int> inlier1;
        inlier1.reserve(lines_homo.size());
        for (int j = 0; j < lines_homo.size(); ++j)
        {
            auto dist = calcdist(M, lines_homo[j]);
            if (dist < thDist)
            {
                inlier1.emplace_back(j);
            }
        }

        //inlrNum(p) = length(inlier1);
        //inlrNum[p] = inlier1.size();
        //if length(inlier1) < thInlr, continue; end
        if (inlier1.size() < thInlr)
        {
            continue;
        }
        else
        {
            ++p;
            inlrNum.emplace_back(inlier1.size());
            //fLib{p} = M;
            fLib.emplace_back(M);
            //end
        }
    }

    // 3. choose the coef with the most inliers
    //[~,idx] = max(inlrNum);
    auto iter = std::max_element(inlrNum.begin(), inlrNum.end());
    auto idx = iter - inlrNum.begin();
    //if ~isempty(fLib{idx})
    if (fLib[idx][2] != 0)
    {
        //M = fLib{idx};
        M = fLib[idx];
        //dist = calcdist(M, lines_homo);
        //inlierIdx = find(dist < thDist);
        for (int j = 0; j < lines_homo.size(); ++j)
        {
            auto dist = calcdist(M, lines_homo[j]);
            if (dist < thDist)
            {
                inlierIdx.emplace_back(j);
            }
        }
    }
    else
    {
        printf("no enough inliers.\n");
        // randomly generate M
        //M = randIndex(lines_homo, ptNum, minPtNum);
        std::vector<int> index;
        randIndex(lines_homo, ptNum, minPtNum, M, index);
        //inlierIdx = [];
        inlierIdx.clear();
        //end
    }
}

#endif //_RANSAC_INTERSECTION_HPP_
