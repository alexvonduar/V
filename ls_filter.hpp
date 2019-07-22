#ifndef _LS_FILTER_HPP_
#define _LS_FILTER_HPP_

#include <opencv2/opencv.hpp>

#include "util/line_size.hpp"

//function seglines = ls_filter(thres_aligned, length_t, lsd)
void ls_filter(const double &thres_aligned,
               const double &length_t,
               const std::vector<cv::Vec4d> &lsd,
               std::vector<cv::Vec4d> &seglines)
{
    auto num_lines = lsd.size();
    std::vector<cv::Vec4d> _seglines = lsd;
    std::vector<double> L0;
    L0.reserve(num_lines);
    std::vector<cv::Vec2d> seglinesC;
    seglinesC.reserve(num_lines);
    std::vector<cv::Vec2d> seglinesJ;
    seglinesJ.reserve(num_lines);
    //L0 = line_size(seglines);
    //seglinesC(1:size(seglines,1),1) = (seglines(:,3)+seglines(:,1))/2;
    //seglinesC(1:size(seglines,1),2) = (seglines(:,4)+seglines(:,2))/2;
    //seglinesJ(1:size(seglines,1),1) = -(seglines(:,4)-seglines(:,2))./L0;
    //seglinesJ(1:size(seglines,1),2) = (seglines(:,3)-seglines(:,1))./L0;
    for (const auto &ls : _seglines)
    {
        auto len = line_size(ls);
        seglinesC.emplace_back(cv::Vec2d{(ls[2] + ls[0]) / 2, (ls[3] + ls[1]) / 2});
        seglinesJ.emplace_back(cv::Vec2d{(ls[1] - ls[3]) / len, (ls[2] - ls[0]) / len});
        L0.emplace_back(len);
    }
    //[L0, I] = sort(L0, 'descend');
    //seglines =  seglines(I,:);
    //seglinesC = seglinesC(I,:);
    //seglinesJ = seglinesJ(I,:);
    {
        std::vector<std::pair<double, int>> index;
        index.reserve(L0.size());
        for (int i = 0; i < L0.size(); ++i)
        {
            index.emplace_back(std::pair<double, int>{L0[i], i});
        }
        std::sort(index.begin(), index.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
        std::vector<cv::Vec4d> seglines_sorted;
        seglines_sorted.reserve(num_lines);
        std::vector<cv::Vec2d> seglinesC_sorted;
        seglinesC_sorted.reserve(num_lines);
        std::vector<cv::Vec2d> seglinesJ_sorted;
        seglinesJ_sorted.reserve(num_lines);
        for (int i = 0; i < index.size(); ++i)
        {
            L0[i] = index[i].first;
            auto id = index[i].second;
            seglines_sorted.emplace_back(_seglines[id]);
            seglinesC_sorted.emplace_back(seglinesC[id]);
            seglinesJ_sorted.emplace_back(seglinesJ[id]);
        }
        _seglines = seglines_sorted;
        seglinesC = seglinesC_sorted;
        seglinesJ = seglinesJ_sorted;
    }
    //i = 1;
    //while i < size(seglines,1)
    //    C = seglinesC(i,1:2);
    //    J = seglinesJ(i,:);
    //    X1 = abs((seglines((i+1):end,1:2)-repmat(C,size(seglines,1)-i,1))*J');
    //    X2 = abs((seglines((i+1):end,3:4)-repmat(C,size(seglines,1)-i,1))*J');
    //    I = intersect(find(X1 < thres_aligned),find(X2 < thres_aligned))+i;
    //    seglines(I,:) = [];
    //    L0(I) = [];
    //    seglinesC(I,:) = [];
    //    seglinesJ(I,:) = [];
    //    i = i+1;
    //end
    //int total_count = 0;
    for (int i = 0; i < _seglines.size(); ++i)
    {
        if (L0[i] == 0)
        {
            continue;
        }
        auto C = seglinesC[i];
        auto J = seglinesJ[i];
        //int count = 0;
        for (int j = i + 1; j < _seglines.size(); ++j)
        {
            if (L0[j] == 0)
            {
                continue;
            }
            cv::Vec2d p1{_seglines[j][0], _seglines[j][1]};
            cv::Vec2d p2{_seglines[j][2], _seglines[j][3]};
            p1[0] -= C[0];
            p1[1] -= C[1];
            p2[0] -= C[0];
            p2[1] -= C[1];
            auto X1 = std::abs(p1[0] * J[0] + p1[1] * J[1]);
            auto X2 = std::abs(p2[0] * J[0] + p2[1] * J[1]);
            if (X1 < thres_aligned and X2 < thres_aligned)
            {
                L0[j] = 0;
                //++count;
            }
        }
        //std::cout << "iter " << i << " removed " << count << " lines" << std::endl;
        //total_count += count;
    }
    //std::cout << "removed " << total_count << " lines" << std::endl;
    //i = 1;
    //while i <= size(seglines,1)
    //    if L0(i) < length_t/2
    //        seglines(i,:) = [];
    //        L0(i) = [];
    //        seglinesJ(i,:) = [];
    //        i = i-1;
    //    end
    //    i=i+1;
    //end
    for (int i = 0; i < _seglines.size(); ++i)
    {
        if (L0[i] >= (length_t / 2))
        {
            seglines.emplace_back(_seglines[i]);
        }
    }
    //end
}

#endif //_LS_FILTER_HPP_
