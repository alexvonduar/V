#include <opencv2/opencv.hpp>

#include "V.h"

int main(int argc, char *argv[])
{
    Params params;
    default_params(params);
    params.include_infinite_hvps = 1; // includes the detection of infinite
    // horizontal VPs. TO BE REMOVED if one want to get the same results as
    // in our ECCV'2018 paper

    cv::Mat img = cv::imread("/home/alex/work/V/test_images/260001.jpg");

    std::vector<cv::Vec2d> hl;
    std::vector<cv::Vec2d> hvps;
    std::vector<std::vector<int>> hvp_groups;
    cv::Vec2d z;
    std::vector<int> z_groups;
    std::vector<cv::Vec4d> ls;
    V(params, img, hl, hvps, hvp_groups, z, z_groups, ls);
}