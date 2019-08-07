#include <cstdio>

#include <opencv2/opencv.hpp>

#include "V.h"

int main(int argc, char *argv[])
{
    std::string fname("/home/alex/work/V/test_images/260001.jpg");
    if (argc > 1)
    {
        fname = std::string(argv[1]);
    }

    Params params;
    default_params(params);
    params.include_infinite_hvps = 1; // includes the detection of infinite
    // horizontal VPs. TO BE REMOVED if one want to get the same results as
    // in our ECCV'2018 paper

    if (1)
    {
        auto last_p = fname.find_last_of('.');
        auto debug_fname = fname.substr(0, last_p);
        debug_fname += "_cpp_debug.txt";
        std::cout << "debug output: " << debug_fname << std::endl;
        params.debug_fileid = fopen(debug_fname.c_str(), "w");
    }

    cv::Mat img = cv::imread(fname);

    std::vector<cv::Vec2d> hl;
    std::vector<cv::Vec2d> hvps;
    std::vector<std::vector<int>> hvp_groups;
    cv::Vec2d z;
    std::vector<int> z_groups;
    std::vector<cv::Vec4d> ls;
    vp(params, img, hl, hvps, hvp_groups, z, z_groups, ls);

    if (params.debug_fileid != nullptr)
    {
        fclose(params.debug_fileid);
    }
}