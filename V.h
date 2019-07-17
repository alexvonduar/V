#ifndef _V_H_
#define _V_H_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"

int V(const Params &params,
      const cv::Mat &img,
      std::vector<cv::Vec2d> &hl,
      std::vector<cv::Vec2d> &hvps,
      std::vector<std::vector<int>> &hvp_groups,
      cv::Vec2d &z,
      std::vector<int> &z_group,
      std::vector<cv::Vec4d> &ls);

#endif //_V_H_
