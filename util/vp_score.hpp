#ifndef _VP_SCORE_HPP_
#define _VP_SCORE_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"

//function [score, horgroup] = vp_score(vp_homo, lines_homo, score_function)
static inline double vp_score(const cv::Vec3d &vp_homo, const std::vector<cv::Vec3d> &lines_homo,
                              const Params &params,
                              //double (*score_function)(const double &, const double &),
                              //double& score,
                              std::vector<int> &horgroup)
{

    //cos_mat = vp_homo'*lines_homo;
    //theta_mat = abs(asind(cos_mat));
    //score_mat = score_function(theta_mat);
    //tmp = mat2cell(score_mat, ones(size(score_mat,1),1), size(score_mat,2));
    //horgroup = cellfun(@find, tmp, 'uniformoutput', false);
    //score = cellfun(@sum, tmp);
    ///std::vector<double> cos_mat;
    ///cos_mat.reserve(lines_homo.size());
    ///std::vector<double> theta_mat;
    ///theta_mat.reserve(lines_homo.size());
    ///std::vector<double> score_mat;
    ///score_mat.reserve(lines_homo.size());
    double score = 0.;
    for (int i = 0; i < lines_homo.size(); ++i)
    {
        auto cos = vp_homo.dot(lines_homo[i]);
        auto theta = std::abs(std::asin(cos) * 180 / CV_PI);
        //theta_mat.emplace_back(theta);
        //cos_mat.emplace_back(cos);
        //score_mat.emplace_back(score_function(params.theta_con, theta));
        auto _score = params.score_function(params.theta_con, theta);
        if (_score > 0)
        {
            horgroup.emplace_back(i);
            score += _score;
        }
    }

    return score;
}


#endif
