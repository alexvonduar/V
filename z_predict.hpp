#ifndef _Z_PREDICT_HPP_
#define _Z_PREDICT_HPP_

#include <opencv2/opencv.hpp>

#include "default_params.hpp"
#include "lines_normal.hpp"
#include "vp_ransac_refinement.hpp"
#include "util/atand.hpp"

#if 0
//function [score, horgroup] = vp_score(vp_homo, lines_homo, score_function)
static inline double vp_score(const cv::Vec<T, 3> &vp_homo, const std::vector<cv::Vec<T, 3>> &lines_homo,
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
        if (_score)
        {
            horgroup.emplace_back(i);
            score += _score;
        }
    }

    return score;
}
#else
#include "util/vp_score.hpp"
#endif

//function [zenith_homo, zengroupIds] = z_predict(lines_homo, zenith_line_homo, params, refine)
template <typename T, typename Dummy = typename std::enable_if<std::is_floating_point<T>::value>::type>
static inline void z_predict(const std::vector<cv::Vec<T, 3>> &lines_homo,
                             const cv::Vec<T, 3> &zenith_line_homo,
                             const Params &params,
                             const bool &refine,
                             cv::Vec<T, 3> &zenith_homo,
                             std::vector<int> &zengroupIds)
{
    zengroupIds.clear();
    // threshold zenith LSs ids

    //lines_tilt_ortho = atand(-lines_homo(2,:) ./ lines_homo(1,:));
    //zenith_tilt = atand(-zenith_line_homo(2,1) ./ zenith_line_homo(1,1));
    //verticalInd = find(abs(lines_tilt_ortho - zenith_tilt) < params.theta_z);
    //zenith_homo_pred = lines_normal(lines_homo(:,verticalInd));
    auto zenith_tilt = atand(-zenith_line_homo[1], zenith_line_homo[0]);
    std::vector<cv::Vec<T, 3>> lines;
    lines.reserve(lines_homo.size());
    std::vector<int> validIds;
    validIds.reserve(lines_homo.size());
    int count = 0;
    for (const auto &l : lines_homo)
    {
        auto lines_tilt_ortho = atand(-l[1], l[0]);
        //std::cout << lines_tilt_ortho << std::endl;
        if (std::abs(lines_tilt_ortho - zenith_tilt) < params.theta_z)
        {
            //std::cout << "line: " << count << " select" << std::endl;
            lines.emplace_back(l);
            validIds.emplace_back(count);
        }
        ++count;
    }

    cv::Vec<T, 3> zenith_homo_pred;
    lines_normal(lines, cv::Mat_<T>(), params, zenith_homo_pred);
    if (params.debug_fileid != nullptr)
    {
        fprintf(params.debug_fileid, "z_predict: zenith tilt %f\n", zenith_tilt);
        fprintf(params.debug_fileid, "z_predict: %ld vertical lines\n", lines.size());
        for (int i = 0; i < validIds.size(); ++i)
        {
            fprintf(params.debug_fileid, "%d ", validIds[i]);
        }
        fprintf(params.debug_fileid, "\n");
    }

    // refinement

    if (refine)
    {
        //[zenith_homo, ~] = vp_ransac_refinement(lines_homo(:,verticalInd), params);
        std::vector<int> inlierId;
        vp_ransac_refinement(lines, params, zenith_homo, inlierId);
    }
    else
    {
        zenith_homo = zenith_homo_pred;
        //end
    }

    // LS grouping

    //[~, groups] = vp_score(zenith_homo, lines_homo, params.score_function);
    //zengroupIds = groups{1};
    auto score = vp_score(zenith_homo, lines_homo, params, zengroupIds);
}

#endif //_Z_PREDICT_HPP_
