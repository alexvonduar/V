#ifndef _LINES_NORMAL_HPP_
#define _LINES_NORMAL_HPP_

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>

#include "default_params.hpp"

//function vp_homo = lines_normal(line_homo, b)
// this function solve the homogeneous least squares prob
// A = line_homo', find x s.t. min(||Ax||)
// if b is given, then add a constraint: b'x = 0, which means x is forced to
// be orthogonal to b. [Zhai et al. 2016]
static inline void lines_normal(const std::vector<cv::Vec3d> &line_homo, const cv::Mat &b, const Params &params, cv::Vec3d &vp_homo)
{
    Eigen::MatrixXd l(3, line_homo.size());
    for (int i = 0; i < line_homo.size(); ++i)
    {
        l(0, i) = line_homo[i][0];
        l(1, i) = line_homo[i][1];
        l(2, i) = line_homo[i][2];
    }
    auto lt = l.transpose();
    Eigen::MatrixXd U;
    //if ~exist('b', 'var')
    if (b.empty())
    {

        //[U, ~, ~] = svd(line_homo*line_homo');
#if 1
        auto A = l * lt;
#else // compute l*l' by hand to verify Eigen matrix multiplication
        Eigen::MatrixXd A(3, 3);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                auto N = line_homo.size();
                long double s = 0;
                for (int k = 0; k < N; ++k)
                {
                    s += line_homo[k][i] * line_homo[k][j];
                }
                A(i, j) = s;
            }
        }
#endif
        if (0 and params.debug_fileid != nullptr)
        {
            fprintf(params.debug_fileid, "lines SVD: \n");
            fprintf(params.debug_fileid, "l: \n");
            for (int i = 0; i < l.rows(); ++i)
            {
                for (int j = 0; j < l.cols(); ++j)
                {
                    fprintf(params.debug_fileid, "%.1079g ", l(i, j));
                }
                fprintf(params.debug_fileid, "\n");
            }
            fprintf(params.debug_fileid, "lt: \n");
            for (int i = 0; i < lt.rows(); ++i)
            {
                for (int j = 0; j < lt.cols(); ++j)
                {
                    fprintf(params.debug_fileid, "%.1079g ", lt(i, j));
                }
                fprintf(params.debug_fileid, "\n");
            }
            fprintf(params.debug_fileid, "A: \n");
            for (int i = 0; i < A.rows(); ++i)
            {
                for (int j = 0; j < A.cols(); ++j)
                {
                    fprintf(params.debug_fileid, "%.13g ", A(i, j));
                }
                fprintf(params.debug_fileid, "\n");
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU);
        U = svd.matrixU();

        if (params.debug_fileid != nullptr)
        {
            auto sign = 1;
            if (U(2, 2) < 0)
            {
                sign = -1;
            }
            fprintf(params.debug_fileid, "U: \n");
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    fprintf(params.debug_fileid, "%.13g ", U(i, j) * sign);
                }
                fprintf(params.debug_fileid, "\n");
            }
        }

        //vp_homo = U(:,3);
        vp_homo = cv::Vec3d{U(0, 2), U(1, 2), U(2, 2)};
    }
    else
    {
        Eigen::MatrixXd B(b.rows, b.cols);
        for (int i = 0; i < b.rows; ++i)
        {
            for (int j = 0; j < b.cols; ++j)
            {
                B(i, j) = b.at<double>(i, j);
            }
        }
        //[U, ~, ~] = svd(b);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeFullU);
        //p = U(:,2); q = U(:,3);
        Eigen::MatrixXd pq(3, 2);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                pq(i, j) = svd.matrixU()(i, j + 1);
            }
        }

        //Am = line_homo'*[p,q];
        auto Am = lt * pq;
        auto Amt = Am.transpose();

        //[U, ~, ~] = svd(Am'*Am);
        Eigen::JacobiSVD<Eigen::MatrixXd> svdAm(Amt * Am, Eigen::ComputeFullU);
        //lambda = U(:,end);
        Eigen::MatrixXd lambda(2, 1);
        lambda << svdAm.matrixU()(0, 1), svdAm.matrixU()(1, 1);
        //vp_homo = [p,q]*lambda;
        auto vph = pq * lambda;
        auto x = vph(0, 0);
        auto y = vph(1, 0);
        auto z = vph(2, 0);
        //vp_homo = cv::Vec3d{x, y, z};
        auto norm = std::sqrt(x * x + y * y + z * z);
        //vp_homo /= cv::norm(vp_homo);
        vp_homo = cv::Vec3d{x / norm, y / norm, z / norm};

        //end
    }

    // force to the z-positive semisphere
    //vp_homo = vp_homo * sign(vp_homo(3)+eps);
    if ((vp_homo[2] + std::nextafter(0., 1.)) < 0)
    {
        //vp_homo *= -1;
        vp_homo[0] *= -1;
        vp_homo[1] *= -1;
        vp_homo[2] *= -1;
    }

    // [0 0 0] -> [0 1 0]
    //if sum(vp_homo) == 0
    if (vp_homo[0] + vp_homo[1] + vp_homo[2] == 0)
    {
        vp_homo = cv::Vec3d{0., 1., 0.};
        //end
    }
}

static inline void lines_normal(const std::vector<cv::Vec3d> &line_homo, const cv::Vec3d &b, const Params &params, cv::Vec3d &vp_homo)
{
    cv::Mat _b(3, 1, CV_64FC1);
    _b.at<double>(0, 0) = b[0];
    _b.at<double>(1, 0) = b[1];
    _b.at<double>(2, 0) = b[2];
    lines_normal(line_homo, _b, params, vp_homo);
}

#endif //_LINES_NORMAL_HPP_
