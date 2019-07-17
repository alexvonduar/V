#ifndef _LINES_NORMAL_HPP_
#define _LINES_NORMAL_HPP_

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>

//function vp_homo = lines_normal(line_homo, b)
// this function solve the homogeneous least squares prob
// A = line_homo', find x s.t. min(||Ax||)
// if b is given, then add a constraint: b'x = 0, which means x is forced to
// be orthogonal to b. [Zhai et al. 2016]
static inline void lines_normal(const std::vector<cv::Vec3d> &line_homo, const cv::Mat &b, cv::Vec3d &vp_homo)
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
        auto A = l * lt;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU);
        U = svd.matrixU();

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
        vp_homo = cv::Vec3d{vph(0, 0), vph(1, 0), vph(2, 0)};
        vp_homo /= cv::norm(vp_homo);

        //end
    }

    // force to the z-positive semisphere
    //vp_homo = vp_homo * sign(vp_homo(3)+eps);
    if ((vp_homo[2] + std::nextafter(0., 1.)) < 0)
    {
        vp_homo *= -1;
    }

    // [0 0 0] -> [0 1 0]
    //if sum(vp_homo) == 0
    if (vp_homo[0] + vp_homo[1] + vp_homo[2] == 0)
    {
        vp_homo = cv::Vec3d{0., 1., 0.};
        //end
    }
}

static inline void lines_normal(const std::vector<cv::Vec3d> &line_homo, const cv::Vec3d &b, cv::Vec3d &vp_homo)
{
    cv::Mat _b(3, 1, CV_64FC1);
    _b.at<double>(0, 0) = b[0];
    _b.at<double>(1, 0) = b[1];
    _b.at<double>(2, 0) = b[2];
    lines_normal(line_homo, _b, vp_homo);
}

#endif //_LINES_NORMAL_HPP_
