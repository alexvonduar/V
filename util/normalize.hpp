#ifndef _NORMALIZE_HPP_
#define _NORMALIZE_HPP_

#include <opencv2/opencv.hpp>

/*
% convert lines or points from cartesian system (image frame) to homogenous system
% according to [Zhai et al., CVPR'2016]
function geometry_homo = normalize(geometry_img, width, height, focal)
center = [width, height] / 2;
if size(geometry_img, 2) == 2
  geometry_normalized = bsxfun(@minus, geometry_img, center) / focal;
else
  geometry_normalized = bsxfun(@minus, geometry_img, [center, center]); % / focal;
end
N = size(geometry_normalized, 1);
if size(geometry_normalized, 2) == 2 
  geometry_homo = [geometry_normalized ones(N, 1)]';
%   geometry_homo = bsxfun(@times, geometry_homo, sign(geometry_homo(3,:) + eps));
else  
  x1 = normalize(geometry_img(:,[1,2]), width, height, focal);
  x2 = normalize(geometry_img(:,[3,4]), width, height, focal);
  geometry_homo = cross(x1, x2);
  geometry_homo = bsxfun(@times, geometry_homo, sign(geometry_homo(2,:) + eps));
end
geometry_homo = bsxfun(@rdivide, geometry_homo, sqrt(sum(geometry_homo.^2, 1)));
*/

static inline void normalize(const cv::Vec2d &point, const int &width, const int &height, const double &focal, cv::Vec3d &npoint)
{
  //cv::Vec2d center{double(width) / 2, double(height) / 2};
  //auto geometry_normalized = point - center;
  //geometry_normalized /= focal;
  npoint = cv::Vec3d{(point[0] - ((double)width / 2.)) / focal, (point[1] - ((double)height / 2.)) / focal, 1};
  auto norm = std::sqrt(npoint[0] * npoint[0] + npoint[1] * npoint[1] + 1);
  //npoint /= norm;
  npoint[0] /= norm;
  npoint[1] /= norm;
  npoint[2] /= norm;
}

static inline void normalize(const cv::Vec4d &line, const int &width, const int &height, const double &focal, cv::Vec3d &nline)
{
  //cv::Vec4d center{double(width) / 2, double(height) / 2, double(width) / 2, double(height) / 2};
  //auto geometry_normalized = line - center;
  //geometry_normalized /= focal;
  cv::Vec3d x1, x2;
  normalize(cv::Vec2d{line[0], line[1]}, width, height, focal, x1);
  normalize(cv::Vec2d{line[2], line[3]}, width, height, focal, x2);
  auto geometry_homo = x1.cross(x2);
  if ((geometry_homo[1] + std::nextafter(0.0, 1.0)) < 0)
  {
    //geometry_homo *= -1;
    geometry_homo[0] *= -1;
    geometry_homo[1] *= -1;
    geometry_homo[2] *= -1;
  }
  auto norm = std::sqrt(geometry_homo[0] * geometry_homo[0] + geometry_homo[1] * geometry_homo[1] + geometry_homo[2] * geometry_homo[2]);
  //geometry_homo /= norm;
  geometry_homo[0] /= norm;
  geometry_homo[1] /= norm;
  geometry_homo[2] /= norm;
  nline = geometry_homo;
}

#endif //_NORMALIZE_HPP_
