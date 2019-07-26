#ifndef _DRAW_HPP_
#define _DRAW_HPP_

#include <opencv2/opencv.hpp>

typedef struct
{
    cv::Mat img;
    int x;
    int y;
    int winW;
    int winH;
    std::string winName;
} DrawData;

template <typename T>
static inline void draw(const std::vector<cv::Vec<T, 4>> &lines, const cv::Scalar *pcolor, const T &scale, cv::Mat &img)
{
    cv::Scalar color;
    if (pcolor == nullptr)
    {
        color = cv::Scalar{255, 0, 0};
    }
    else
    {
        color = *pcolor;
    }

    for (const auto &line : lines)
    {
        auto p1 = cv::Point_<T>{line[0], line[1]};
        auto p2 = cv::Point_<T>{line[2], line[3]};
        if (scale != 1)
        {
            p1 *= scale;
            p2 *= scale;
        }
        cv::line(img, p1, p2, color);
    }
}

void on_hscroll(int pos, void *userdata)
{
    DrawData *p = (DrawData *)userdata;
    p->x = pos;
    if (pos + p->winW > p->img.cols)
    {
        p->x = p->img.cols - p->winW;
    }
    cv::Mat winImage = p->img(cv::Rect(p->x, p->y, p->winW, p->winH));
    cv::imshow(p->winName, winImage);
}

void on_vscroll(int pos, void *userdata)
{
    DrawData *p = (DrawData *)userdata;
    p->y = pos;
    if (pos + p->winH > p->img.rows)
    {
        p->y = p->img.rows - p->winH;
    }
    cv::Mat winImage = p->img(cv::Rect(p->x, p->y, p->winW, p->winH));
    cv::imshow(p->winName, winImage);
}

static inline void draw(const std::string &win_name, const cv::Mat &img, const std::vector<cv::Vec4d> &lines, const cv::Scalar& _color, const bool &do_scale = false)
{
    DrawData data;
    data.winName = win_name;
    img.copyTo(data.img);
    data.winH = 900;
    data.winW = 600;
    cv::Scalar color = _color;
    if (do_scale)
    {
        double scale = 1;
        if (data.winH < data.img.rows)
        {
            scale = (double)(900) / data.img.rows;
        }
        cv::resize(data.img, data.img, cv::Size(), scale, scale, cv::INTER_AREA);
        draw(lines, &color, scale, data.img);
        cv::imshow(win_name, data.img);
    }
    else
    {
        draw(lines, &color, 1., data.img);
        cv::namedWindow(data.winName, cv::WINDOW_AUTOSIZE);

        if (data.winH >= data.img.rows)
            data.winH = data.img.rows - 1;
        if (data.winW >= data.img.cols)
            data.winW = data.img.cols - 1;

        int scrolHight = 0;
        data.y = scrolHight;
        int scrolWidth = 0;
        data.x = scrolWidth;
        cv::createTrackbar("Vscroll", data.winName, &scrolHight, (data.img.rows - data.winH), on_vscroll, &data);
        cv::createTrackbar("Hscroll", data.winName, &scrolWidth, (data.img.cols - data.winW), on_hscroll, &data);
        on_vscroll(0, &data);
    }
    while (cv::waitKey(0) != 'q')
    {

    } //while
    cv::destroyWindow(data.winName);
}

#endif //_DRAW_HPP_
