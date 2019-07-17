#ifndef _MNF_MODES_HPP_
#define _MNF_MODES_HPP_

#include <opencv2/opencv.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <omp.h>
//#include <mex.h>

#define dprint(expr) printf(#expr " = %g \n", (double)expr)

static inline double relative_entropy(double r, double p)
{
    return r * log(r / p) + (1 - r) * log((1 - r) / (1 - p));
}

static inline void mnf_intervals(double *histo, int L, double epsilon, int **intervals, int *Nintervals, double **H)
{
    double *hcum = (double *)malloc(sizeof(double) * L);
    int i, M, a, b;

    *intervals = (int *)malloc(sizeof(int) * L * L * 2);
    *H = (double *)malloc(sizeof(double) * L * L);
    *Nintervals = 0;

    hcum[0] = histo[0];
    for (i = 1; i < L; i++)
    {
        hcum[i] = hcum[i - 1] + histo[i];
    }
    M = hcum[L - 1];

#pragma omp parallel for schedule(dynamic) private(a, b) shared(intervals, H, Nintervals)
    for (a = 0; a < L; a++)
        for (b = a; b < L; b++)
        {
            double p = (double)(b - a + 1) / (double)L, r;
            int nsamples;

            if (a > 0)
                nsamples = hcum[b] - hcum[a - 1];
            else
                nsamples = hcum[b];
            r = (double)nsamples / (double)M;
            if (r > p)
            {
                double e = relative_entropy(r, p);
                if (e > log((double)L * (L + 1) / 2 / epsilon) / (double)M)
                {
#pragma omp critical
                    {
                        (*intervals)[2 * (*Nintervals) + 0] = a + 1;
                        (*intervals)[2 * (*Nintervals) + 1] = b + 1;
                        (*H)[*Nintervals] = e;
                        *Nintervals = *Nintervals + 1;
                    }
                }
            }
        }

    free(hcum);
}

static inline void mnf_gaps(double *histo, int L, double epsilon, int **gaps, int *Ngaps, double **H)
{
    double *hcum = (double *)malloc(sizeof(double) * L);
    int i, M, a, b;

    *gaps = (int *)malloc(sizeof(int) * L * L * 2);
    *H = (double *)malloc(sizeof(double) * L * L);
    *Ngaps = 0;

    hcum[0] = histo[0];
    for (i = 1; i < L; i++)
    {
        hcum[i] = hcum[i - 1] + histo[i];
    }
    M = hcum[L - 1];

#pragma omp parallel for schedule(dynamic) private(a, b) shared(gaps, H, Ngaps)
    for (a = 0; a < L; a++)
        for (b = a; b < L; b++)
        {
            double p = (double)(b - a + 1) / (double)L, r;
            int nsamples;

            if (a > 0)
                nsamples = hcum[b] - hcum[a - 1];
            else
                nsamples = hcum[b];
            r = (double)nsamples / (double)M;
            if (r < p)
            {
                double e = relative_entropy(r, p);
                if (e > log((double)L * (L + 1) / 2 / epsilon) / (double)M)
                {
#pragma omp critical
                    {
                        (*gaps)[2 * (*Ngaps) + 0] = a + 1;
                        (*gaps)[2 * (*Ngaps) + 1] = b + 1;
                        (*H)[*Ngaps] = e;
                        *Ngaps = *Ngaps + 1;
                    }
                }
            }
        }

    free(hcum);
}

static inline void mnf_modes(int *intervals, int Nintervals, double *Hintervals, int *gaps, int Ngaps, int **modes, int *Nmodes, double **Hmodes)
{
    int i;

    *modes = (int *)malloc(sizeof(int) * Nintervals * 2);
    *Hmodes = (double *)malloc(sizeof(double) * Nintervals);
    *Nmodes = 0;

#pragma omp parallel for schedule(dynamic) private(i) shared(modes, Hmodes, Nmodes)
    for (i = 0; i < Nintervals; i++)
    {
        int found = 0;
        int j = 0;
        while (j < Ngaps && found == 0)
        {
            if (gaps[2 * j + 0] >= intervals[2 * i + 0] && gaps[2 * j + 1] <= intervals[2 * i + 1])
                found = 1;
            j = j + 1;
        }
        if (found == 0)
        {
#pragma omp critical
            {
                (*modes)[2 * (*Nmodes) + 0] = intervals[2 * i + 0];
                (*modes)[2 * (*Nmodes) + 1] = intervals[2 * i + 1];
                (*Hmodes)[*Nmodes] = Hintervals[i];
                *Nmodes = *Nmodes + 1;
            }
        }
    }
}

static inline void max_mnf_modes(int *modes, int Nmodes, double *Hmodes, int **max_modes, int *Nmax_modes, double **Hmax_modes)
{
    int i;

    *max_modes = (int *)malloc(sizeof(int) * Nmodes * 2);
    *Hmax_modes = (double *)malloc(sizeof(double) * Nmodes);
    *Nmax_modes = 0;

#pragma omp parallel for schedule(dynamic) private(i) shared(max_modes, Hmax_modes, Nmax_modes)
    for (i = 0; i < Nmodes; i++)
    {
        int found = 0;
        int j = 0;
        while (j < Nmodes && found == 0)
        {
            if (j != i && modes[2 * j + 0] >= modes[2 * i + 0] && modes[2 * j + 1] <= modes[2 * i + 1] && Hmodes[j] > Hmodes[i])
                found = 1;
            j = j + 1;
        }
        if (found == 0)
        {
            int j = 0;
            while (j < Nmodes && found == 0)
            {
                if (j != i && modes[2 * i + 0] >= modes[2 * j + 0] && modes[2 * i + 1] <= modes[2 * j + 1] && Hmodes[j] > Hmodes[i])
                    found = 1;
                j = j + 1;
            }
            if (found == 0)
            {
#pragma omp critical
                {
                    (*max_modes)[2 * (*Nmax_modes) + 0] = modes[2 * i + 0];
                    (*max_modes)[2 * (*Nmax_modes) + 1] = modes[2 * i + 1];
                    (*Hmax_modes)[*Nmax_modes] = Hmodes[i];
                    *Nmax_modes = *Nmax_modes + 1;
                }
            }
        }
    }
}

//void mexFunction(int nlhs, mxArray *plhs[],
//                 int nrhs, const mxArray *prhs[])
void mnf_modes(const cv::Mat &histogram, const double &epsilon, std::vector<std::pair<int, int>> &m, std::vector<double> &e)
{
    /* 
     input: 1xN histogram; epsilon
     output: Noutx2 array of boundaries of the maximum meaningful modes, Noutx1 array of entropies inside the boundaries
  */

    double *outArray1, *outArray2;
    int i, N, Nout, Nintervals, Ngaps, Nmodes, Nmax_modes;
    //double epsilon;
    double *histo;
    int *intervals, *gaps, *modes, *max_modes;
    double *Hintervals, *Hgaps, *Hmodes, *Hmax_modes;

    ///* check that number of rows in first input argument is 1 */
    //if (mxGetM(prhs[0]) != 1)
    //{
    //    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector",
    //                      "First input must be a row vector.");
    //}

    ///* make sure the second input argument is type double */
    //if (!mxIsDouble(prhs[1]) ||
    //    mxIsComplex(prhs[1]))
    //{
    //    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble",
    //                      "Second input must be type double.");
    //}

    //N = mxGetN(prhs[0]);
    N = histogram.cols;
    //histo = mxGetPr(prhs[0]);
    histo = (double *)histogram.ptr();
    //epsilon = mxGetScalar(prhs[1]);

    mnf_intervals(histo, N, epsilon, &intervals, &Nintervals, &Hintervals);
    /* dprint(Nintervals);*/
    mnf_gaps(histo, N, epsilon, &gaps, &Ngaps, &Hgaps);
    /*dprint(Ngaps);*/
    mnf_modes(intervals, Nintervals, Hintervals, gaps, Ngaps, &modes, &Nmodes, &Hmodes);
    /*dprint(Nmodes);*/
    max_mnf_modes(modes, Nmodes, Hmodes, &max_modes, &Nmax_modes, &Hmax_modes);
    /*dprint(Nmax_modes);*/

    Nout = Nmax_modes;

    //plhs[0] = mxCreateDoubleMatrix(Nout, 2, mxREAL);
    //plhs[1] = mxCreateDoubleMatrix(Nout, 1, mxREAL);
    //outArray1 = mxGetPr(plhs[0]);
    //outArray2 = mxGetPr(plhs[1]);
    //for (i = 0; i < Nout; i++)
    //{
    //    outArray1[0 * Nout + i] = max_modes[2 * i + 0];
    //    outArray1[1 * Nout + i] = max_modes[2 * i + 1];
    //}
    //for (i = 0; i < Nout; i++)
    //{
    //    outArray2[i] = Hmax_modes[i];
    //}
    m.reserve(Nout);
    e.reserve(Nout);
    for (int i = 0; i < Nout; ++i)
    {
        m.emplace_back(std::pair<int, int>{max_modes[2 * i] - 1, max_modes[2 * i + 1] - 1});
        e.emplace_back(Hmax_modes[i]);
    }

    free(intervals);
    free(Hintervals);
    free(gaps);
    free(Hgaps);
    free(modes);
    free(Hmodes);
    free(max_modes);
    free(Hmax_modes);
}

static inline void mnf_modes_sort(std::vector<std::pair<int, int>> &m, std::vector<double> &e)
{
    assert(m.size() == e.size());
    if (e.size() == 0)
    {
        return;
    }
    std::vector<std::pair<double, int>> index;
    index.reserve(e.size());
    for (int i = 0; i < e.size(); ++i)
    {
        index.emplace_back(std::pair<double, int>{e[i], i});
    }

    std::sort(index.begin(), index.end(), [](const auto &a, const auto &b) { return a.first > b.first; });

    std::vector<std::pair<int, int>> _m;
    _m.reserve(e.size());
    for (int i = 0; i < e.size(); ++i)
    {
        e[i] = index[i].first;
        _m.emplace_back(m[index[i].second]);
    }
    m = _m;
}

#endif //_MNF_MODES_HPP_
