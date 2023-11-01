//
// Created by Kai Zhao on 4/20/20.
//

#ifndef SZ_STATISTIC_HPP
#define SZ_STATISTIC_HPP

#include "Config.hpp"

namespace SZ3 {
    template<class T>
    T data_range(const T *data, size_t num) {
        T max = data[0];
        T min = data[0];
        for (size_t i = 1; i < num; i++) {
            if (max < data[i]) max = data[i];
            if (min > data[i]) min = data[i];
        }
        return max - min;
    }

    int factorial(int n) {
        return (n == 0) || (n == 1) ? 1 : n * factorial(n - 1);
    }

    double computeABSErrBoundFromPSNR(double psnr, double threshold, double value_range) {
        double v1 = psnr + 10 * log10(1 - 2.0 / 3.0 * threshold);
        double v2 = v1 / (-20);
        double v3 = pow(10, v2);
        return value_range * v3;
    }

    template<class T>
    void calAbsErrorBound(Config &conf, const T *data, T range = 0) {
        if (conf.errorBoundMode != EB_ABS) {
            if (conf.errorBoundMode == EB_REL) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = conf.relErrorBound * ((range > 0) ? range : data_range(data, conf.num));
            } else if (conf.errorBoundMode == EB_PSNR) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = computeABSErrBoundFromPSNR(conf.psnrErrorBound, 0.99, ((range > 0) ? range : data_range(data, conf.num)));
            } else if (conf.errorBoundMode == EB_L2NORM) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = sqrt(3.0 / conf.num) * conf.l2normErrorBound;
            } else if (conf.errorBoundMode == EB_ABS_AND_REL) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = std::min(conf.absErrorBound, conf.relErrorBound * ((range > 0) ? range : data_range(data, conf.num)));
            } else if (conf.errorBoundMode == EB_ABS_OR_REL) {
                conf.errorBoundMode = EB_ABS;
                conf.absErrorBound = std::max(conf.absErrorBound, conf.relErrorBound * ((range > 0) ? range : data_range(data, conf.num)));
            } else {
                printf("Error, error bound mode not supported\n");
                exit(0);
            }
        }
    }

    template<typename Type>
    double autocorrelation1DLag1(const Type *data, size_t numOfElem, Type avg) {
        double cov = 0;
        for (size_t i = 0; i < numOfElem; i++) {
            cov += (data[i] - avg) * (data[i] - avg);
        }
        cov = cov / numOfElem;

        if (cov == 0) {
            return 0;
        } else {
            int delta = 1;
            double sum = 0;

            for (size_t i = 0; i < numOfElem - delta; i++) {
                sum += (data[i] - avg) * (data[i + delta] - avg);
            }
            return sum / (numOfElem - delta) / cov;
        }
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements, double &psnr, double &nrmse, double &max_diff) {
        size_t i = 0;
        double Max = ori_data[0];
        double Min = ori_data[0];
        max_diff = fabs(data[0] - ori_data[0]);
        double diff_sum = 0;
        double maxpw_relerr = 0;
        double sum1 = 0, sum2 = 0, l2sum = 0;
        for (i = 0; i < num_elements; i++) {
            sum1 += ori_data[i];
            sum2 += data[i];
            l2sum += data[i] * data[i];
        }
        double mean1 = sum1 / num_elements;
        double mean2 = sum2 / num_elements;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0, relerr = 0;

        double *diff = (double *) malloc(num_elements * sizeof(double));

        for (i = 0; i < num_elements; i++) {
            diff[i] = data[i] - ori_data[i];
            diff_sum += data[i] - ori_data[i];
            if (Max < ori_data[i]) Max = ori_data[i];
            if (Min > ori_data[i]) Min = ori_data[i];
            double err = fabs(data[i] - ori_data[i]);
            if (ori_data[i] != 0) {
                relerr = err / fabs(ori_data[i]);
                if (maxpw_relerr < relerr)
                    maxpw_relerr = relerr;
            }

            if (max_diff < err)
                max_diff = err;
            prodSum += (ori_data[i] - mean1) * (data[i] - mean2);
            sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
            sum4 += (data[i] - mean2) * (data[i] - mean2);
            sum += err * err;
        }
        double std1 = sqrt(sum3 / num_elements);
        double std2 = sqrt(sum4 / num_elements);
        double ee = prodSum / num_elements;
        double acEff = ee / std1 / std2;

        double mse = sum / num_elements;
        double range = Max - Min;
        psnr = 20 * log10(range) - 10 * log10(mse);
        nrmse = sqrt(mse) / range;

        double normErr = sqrt(sum);
        double normErr_norm = normErr / sqrt(l2sum);

        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.2G\n", max_diff);
        printf("Max relative error = %.2G\n", max_diff / (Max - Min));
        printf("Max pw relative error = %.2G\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.10G\n", psnr, nrmse);
        printf("normError = %f, normErr_norm = %f\n", normErr, normErr_norm);
        printf("acEff=%f\n", acEff);
//        printf("errAutoCorr=%.10f\n", autocorrelation1DLag1<double>(diff, num_elements, diff_sum / num_elements));
        free(diff);
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements) {
        double psnr, nrmse, max_diff;
        verify(ori_data, data, num_elements, psnr, nrmse, max_diff);
    }

    template<typename Type>
    void verify(Type *ori_data, Type *data, size_t num_elements, double &psnr, double &nrmse) {
        double max_diff;
        verify(ori_data, data, num_elements, psnr, nrmse, max_diff);
    }
};


#endif
