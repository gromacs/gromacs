#ifndef SZ_optimize_quant_intervals_hpp
#define SZ_optimize_quant_intervals_hpp

#include <vector>
#include "SZ3/predictor/MetaLorenzoPredictor.hpp"

namespace SZ3 {

#define QuantIntvMeanCapacity 8192
#define QuantIntvSampleDistance 100
#define QuantIntvSampleCapacity 32768
#define QuantIntvAccThreshold 0.999

// copied from conf.c
    unsigned int
    round_up_power_of_2(unsigned int base) {
        base -= 1;
        base = base | (base >> 1);
        base = base | (base >> 2);
        base = base | (base >> 4);
        base = base | (base >> 8);
        base = base | (base >> 16);
        return base + 1;
    }

    unsigned int
    static estimate_quantization_intervals(const std::vector<size_t> &intervals, size_t sample_count) {
        size_t target = sample_count * QuantIntvAccThreshold;
        size_t sum = 0;
        size_t i = 0;
        for (i = 0; i < intervals.size(); i++) {
            sum += intervals[i];
            if (sum > target)
                break;
        }
        if (i == intervals.size()) i = intervals.size() - 1;
        unsigned int accIntervals = 2 * (i + 1);
        unsigned int num_intervals = 2 * round_up_power_of_2(accIntervals);
        return (num_intervals > 32) ? num_intervals : 32;
    }

    template<typename T>
    float estimate_mean_freq_and_position(const std::vector<size_t> &freq_intervals, double precision, size_t sample_count,
                                    T &mean_guess) {
        size_t max_sum = 0;
        size_t max_index = 0;
        size_t tmp_sum = 0;
        for (size_t i = 1; i < freq_intervals.size() - 2; i++) {
            tmp_sum = freq_intervals[i] + freq_intervals[i + 1];
            if (tmp_sum > max_sum) {
                max_sum = tmp_sum;
                max_index = i;
            }
        }
        mean_guess += precision * (ptrdiff_t) (max_index + 1 - (freq_intervals.size() >> 1));
        return max_sum * 1.0 / sample_count;
    }

    template<typename T>
    float
    sample_rough_mean_3d(const T *data, size_t r1, size_t r2, size_t r3, size_t sample_distance) {
        double mean = 0;
        size_t len = r1 * r2 * r3;
        const T *data_pos = data;
        size_t offset_count = 0;
        size_t offset_count_2 = 0;
        size_t mean_count = 0;
        while (data_pos - data < len) {
            mean += *data_pos;
            mean_count++;
            data_pos += sample_distance;
            offset_count += sample_distance;
            offset_count_2 += sample_distance;
            if (offset_count >= r3) {
                offset_count = 0;
                data_pos -= 1;
            }
            if (offset_count_2 >= r2 * r3) {
                offset_count_2 = 0;
                data_pos -= 1;
            }
        }
        if (mean_count > 0) mean /= mean_count;
        return mean;
    }

    template<typename T>
    int optimize_quant_invl_3d(const T *data, size_t r1, size_t r2, size_t r3, double precision, float &pred_freq, float &mean_freq, T &mean_guess) {
        float mean_rough = sample_rough_mean_3d(data, r1, r2, r3, sqrt(r1 * r2 * r3));
        std::vector<size_t> intervals = std::vector<size_t>(QuantIntvSampleCapacity, 0);
        std::vector<size_t> freq_intervals = std::vector<size_t>(QuantIntvMeanCapacity, 0);
        size_t freq_count = 0;
        size_t sample_count = 0;
        size_t sample_distance = QuantIntvSampleDistance;
        size_t offset_count = sample_distance - 2; // count r3 offset
        size_t offset_count_2 = 0;
        size_t r23 = r2 * r3;
        size_t len = r1 * r23;
        const T *data_pos = data + r23 + r3 + offset_count;
        size_t n1_count = 1, n2_count = 1; // count i,j sum
        T pred_value = 0;
        double mean_diff = 0;
        ptrdiff_t freq_index = 0;
        size_t pred_index = 0;
        float pred_err = 0;
        int radius = (QuantIntvMeanCapacity >> 1);
        while (data_pos - data < len) {
            pred_value = SZMETA::lorenzo_predict_3d(data_pos, r23, r3);
            pred_err = fabs(pred_value - *data_pos);
            if (pred_err < precision) freq_count++;
            pred_index = (pred_err / precision + 1) / 2;
            if (pred_index >= intervals.size()) {
                pred_index = intervals.size() - 1;
            }
            intervals[pred_index]++;

            mean_diff = *data_pos - mean_rough;
            if (mean_diff > 0) freq_index = (ptrdiff_t) (mean_diff / precision) + radius;
            else freq_index = (ptrdiff_t) (mean_diff / precision) - 1 + radius;
            if (freq_index <= 0) {
                freq_intervals[0]++;
            } else if (freq_index >= freq_intervals.size()) {
                freq_intervals[freq_intervals.size() - 1]++;
            } else {
                freq_intervals[freq_index]++;
            }
            offset_count += sample_distance;
            if (offset_count >= r3) {
                n2_count++;
                if (n2_count == r2) {
                    n1_count++;
                    n2_count = 1;
                    data_pos += r3;
                }
                offset_count_2 = (n1_count + n2_count) % sample_distance;
                data_pos += (r3 + sample_distance - offset_count) + (sample_distance - offset_count_2);
                offset_count = (sample_distance - offset_count_2);
                if (offset_count == 0) offset_count++;
            } else data_pos += sample_distance;
            sample_count++;
        }
        pred_freq = freq_count * 1.0 / sample_count;
        mean_guess = mean_rough;
        mean_freq = estimate_mean_freq_and_position(freq_intervals, precision, sample_count, mean_guess);
        return estimate_quantization_intervals(intervals, sample_count);
    }
}

#endif