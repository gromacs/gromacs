//
// Created by Kai Zhao on 12/9/19.
//

#ifndef _SZ_KMEANS_UTIL
#define _SZ_KMEANS_UTIL

#include <algorithm>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <random>
#include <vector>
#include <unordered_set>
#include <SZ3/utils/Timer.hpp>
#include <SZ3/def.hpp>


namespace SZ3 {

/*
 *  Internal implementation of the SMAWK algorithm.
 */
    template<typename T>
    void _smawk(
            const std::vector<size_t> &rows,
            const std::vector<size_t> &cols,
            const std::function<T(size_t, size_t)> &lookup,
            std::vector<size_t> *result) {
        // Recursion base case
        if (rows.size() == 0) return;

        // ********************************
        // * REDUCE
        // ********************************

        std::vector<size_t> _cols;  // Stack of surviving columns
        for (size_t col : cols) {
            while (true) {
                if (_cols.size() == 0) break;
                size_t row = rows[_cols.size() - 1];
                if (lookup(row, col) >= lookup(row, _cols.back()))
                    break;
                _cols.pop_back();
            }
            if (_cols.size() < rows.size())
                _cols.push_back(col);
        }

        // Call recursively on odd-indexed rows
        std::vector<size_t> odd_rows;
        for (size_t i = 1; i < rows.size(); i += 2) {
            odd_rows.push_back(rows[i]);
        }
        _smawk(odd_rows, _cols, lookup, result);

        std::unordered_map<size_t, size_t> col_idx_lookup;
        for (size_t idx = 0; idx < _cols.size(); ++idx) {
            col_idx_lookup[_cols[idx]] = idx;
        }

        // ********************************
        // * INTERPOLATE
        // ********************************

        // Fill-in even-indexed rows
        size_t start = 0;
        for (size_t r = 0; r < rows.size(); r += 2) {
            size_t row = rows[r];
            size_t stop = _cols.size() - 1;
            if (r < rows.size() - 1)
                stop = col_idx_lookup[(*result)[rows[r + 1]]];
            size_t argmin = _cols[start];
            T min = lookup(row, argmin);
            for (size_t c = start + 1; c <= stop; ++c) {
                T value = lookup(row, _cols[c]);
                if (c == start || value < min) {
                    argmin = _cols[c];
                    min = value;
                }
            }
            (*result)[row] = argmin;
            start = stop;
        }
    }

/*
 *  Interface for the SMAWK algorithm, for finding the minimum value in each row
 *  of an implicitly-defined totally monotone matrix.
 */
    template<typename T>
    std::vector<size_t> smawk(
            const size_t num_rows,
            const size_t num_cols,
            const std::function<T(size_t, size_t)> &lookup) {
        std::vector<size_t> result;
        result.resize(num_rows);
        std::vector<size_t> rows(num_rows);
        iota(begin(rows), end(rows), 0);
        std::vector<size_t> cols(num_cols);
        iota(begin(cols), end(cols), 0);
        _smawk<T>(rows, cols, lookup, &result);
        return result;
    }

/*
 *  Calculates cluster costs in O(1) using prefix sum arrays.
 */
    template<class DT>
    class CostCalculator {
        std::vector<double> cumsum;
        std::vector<double> cumsum2;

    public:
        CostCalculator(const std::vector<DT> &vec, size_t n) {
            cumsum.push_back(0.0);
            cumsum2.push_back(0.0);
            for (size_t i = 0; i < n; ++i) {
                double x = vec[i];
                cumsum.push_back(x + cumsum[i]);
                cumsum2.push_back(x * x + cumsum2[i]);
            }
        }

        double calc(size_t i, size_t j) {
            if (j < i) return 0.0;
            double mu = (cumsum[j + 1] - cumsum[i]) / (j - i + 1);
            double result = cumsum2[j + 1] - cumsum2[i];
            result += (j - i + 1) * (mu * mu);
            result -= (2 * mu) * (cumsum[j + 1] - cumsum[i]);
            return result;
        }
    };

    template<typename T>
    class Matrix {
        std::vector<T> data;
        size_t num_rows;
        size_t num_cols;

    public:
        Matrix(size_t num_rows, size_t num_cols) {
            this->num_rows = num_rows;
            this->num_cols = num_cols;
            data.resize(num_rows * num_cols);
        }

        inline T get(size_t i, size_t j) {
            return data[i * num_cols + j];
        }

        inline void set(size_t i, size_t j, T value) {
            data[i * num_cols + j] = value;
        }
    };

    template<class DT>
    void cluster(
            DT *array,
            size_t n,
            int &k,
            size_t *clusters,
            DT *centroids) {
        // ***************************************************
        // * Sort input array and save info for de-sorting
        // ***************************************************

        std::vector<size_t> sort_idxs(n);
        iota(sort_idxs.begin(), sort_idxs.end(), 0);
        sort(
                sort_idxs.begin(),
                sort_idxs.end(),
                [&array](size_t a, size_t b) { return array[a] < array[b]; });
//    vector<size_t> undo_sort_lookup(n);
        std::vector<DT> sorted_array(n);
        for (size_t i = 0; i < n; ++i) {
            sorted_array[i] = array[sort_idxs[i]];
//        undo_sort_lookup[sort_idxs[i]] = i;
        }

        // ***************************************************
        // * Set D and T using dynamic programming algorithm
        // ***************************************************

        // Algorithm as presented in section 2.2 of (Gronlund et al., 2017).

        CostCalculator<DT> cost_calculator(sorted_array, n);
        Matrix<DT> D(k, n);
        Matrix<size_t> T(k, n);

        for (size_t i = 0; i < n; ++i) {
            D.set(0, i, cost_calculator.calc(0, i));
            T.set(0, i, 0);
        }

        double ratio_avg = 0;
        bool findk = false;
        size_t bestk = 0;
        for (size_t k_ = 1; k_ < k; ++k_) {
            auto C = [&D, &k_, &cost_calculator](size_t i, size_t j) -> DT {
                size_t col = i < j - 1 ? i : j - 1;
                return D.get(k_ - 1, col) + cost_calculator.calc(j, i);
            };
            std::vector<size_t> row_argmins = smawk<DT>(n, n, C);
            for (size_t i = 0; i < row_argmins.size(); ++i) {
                size_t argmin = row_argmins[i];
                DT min = C(i, argmin);
                D.set(k_, i, min);
                T.set(k_, i, argmin);
            }
            float ratio = D.get(k_ - 1, n - 1) / D.get(k_, n - 1);
            ratio_avg = (ratio_avg * (k_ - 1) + ratio) / (k_);
//        std::cout << k_ + 1 << " , " << D.get(k_, n - 1) << " , " << ratio << " , " << ratio_avg << std::endl;
            if (ratio / ratio_avg > 1.5) {
                bestk = k_ + 1;
                findk = true;
            } else {
                if (findk) {
                    break;
                }
            }
        }

        if (!findk) {
            return;
        }
        k = bestk;
//        std::cout << "# groups = " << k << std::endl;

        // ***************************************************
        // * Extract cluster assignments by backtracking
        // ***************************************************

        // Note:  This step requires O(kn) memory usage due to saving the entire
        //       T matrix. However, it can be modified so that the memory usage is O(n).
        //       D and T would not need to be retained in full (D already doesn't need
        //       to be fully retained, although it currently is).
        //       Details are in section 3 of (Grønlund et al., 2017).

//    vector<DT> sorted_clusters(n);

        size_t t = n;
        size_t k_ = k - 1;
        size_t n_ = n - 1;
        // The do/while loop was used in place of:
        //   for (k_ = k - 1; k_ >= 0; --k_)
        // to avoid wraparound of an unsigned type.
        do {
            size_t t_ = t;
            t = T.get(k_, n_);
            DT centroid = 0.0;
            for (size_t i = t; i < t_; ++i) {
//            sorted_clusters[i] = k_;
                centroid += (sorted_array[i] - centroid) / (i - t + 1);
            }
            centroids[k_] = centroid;
            k_ -= 1;
            n_ = t - 1;
        } while (t > 0);

        // ***************************************************
        // * Order cluster assignments to match de-sorted
        // * ordering
        // ***************************************************

//    for (size_t i = 0; i < n; ++i) {
//        clusters[i] = sorted_clusters[undo_sort_lookup[i]];
//    }
    }

    template<class T>
    int f(T data, double start_position, double offset) {
        return round((data - start_position) / offset);
    }

    template<class T>
    int f1(T data, T *boundary, int n, double start_position, double offset) {
        return round((data - start_position) / offset);
    }

    template<class T>
    int f2(T data, T *boundary, int n, double start_position, double offset) {
        int low = 0, high = n; // numElems is the size of the array i.e arr.size()
        while (low != high) {
            int mid = (low + high) / 2; // Or a fancy way to avoid int overflow
            if (boundary[mid] <= data) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        return high;
    }

    template<class T>
    int f3(T data, T *boundary, int n, double start_position, double offset) {
        for (int i = 1; i < n; i++) {
            if (boundary[i] > data) {
                return i - 1;
            }
        }
        return n - 1;
    }


    template<class T>
    void get_cluster(T *data, size_t num, float &level_start, float &level_offset, int &level_num,
                     size_t sample_num) {
        T max = *std::max_element(data, data + num);
        T min = *std::min_element(data, data + num);
        SZ3::Timer timer;
        timer.start();
        std::vector<T> sample;
        if (num == sample_num) {
            sample = std::vector<T>(data, data + num);
        } else {
            sample.reserve(sample_num);
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//        std::uniform_int_distribution<> dis(0, 2 * sample_rate);
            size_t input_idx = 0;
//        for (int i = 0; i < num / sample_rate; i++) {
//            input_idx = (input_idx + dis(gen)) % num;
//            std::cout << input_idx << " ";
//            sample[i] = input[input_idx];
//        }
//        std::cout << std::endl;
            std::uniform_int_distribution<> dis2(0, num);
            std::unordered_set<size_t> sampledkeys;
//            printf("total_num=%lu, sample_num=%lu\n", num, sample_num);
            for (size_t i = 0; i < sample_num; i++) {
                do {
                    input_idx = dis2(gen);
                } while (sampledkeys.find(input_idx) != sampledkeys.end());
//            std::cout << input_idx << " ";
                sample[i] = data[input_idx];
            }
//        std::cout << std::endl;

//    std::sample(input.begin(), input.end(),
//                sample.begin(),
//                num / sample_rate,
//                std::mt19937{std::random_device{}()});
        }

        timer.stop("random sample");
//    sample = input;
        std::vector<size_t> idx(num);

        timer.start();
        int k = 150;
        std::vector<float> cents(k);
        cluster(sample.data(), sample_num, k, idx.data(), cents.data());
//    cluster(input.get(), num, 16, idx.data(), cents.data());
        timer.stop("kmeans1d");
        if (k == 150) {
//            std::cout << "No clusters are found." << std::endl;
            level_num = 0;
            return;
//            exit(0);
        }

//        std::cout << "centers : ";
//        for (size_t i = 0; i < k; i++) {
//            std::cout << cents[i] << " ";
//        }
//        std::cout << std::endl;

//        std::cout << "center diff : ";
//        for (size_t i = 1; i < k; i++) {
//            std::cout << cents[i] - cents[i - 1] << " ";
//        }
//        std::cout << std::endl;

//        std::vector<float> boundary(k);
//        boundary[0] = std::numeric_limits<float>::min();
//        for (size_t i = 1; i < k; i++) {
//            boundary[i] = (cents[i - 1] + cents[i]) / 2;
//        }

        level_offset = (cents[k - 1] - cents[0]) / (k - 1);
        level_start = cents[0];
        for (size_t i = 1; i < k; i++) {
            level_start += cents[i] - i * level_offset;
        }
        level_start /= k;
        level_num = f(max, level_start, level_offset) + 1;
    }
}

#endif