//
// Created by Kai Zhao on 4/20/20.
//

#ifndef SZ_EXTRACTION_HPP
#define SZ_EXTRACTION_HPP


namespace SZ3 {

    template<uint N>
    float cal_sampling_ratio(size_t block, size_t n, size_t dmin, std::vector<size_t> dims) {
        size_t sample_n = 1;
        for (auto dim: dims) {
            sample_n *= dim / dmin * 2 * block;
        }
        return sample_n * 1.0 / n;
    }


    template<class T, uint N>
    inline typename std::enable_if<N == 4, std::vector<T>>::type
    sampling(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, size_t &sampling_block) {
        assert(dims.size() == N);
        assert(sample_dims.size() == N);
        Timer timer(true);
        size_t num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());

        size_t dmin = *std::min_element(dims.begin(), dims.end());
        sampling_block = dmin;
        while (cal_sampling_ratio<N>(sampling_block, num, dmin, dims) > 0.035) {
            sampling_block--;
        }
        if (sampling_block * 2 > dmin) {
            sampling_block = dmin / 2;
        }
        if (sampling_block < 9) {
            sample_dims = dims;
            sample_num = num;
            return std::vector<T>();
        }
        size_t b0 = dims[0] / dmin;
        size_t b1 = dims[1] / dmin;
        size_t b2 = dims[2] / dmin;
        size_t b3 = dims[3] / dmin;
        sample_dims[0] = b0 * 2 * sampling_block;
        sample_dims[1] = b1 * 2 * sampling_block;
        sample_dims[2] = b2 * 2 * sampling_block;
        sample_dims[3] = b3 * 2 * sampling_block;
        size_t di, dj, dk, dt;
        sample_num = sample_dims[0] * sample_dims[1] * sample_dims[2] * sample_dims[3];
        std::vector<T> sampling_data(sample_num, 0);

        for (size_t bi = 0; bi < b0; bi++) {
            for (size_t bj = 0; bj < b1; bj++) {
                for (size_t bk = 0; bk < b2; bk++) {
                    for (size_t bt = 0; bt < b3; bt++) {
                        for (size_t i = 0; i < 2 * sampling_block; i++) {
                            for (size_t j = 0; j < 2 * sampling_block; j++) {
                                for (size_t k = 0; k < 2 * sampling_block; k++) {
                                    for (size_t t = 0; t < 2 * sampling_block; t++) {
                                        di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                                        dj = j < sampling_block ? j + sampling_block : dmin - 3 * sampling_block + j;
                                        dk = k < sampling_block ? k + sampling_block : dmin - 3 * sampling_block + k;
                                        dt = t < sampling_block ? t + sampling_block : dmin - 3 * sampling_block + t;
                                        size_t idx = (bi * dmin + di) * dims[1] * dims[2] * dims[3]
                                                     + (bj * dmin + dj) * dims[2] * dims[3]
                                                     + (bk * dmin + dk) * dims[3]
                                                     + bt * dmin + dt;
                                        sampling_data[(bi * 2 * sampling_block + i) * sample_dims[1] * sample_dims[2] * sample_dims[3]
                                                      + (bj * 2 * sampling_block + j) * sample_dims[1] * sample_dims[2]
                                                      + (bk * 2 * sampling_block + k) * sample_dims[2]
                                                      + bt * 2 * sampling_block + t] = data[idx];
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
        return sampling_data;
    }


    template<class T, uint N>
    inline typename std::enable_if<N == 3, std::vector<T>>::type
    sampling(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, size_t &sampling_block) {
        assert(dims.size() == N);
        assert(sample_dims.size() == N);
        size_t num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());

        size_t dmin = *std::min_element(dims.begin(), dims.end());
        sampling_block = dmin;
        while (cal_sampling_ratio<N>(sampling_block, num, dmin, dims) > 0.035) {
            sampling_block--;
        }
        if (sampling_block * 2 > dmin) {
            sampling_block = dmin / 2;
        }
        if (sampling_block < 9) {
            sample_dims = dims;
            sample_num = num;
            return std::vector<T>();
        }
        size_t b0 = dims[0] / dmin;
        size_t b1 = dims[1] / dmin;
        size_t b2 = dims[2] / dmin;
        sample_dims[0] = b0 * 2 * sampling_block;
        sample_dims[1] = b1 * 2 * sampling_block;
        sample_dims[2] = b2 * 2 * sampling_block;
        size_t di, dj, dk;
        sample_num = sample_dims[0] * sample_dims[1] * sample_dims[2];
        std::vector<T> sampling_data(sample_num, 0);

        for (size_t bi = 0; bi < b0; bi++) {
            for (size_t bj = 0; bj < b1; bj++) {
                for (size_t bk = 0; bk < b2; bk++) {
                    for (size_t i = 0; i < 2 * sampling_block; i++) {
                        for (size_t j = 0; j < 2 * sampling_block; j++) {
                            for (size_t k = 0; k < 2 * sampling_block; k++) {
                                di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                                dj = j < sampling_block ? j + sampling_block : dmin - 3 * sampling_block + j;
                                dk = k < sampling_block ? k + sampling_block : dmin - 3 * sampling_block + k;
                                auto d = data[(bi * dmin + di) * dims[1] * dims[2] + (bj * dmin + dj) * dims[2] +
                                              bk * dmin + dk];
                                sampling_data[(bi * 2 * sampling_block + i) * sample_dims[1] * sample_dims[2]
                                              + (bj * 2 * sampling_block + j) * sample_dims[2]
                                              + bk * 2 * sampling_block + k] = d;
                            }
                        }
                    }
                }
            }
        }
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
        return sampling_data;
    }

    template<class T, uint N>
    inline typename std::enable_if<N == 2, std::vector<T>>::type
    sampling(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, size_t &sampling_block) {
        assert(dims.size() == N);
        assert(sample_dims.size() == N);
        Timer timer(true);
        size_t num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());

        size_t dmin = *std::min_element(dims.begin(), dims.end());
        sampling_block = dmin;
        while (cal_sampling_ratio<N>(sampling_block, num, dmin, dims) > 0.035) {
            sampling_block--;
        }
        if (sampling_block * 2 > dmin) {
            sampling_block = dmin / 2;
        }
        if (sampling_block < 9) {
            sample_dims = dims;
            sample_num = num;
            return std::vector<T>();
        }
        size_t b0 = dims[0] / dmin;
        size_t b1 = dims[1] / dmin;
        sample_dims[0] = b0 * 2 * sampling_block;
        sample_dims[1] = b1 * 2 * sampling_block;
        size_t di, dj;
        sample_num = sample_dims[0] * sample_dims[1];
        std::vector<T> sampling_data(sample_num, 0);

        for (size_t bi = 0; bi < b0; bi++) {
            for (size_t bj = 0; bj < b1; bj++) {
                for (size_t i = 0; i < 2 * sampling_block; i++) {
                    for (size_t j = 0; j < 2 * sampling_block; j++) {
                        di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                        dj = j < sampling_block ? j + sampling_block : dmin - 3 * sampling_block + j;
                        auto d = data[(bi * dmin + di) * dims[1] + bj * dmin + dj];
                        sampling_data[(bi * 2 * sampling_block + i) * sample_dims[1]
                                      + bj * 2 * sampling_block + j] = d;
                    }
                }
            }
        }
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
        return sampling_data;
    }

    template<class T, uint N>
    inline typename std::enable_if<N == 1, std::vector<T>>::type
    sampling(T *data, std::vector<size_t> dims, size_t &sample_num, std::vector<size_t> &sample_dims, size_t &sampling_block) {
        assert(dims.size() == N);
        assert(sample_dims.size() == N);
        Timer timer(true);
        size_t num = std::accumulate(dims.begin(), dims.end(), (size_t) 1, std::multiplies<size_t>());

        size_t dmin = *std::min_element(dims.begin(), dims.end());
        sampling_block = dmin;
        while (cal_sampling_ratio<N>(sampling_block, num, dmin, dims) > 0.035) {
            sampling_block--;
        }
        if (sampling_block * 2 > dmin) {
            sampling_block = dmin / 2;
        }
        if (sampling_block < 9) {
            sample_dims = dims;
            sample_num = num;
            return std::vector<T>();
        }
        size_t b0 = dims[0] / dmin;
        sample_dims[0] = b0 * 2 * sampling_block;
        size_t di;
        sample_num = sample_dims[0];
        std::vector<T> sampling_data(sample_num, 0);

        for (size_t bi = 0; bi < b0; bi++) {
            for (size_t i = 0; i < 2 * sampling_block; i++) {
                di = i < sampling_block ? i + sampling_block : dmin - 3 * sampling_block + i;
                auto d = data[bi * dmin + di];
                sampling_data[bi * 2 * sampling_block + i] = d;
            }
        }
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
        return sampling_data;
    }
};


#endif
