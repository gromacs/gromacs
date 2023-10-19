#ifndef SZ3_IMPL_SZDISPATCHER_OMP_HPP
#define SZ3_IMPL_SZDISPATCHER_OMP_HPP

#include "SZ3/api/impl/SZDispatcher.hpp"
#include <cmath>
#include <memory>


#include "omp.h"

template<class T, SZ::uint N>
char *SZ_compress_OMP(SZ::Config &conf, const T *data, size_t &outSize) {
    unsigned char *buffer, *buffer_pos;

#ifdef _OPENMP

    assert(N == conf.N);

    std::vector<char *> compressed_t;
    std::vector<size_t> cmp_size_t, cmp_start_t;
    std::vector<T> min_t, max_t;
    std::vector<SZ::Config> conf_t;
//    SZ::Timer timer(true);
    int nThreads = 1;
    double eb;
#pragma omp parallel
    {
#pragma omp single
        {
            nThreads = omp_get_num_threads();
            if (conf.dims[0] < nThreads) {
                nThreads = conf.dims[0];
            }
            printf("OpenMP threads = %d\n", nThreads);
            compressed_t.resize(nThreads);
            cmp_size_t.resize(nThreads + 1);
            cmp_start_t.resize(nThreads + 1);
            conf_t.resize(nThreads);
            min_t.resize(nThreads);
            max_t.resize(nThreads);
        }


        int tid = omp_get_thread_num();

        auto dims_t = conf.dims;
        int lo = tid * conf.dims[0] / nThreads;
        int hi = (tid + 1) * conf.dims[0] / nThreads;
        dims_t[0] = hi - lo;
        auto it = dims_t.begin();
        size_t num_t_base = std::accumulate(++it, dims_t.end(), (size_t) 1, std::multiplies<size_t>());
        size_t num_t = dims_t[0] * num_t_base;

//        T *data_t = data + lo * num_t_base;
        std::vector<T> data_t(data + lo * num_t_base, data + lo * num_t_base + num_t);
        if (conf.errorBoundMode != SZ::EB_ABS) {
            auto minmax = std::minmax_element(data_t.begin(), data_t.end());
            min_t[tid] = *minmax.first;
            max_t[tid] = *minmax.second;
#pragma omp barrier
#pragma omp single
            {
                T range = *std::max_element(max_t.begin(), max_t.end()) - *std::min_element(min_t.begin(), min_t.end());
                SZ::calAbsErrorBound<T>(conf, data, range);
//                timer.stop("OMP init");
//                timer.start();
            }
        }

        conf_t[tid] = conf;
        conf_t[tid].setDims(dims_t.begin(), dims_t.end());
        compressed_t[tid] = SZ_compress_dispatcher<T, N>(conf_t[tid], data_t.data(), cmp_size_t[tid]);

#pragma omp barrier
#pragma omp single
        {
//            timer.stop("OMP compression");
//            timer.start();
            cmp_start_t[0] = 0;
            for (int i = 1; i <= nThreads; i++) {
                cmp_start_t[i] = cmp_start_t[i - 1] + cmp_size_t[i - 1];
            }
            size_t bufferSize = sizeof(int) + (nThreads + 1) * SZ::Config::size_est() + cmp_start_t[nThreads];
            buffer = new SZ::uchar[bufferSize];
            buffer_pos = buffer;
            SZ::write(nThreads, buffer_pos);
            for (int i = 0; i < nThreads; i++) {
                conf_t[i].save(buffer_pos);
            }
            SZ::write(cmp_size_t.data(), nThreads, buffer_pos);
        }

        memcpy(buffer_pos + cmp_start_t[tid], compressed_t[tid], cmp_size_t[tid]);
        delete[] compressed_t[tid];
    }

    outSize = buffer_pos - buffer + cmp_start_t[nThreads];
//    timer.stop("OMP memcpy");

#endif
    return (char *) buffer;
}


template<class T, SZ::uint N>
void SZ_decompress_OMP(const SZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
#ifdef _OPENMP

    const unsigned char *cmpr_data_pos = (unsigned char *) cmpData;
    int nThreads = 1;
    SZ::read(nThreads, cmpr_data_pos);
    omp_set_num_threads(nThreads);
    printf("OpenMP threads = %d\n", nThreads);

    std::vector<SZ::Config> conf_t(nThreads);
    for (int i = 0; i < nThreads; i++) {
        conf_t[i].load(cmpr_data_pos);
    }

    std::vector<size_t> cmp_start_t, cmp_size_t;
    cmp_size_t.resize(nThreads);
    SZ::read(cmp_size_t.data(), nThreads, cmpr_data_pos);
    char *cmpr_data_p = cmpData + (cmpr_data_pos - (unsigned char *) cmpData);

    cmp_start_t.resize(nThreads + 1);
    cmp_start_t[0] = 0;
    for (int i = 1; i <= nThreads; i++) {
        cmp_start_t[i] = cmp_start_t[i - 1] + cmp_size_t[i - 1];
    }

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto dims_t = conf.dims;
        int lo = tid * conf.dims[0] / nThreads;
        int hi = (tid + 1) * conf.dims[0] / nThreads;
        dims_t[0] = hi - lo;
        auto it = dims_t.begin();
        size_t num_t_base = std::accumulate(++it, dims_t.end(), (size_t) 1, std::multiplies<size_t>());

        SZ_decompress_dispatcher<T, N>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid], decData + lo * num_t_base);
    }
#endif
}


#endif
