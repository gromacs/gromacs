#ifndef SZ3_IMPL_SZDISPATCHER_OMP_HPP
#define SZ3_IMPL_SZDISPATCHER_OMP_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <memory>
#include <new>
#include <stdexcept>

#include "SZ3/api/impl/SZDispatcher.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"

#ifdef _OPENMP

#include <omp.h>

#endif
namespace SZ3 {
// Without OpenMP the pragmas below drop out and the same code runs as one thread, so a build
// without OpenMP still reads this chunked layout. Keep the code outside the pragmas.
template <class T, uint N>
size_t SZ_compress_OMP(Config& conf, const T* data, uchar* cmpData, size_t cmpCap) {
    int nThreads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    {
        nThreads = omp_get_num_threads();
    }
#endif
    // ALGO_BIOMD predicts each frame from the previous one, so {frames, atoms, 3} is split along the atoms and every
    // chunk keeps all frames; other data along its first dimension. A chunk is outer runs of (hi - lo) * inner values.
    const size_t split = conf.cmprAlgo == ALGO_BIOMD && conf.N == 3;
    const size_t ext = conf.dims[split], outer = split ? conf.dims[0] : 1, inner = conf.num / (outer * ext);
    if (ext < static_cast<size_t>(nThreads)) {
        nThreads = static_cast<int>(ext);
    }
    std::vector<T> min_t(nThreads), max_t(nThreads);
    std::vector<Config> conf_t(nThreads);
    std::vector<std::unique_ptr<uchar[]>> compressed_t(nThreads);
    std::vector<size_t> cmp_size_t(nThreads), cmp_start_t(nThreads + 1);

    // Each loop below hands out chunks, not thread ids: a region may get fewer threads than asked
    // for (nested parallelism, OMP_THREAD_LIMIT, dynamic adjustment), and every chunk still runs.
    if (conf.errorBoundMode != EB_ABS) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif
        for (int i = 0; i < nThreads; i++) {
            size_t lo = static_cast<size_t>(i) * conf.num / nThreads;
            size_t hi = static_cast<size_t>(i + 1) * conf.num / nThreads;
            auto minmax = std::minmax_element(data + lo, data + hi);
            min_t[i] = *minmax.first;
            max_t[i] = *minmax.second;
        }
        T range = *std::max_element(max_t.begin(), max_t.end()) - *std::min_element(min_t.begin(), min_t.end());
        calAbsErrorBound<T>(conf, data, range);
    }

    // An exception that leaves an OpenMP region does not unwind to the caller -- the runtime calls
    // std::terminate -- so each chunk stores its own and the first one is rethrown below.
    std::exception_ptr failure;
    // num_threads(nThreads) applies to this region alone. omp_set_num_threads(nThreads) would write
    // the process-wide nthreads-var instead, and that outlives the call: an application running 32
    // threads that compresses one chunk here would go on running 8 afterwards. SZ3 is a library
    // inside someone else's program and has no business changing that.
#ifdef _OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif
    for (int i = 0; i < nThreads; i++) try {
        size_t lo = static_cast<size_t>(i) * ext / nThreads;
        size_t hi = static_cast<size_t>(i + 1) * ext / nThreads;
        auto dims_t = conf.dims;
        dims_t[split] = hi - lo;
        conf_t[i] = conf;
        conf_t[i].setDims(dims_t.begin(), dims_t.end());
        // Room for the size header Lossless_zstd::compress writes ahead of the zstd stream.
        size_t cmp_size_cap = sizeof(uint64_t) + Lossless_zstd::compress_bound(conf_t[i].num * sizeof(T));
        compressed_t[i].reset(new uchar[cmp_size_cap]);
        const T* data_t = data + lo * inner;
        std::vector<T> gathered;
        if (outer > 1) {
            gathered.resize(conf_t[i].num);
            const size_t run = (hi - lo) * inner;
            for (size_t o = 0; o < outer; o++) std::copy_n(data + (o * ext + lo) * inner, run, &gathered[o * run]);
            data_t = gathered.data();
        }
        // we have to use conf_t[i].N instead of N since each chunk may be a slice of the original data
        if (conf_t[i].N == 1) {
            cmp_size_t[i] = SZ_compress_dispatcher<T, 1>(conf_t[i], data_t, compressed_t[i].get(), cmp_size_cap);
        } else if (conf_t[i].N == 2) {
            cmp_size_t[i] = SZ_compress_dispatcher<T, 2>(conf_t[i], data_t, compressed_t[i].get(), cmp_size_cap);
        } else if (conf_t[i].N == 3) {
            cmp_size_t[i] = SZ_compress_dispatcher<T, 3>(conf_t[i], data_t, compressed_t[i].get(), cmp_size_cap);
        } else if (conf_t[i].N == 4) {
            cmp_size_t[i] = SZ_compress_dispatcher<T, 4>(conf_t[i], data_t, compressed_t[i].get(), cmp_size_cap);
        } else {
            throw std::invalid_argument("Unsupported N");
        }
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) {
        std::rethrow_exception(failure);
    }

    size_t header_size = sizeof(int) + nThreads * sizeof(uint64_t);
    cmp_start_t[0] = 0;
    for (int i = 0; i < nThreads; i++) {
        header_size += conf_t[i].size_est();
        cmp_start_t[i + 1] = cmp_start_t[i] + cmp_size_t[i];
    }
    // cmpCap may come from a bound for fewer chunks, computed by a copy of this header built without OpenMP.
    if (header_size + cmp_start_t[nThreads] > cmpCap) {
        throw std::invalid_argument(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
    }
    uchar* buffer_pos = cmpData;
    write(nThreads, buffer_pos);
    for (int i = 0; i < nThreads; i++) {
        conf_t[i].save(buffer_pos);
    }
    for (int i = 0; i < nThreads; i++) write<uint64_t>(cmp_size_t[i], buffer_pos);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nThreads)
#endif
    for (int i = 0; i < nThreads; i++) {
        memcpy(buffer_pos + cmp_start_t[i], compressed_t[i].get(), cmp_size_t[i]);
    }
    return buffer_pos - cmpData + cmp_start_t[nThreads];
}

template <class T, uint N>
void SZ_decompress_OMP(Config& conf, const uchar* cmpData, size_t cmpSize, T* decData) {
    auto cmpr_data_pos = cmpData;
    const uchar* const cmp_end = cmpData + cmpSize;
    int nThreads = 1;
    if (static_cast<size_t>(cmp_end - cmpr_data_pos) < sizeof(nThreads))
        throw std::out_of_range("SZ3 OMP: truncated thread count");
    read(nThreads, cmpr_data_pos);
    // Each thread contributes at least a config and a size, so the count cannot exceed the buffer size.
    if (nThreads <= 0 || static_cast<size_t>(nThreads) > cmpSize)
        throw std::out_of_range("SZ3 OMP: invalid thread count");

    // Split as SZ_compress_OMP split it, from the same Config.
    const size_t split = conf.cmprAlgo == ALGO_BIOMD && conf.N == 3;
    const size_t ext = conf.dims[split], outer = split ? conf.dims[0] : 1, inner = conf.num / (outer * ext);

    std::vector<Config> conf_t(nThreads);
    for (int i = 0; i < nThreads; i++) {
        size_t confRemaining = static_cast<size_t>(cmp_end - cmpr_data_pos);
        conf_t[i].load(cmpr_data_pos, confRemaining);
    }

    if (conf_t[0].sz3MagicNumber != SZ3_MAGIC_NUMBER) {
        throw std::invalid_argument("magic number mismatch, the input data is not compressed by SZ3");
    }
    {
        // backward compatibility for v3.3.2; Config::load() does not read the version the dispatcher checks.
        for (auto& c : conf_t) c.sz3DataVer = conf.sz3DataVer;
    }

    std::vector<size_t> cmp_start_t, cmp_size_t;
    cmp_size_t.resize(nThreads);
    if (static_cast<size_t>(cmp_end - cmpr_data_pos) < static_cast<size_t>(nThreads) * sizeof(uint64_t))
        throw std::out_of_range("SZ3 OMP: truncated per-thread sizes");
    for (auto& s : cmp_size_t) {
        uint64_t v = 0;
        read(v, cmpr_data_pos);
        s = v;
    }
    auto cmpr_data_p = cmpr_data_pos;

    cmp_start_t.resize(nThreads + 1);
    cmp_start_t[0] = 0;
    // The payloads follow back-to-back; build the offsets without overflowing so no slice points out.
    const size_t payload_avail = static_cast<size_t>(cmp_end - cmpr_data_p);
    for (int i = 1; i <= nThreads; i++) {
        if (cmp_size_t[i - 1] > payload_avail - cmp_start_t[i - 1])
            throw std::out_of_range("SZ3 OMP: per-thread compressed sizes exceed the buffer");
        cmp_start_t[i] = cmp_start_t[i - 1] + cmp_size_t[i - 1];
    }

    std::exception_ptr failure;
#ifdef _OPENMP
#pragma omp parallel num_threads(nThreads)
#endif
    try {
        // nThreads is how many chunks the writer split the data into -- it came out of the stream,
        // not from this machine. num_threads asks for that many threads but nothing guarantees
        // them: OMP_THREAD_LIMIT caps the whole program, dynamic adjustment may lower it, and a
        // region nested inside the caller's own gets exactly one thread unless nested parallelism
        // was turned on, which it is not by default -- the case a filter called from inside an
        // application's parallel region lands in. Indexing the chunks by thread id would then skip
        // every chunk whose id no thread has, leaving that span of the output untouched and
        // reporting success, so each thread walks the chunks in strides of however many arrived.
#ifdef _OPENMP
        const int actual = omp_get_num_threads();
        const int first = omp_get_thread_num();
#else
        const int actual = 1;
        const int first = 0;
#endif
        for (int tid = first; tid < nThreads; tid += actual) {
            size_t lo = static_cast<size_t>(tid) * ext / nThreads;
            size_t hi = static_cast<size_t>(tid + 1) * ext / nThreads;
            const size_t run = (hi - lo) * inner;
            // A chunk's own Config decides how much is written into its share of decData.
            if (conf_t[tid].num != outer * run) {
                throw std::invalid_argument("SZ3: a chunk of the OpenMP stream does not match its share of the data");
            }
            T* dec_t = decData + lo * inner;
            std::vector<T> gathered;
            if (outer > 1) {
                gathered.resize(conf_t[tid].num);
                dec_t = gathered.data();
            }
            if (conf_t[tid].N == 1) {
                SZ_decompress_dispatcher<T, 1>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid], dec_t);
            } else if (conf_t[tid].N == 2) {
                SZ_decompress_dispatcher<T, 2>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid], dec_t);
            } else if (conf_t[tid].N == 3) {
                SZ_decompress_dispatcher<T, 3>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid], dec_t);
            } else if (conf_t[tid].N == 4) {
                SZ_decompress_dispatcher<T, 4>(conf_t[tid], cmpr_data_p + cmp_start_t[tid], cmp_size_t[tid], dec_t);
            } else {
                throw std::invalid_argument("Unsupported N");
            }
            if (outer > 1) {
                for (size_t o = 0; o < outer; o++)
                    std::copy_n(&gathered[o * run], run, decData + (o * ext + lo) * inner);
            }
        }
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) {
        std::rethrow_exception(failure);
    }
}

template <class T>
size_t SZ_compress_size_bound_omp(const Config& conf) {
    // A region without num_threads has at most omp_get_max_threads() threads, so this bounds the
    // chunk count SZ_compress_OMP picks, however many threads its region then gets.
    int nThreads = 1;
#ifdef _OPENMP
    nThreads = omp_get_max_threads();
#endif
    const size_t ext = conf.dims[conf.cmprAlgo == ALGO_BIOMD && conf.N == 3];
    if (ext < static_cast<size_t>(nThreads)) {
        nThreads = static_cast<int>(ext);
    }
    // the largest chunk: ceil(ext / nThreads) indices of the split dimension
    size_t chunk_size = (ext + nThreads - 1) / nThreads * (conf.num / ext);
    // for each thread, we save conf, compressed size, and compressed data
    // the per-chunk compressed data may carry the size header written by Lossless_zstd::compress
    return sizeof(int) + nThreads * conf.size_est() + 2 * nThreads * sizeof(uint64_t) +
           nThreads * Lossless_zstd::compress_bound(chunk_size * sizeof(T));
}
} // namespace SZ3

#endif
