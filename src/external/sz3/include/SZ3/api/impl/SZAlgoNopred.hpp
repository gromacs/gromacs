#ifndef SZ3_SZALGO_NOPRED_HPP
#define SZ3_SZALGO_NOPRED_HPP

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/NoPredictionDecomposition.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_nopred(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_NOPRED);
    calAbsErrorBound(conf, data);

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_noprediction<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
    //        return cmpData;
}

template <class T, uint N>
void SZ_decompress_nopred(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_NOPRED);
    auto cmpDataPos = cmpData;
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_noprediction<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

}  // namespace SZ3
#endif
