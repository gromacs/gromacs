#ifndef SZ3_SZALGO_LORENZO_REG_HPP
#define SZ3_SZALGO_LORENZO_REG_HPP

#include <cmath>
#include <memory>

#include "SZ3/compressor/SZGenericCompressor.hpp"
#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/QuantOptimization.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {
template <class T, uint N, class Quantizer, class Encoder, class Lossless>
std::shared_ptr<concepts::CompressorInterface<T>> make_compressor_lorenzo_regression(const Config &conf,
                                                                                     Quantizer quantizer,
                                                                                     Encoder encoder,
                                                                                     Lossless lossless) {
    std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors;

    int methodCnt = (conf.lorenzo + conf.lorenzo2 + conf.regression);
    int use_single_predictor = (methodCnt == 1);
    if (methodCnt == 0) {
        throw std::invalid_argument("All lorenzo and regression methods are disabled.");
    }
    if (conf.lorenzo) {
        if (use_single_predictor) {
            return make_compressor_sz_generic<T, N>(
                make_decomposition_blockwise<T, N>(conf, LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer),
                encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        }
    }
    if (conf.lorenzo2) {
        if (use_single_predictor) {
            return make_compressor_sz_generic<T, N>(
                make_decomposition_blockwise<T, N>(conf, LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer),
                encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return make_compressor_sz_generic<T, N>(
                make_decomposition_blockwise<T, N>(
                    conf, RegressionPredictor<T, N, Encoder>(conf.blockSize, conf.absErrorBound), quantizer),
                encoder, lossless);
        } else {
            predictors.push_back(
                std::make_shared<RegressionPredictor<T, N, Encoder>>(conf.blockSize, conf.absErrorBound));
        }
    }

    return make_compressor_sz_generic<T, N>(
        make_decomposition_blockwise<T, N>(conf, ComposedPredictor<T, N, Encoder>(predictors), quantizer), encoder,
        lossless);
}

template <class T, uint N>
size_t SZ_compress_LorenzoReg(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_LORENZO_REG);
    calAbsErrorBound(conf, data);

    auto quantizer = LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
    auto sz = make_compressor_lorenzo_regression<T, N>(conf, quantizer, HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N, class Encoder = HuffmanEncoder<int>>
void SZ_decompress_LorenzoReg(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_LORENZO_REG);
    auto cmpDataPos = cmpData;
    LinearQuantizer<T> quantizer;
    auto sz = make_compressor_lorenzo_regression<T, N>(conf, quantizer, Encoder(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}
}  // namespace SZ3
#endif
