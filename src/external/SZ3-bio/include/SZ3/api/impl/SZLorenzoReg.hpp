#ifndef SZ3_SZ_LORENZO_REG_HPP
#define SZ3_SZ_LORENZO_REG_HPP

#include "SZ3/compressor/SZGeneralCompressor.hpp"
#include "SZ3/frontend/SZFastFrontend.hpp"
#include "SZ3/frontend/SZGeneralFrontend.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/predictor/ComposedPredictor.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/RegressionPredictor.hpp"
#include "SZ3/predictor/PolyRegressionPredictor.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/def.hpp"
#include <cmath>
#include <memory>

namespace SZ3 {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    std::shared_ptr<concepts::CompressorInterface<T>>
    make_lorenzo_regression_compressor(const Config &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) {
        std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors;

        int methodCnt = (conf.lorenzo + conf.lorenzo2 + conf.regression + conf.regression2);
        int use_single_predictor = (methodCnt == 1);
        if (methodCnt == 0) {
            printf("All lorenzo and regression methods are disabled.\n");
            exit(0);
        }
        if (conf.lorenzo) {
            if (use_single_predictor) {
                return make_sz_general_compressor<T, N>(
                        make_sz_general_frontend<T, N>(conf, LorenzoPredictor<T, N, 1>(conf.absErrorBound), quantizer),
                        encoder, lossless);
            } else {
                predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
            }
        }
        if (conf.lorenzo2) {
            if (use_single_predictor) {
                return make_sz_general_compressor<T, N>(
                        make_sz_general_frontend<T, N>(conf, LorenzoPredictor<T, N, 2>(conf.absErrorBound), quantizer),
                        encoder, lossless);
            } else {
                predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
            }
        }
        if (conf.regression) {
            if (use_single_predictor) {
                return make_sz_general_compressor<T, N>(
                        make_sz_general_frontend<T, N>(conf, RegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
            } else {
                predictors.push_back(std::make_shared<RegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
            }
        }

        if (conf.regression2) {
            if (use_single_predictor) {
                return make_sz_general_compressor<T, N>(
                        make_sz_general_frontend<T, N>(conf, PolyRegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
            } else {
                predictors.push_back(std::make_shared<PolyRegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
            }
        }
        return make_sz_general_compressor<T, N>(
                make_sz_general_frontend<T, N>(conf, ComposedPredictor<T, N>(predictors),
                                               quantizer), encoder, lossless);
    }


    template<class T, uint N>
    char *SZ_compress_LorenzoReg(Config &conf, T *data, size_t &outSize) {

        assert(N == conf.N);
        assert(conf.cmprAlgo == ALGO_LORENZO_REG);
        calAbsErrorBound(conf, data);

        char *cmpData;
        auto quantizer = LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
        if (N == 3 && !conf.regression2 || (N == 1 && !conf.regression && !conf.regression2)) {
            // use fast version for 3D
            auto sz = make_sz_general_compressor<T, N>(make_sz_fast_frontend<T, N>(conf, quantizer), HuffmanEncoder<int>(),
                                                       Lossless_zstd());
            cmpData = (char *) sz->compress(conf, data, outSize);
        } else {
            auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, HuffmanEncoder<int>(), Lossless_zstd());
            cmpData = (char *) sz->compress(conf, data, outSize);
        }
        return cmpData;
    }


    template<class T, uint N>
    void SZ_decompress_LorenzoReg(const Config &conf, char *cmpData, size_t cmpSize, T *decData) {
        assert(conf.cmprAlgo == ALGO_LORENZO_REG);

        uchar const *cmpDataPos = (uchar *) cmpData;
        LinearQuantizer<T> quantizer;
        if (N == 3 && !conf.regression2 || (N == 1 && !conf.regression && !conf.regression2)) {
            // use fast version for 3D
            auto sz = make_sz_general_compressor<T, N>(make_sz_fast_frontend<T, N>(conf, quantizer),
                                                       HuffmanEncoder<int>(), Lossless_zstd());
            sz->decompress(cmpDataPos, cmpSize, decData);
            return;

        } else {
            auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, HuffmanEncoder<int>(), Lossless_zstd());
            sz->decompress(cmpDataPos, cmpSize, decData);
            return;
        }

    }
}
#endif