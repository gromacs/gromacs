#ifndef SZ3_SZINTERP_HPP
#define SZ3_SZINTERP_HPP

#include "SZ3/compressor/SZInterpolationCompressor.hpp"
#include "SZ3/compressor/deprecated/SZBlockInterpolationCompressor.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimizatioin.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/impl/SZLorenzoReg.hpp"
#include <cmath>
#include <memory>

namespace SZ3 {
    template<class T, uint N>
    char *SZ_compress_Interp(Config &conf, T *data, size_t &outSize) {


        assert(N == conf.N);
        assert(conf.cmprAlgo == ALGO_INTERP);
        calAbsErrorBound(conf, data);

        auto sz = SZInterpolationCompressor<T, N, LinearQuantizer<T>, HuffmanEncoder<int>, Lossless_zstd>(
                LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2),
                HuffmanEncoder<int>(),
                Lossless_zstd());
        char *cmpData = (char *) sz.compress(conf, data, outSize);
        return cmpData;
    }


    template<class T, uint N>
    void SZ_decompress_Interp(const Config &conf, char *cmpData, size_t cmpSize, T *decData) {
        assert(conf.cmprAlgo == ALGO_INTERP);
        uchar const *cmpDataPos = (uchar *) cmpData;
        auto sz = SZInterpolationCompressor<T, N, LinearQuantizer<T>, HuffmanEncoder<int>, Lossless_zstd>(
                LinearQuantizer<T>(),
                HuffmanEncoder<int>(),
                Lossless_zstd());
        sz.decompress(cmpDataPos, cmpSize, decData);
    }


    template<class T, uint N>
    double do_not_use_this_interp_compress_block_test(T *data, std::vector<size_t> dims, size_t num,
                                                      double eb, int interp_op, int direction_op, int block_size) {

        std::vector<T> data1(data, data + num);
        size_t outSize = 0;

        Config conf;
        conf.absErrorBound = eb;
        conf.setDims(dims.begin(), dims.end());
        conf.blockSize = block_size;
        conf.interpAlgo = interp_op;
        conf.interpDirection = direction_op;
        auto sz = SZBlockInterpolationCompressor<T, N, LinearQuantizer<T>, HuffmanEncoder<int>, Lossless_zstd>(
                LinearQuantizer<T>(eb),
                HuffmanEncoder<int>(),
                Lossless_zstd());
        char *cmpData = (char *) sz.compress(conf, data1.data(), outSize);
        delete[]cmpData;
        auto compression_ratio = num * sizeof(T) * 1.0 / outSize;
        return compression_ratio;
    }

    template<class T, uint N>
    char *SZ_compress_Interp_lorenzo(Config &conf, T *data, size_t &outSize) {
        assert(conf.cmprAlgo == ALGO_INTERP_LORENZO);

        Timer timer(true);

        calAbsErrorBound(conf, data);

        size_t sampling_num, sampling_block;
        std::vector<size_t> sample_dims(N);
        std::vector<T> sampling_data = sampling<T, N>(data, conf.dims, sampling_num, sample_dims, sampling_block);
        if (sampling_num == conf.num) {
            conf.cmprAlgo = ALGO_INTERP;
            return SZ_compress_Interp<T, N>(conf, data, outSize);
        }

        double best_lorenzo_ratio = 0, best_interp_ratio = 0, ratio;
        size_t sampleOutSize;
        char *cmprData;
        Config lorenzo_config = conf;
        {
            //test lorenzo
            lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
            lorenzo_config.lorenzo = true;
            lorenzo_config.lorenzo2 = true;
            lorenzo_config.regression = false;
            lorenzo_config.regression2 = false;
            lorenzo_config.openmp = false;
            lorenzo_config.blockSize = 5;
//        lorenzo_config.quantbinCnt = 65536 * 2;
            std::vector<T> data1(sampling_data);
            cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, data1.data(), sampleOutSize);
            delete[]cmprData;
//    printf("Lorenzo ratio = %.2f\n", ratio);
            best_lorenzo_ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
        }

        {
            //tune interp
            for (auto &interp_op: {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC}) {
                ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                         interp_op, conf.interpDirection, sampling_block);
                if (ratio > best_interp_ratio) {
                    best_interp_ratio = ratio;
                    conf.interpAlgo = interp_op;
                }
            }

            int direction_op = factorial(N) - 1;
            ratio = do_not_use_this_interp_compress_block_test<T, N>(sampling_data.data(), sample_dims, sampling_num, conf.absErrorBound,
                                                                     conf.interpAlgo, direction_op, sampling_block);
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpDirection = direction_op;
            }
        }

        bool useInterp = !(best_lorenzo_ratio > best_interp_ratio && best_lorenzo_ratio < 80 && best_interp_ratio < 80);

        if (useInterp) {
            conf.cmprAlgo = ALGO_INTERP;
            double tuning_time = timer.stop();
            return SZ_compress_Interp<T, N>(conf, data, outSize);
        } else {
            //further tune lorenzo
            if (N == 3) {
                float pred_freq, mean_freq;
                T mean_guess;
                lorenzo_config.quantbinCnt = optimize_quant_invl_3d<T>(data, conf.dims[0], conf.dims[1], conf.dims[2],
                                                                            conf.absErrorBound, pred_freq, mean_freq, mean_guess);
                lorenzo_config.pred_dim = 2;
                cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.pred_dim = 3;
                }
            }

            if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
                auto quant_num = lorenzo_config.quantbinCnt;
                lorenzo_config.quantbinCnt = 16384;
                cmprData = SZ_compress_LorenzoReg<T, N>(lorenzo_config, sampling_data.data(), sampleOutSize);
                delete[]cmprData;
                ratio = sampling_num * 1.0 * sizeof(T) / sampleOutSize;
                if (ratio > best_lorenzo_ratio * 1.02) {
                    best_lorenzo_ratio = ratio;
                } else {
                    lorenzo_config.quantbinCnt = quant_num;
                }
            }
            lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
            conf = lorenzo_config;
            double tuning_time = timer.stop();
            return SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
        }


    }
}
#endif
