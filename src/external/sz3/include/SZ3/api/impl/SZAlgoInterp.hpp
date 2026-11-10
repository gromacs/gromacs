#ifndef SZ3_SZALGO_INTERP_HPP
#define SZ3_SZALGO_INTERP_HPP

#include <memory>

#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/decomposition/BlockwiseDecomposition.hpp"
#include "SZ3/decomposition/InterpolationDecomposition.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Extraction.hpp"
#include "SZ3/utils/QuantOptimization.hpp"
#include "SZ3/utils/Sample.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_Interp(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    assert(conf.cmprAlgo == ALGO_INTERP);
    calAbsErrorBound(conf, data);
    if (conf.interpAnchorStride < 0) {  // set default anchor stride
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interpAnchorStride = anchor_strides[N - 1];
    }

    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        HuffmanEncoder<int>(), Lossless_zstd());
    return sz->compress(conf, data, cmpData, cmpCap);
}

template <class T, uint N, class Encoder = HuffmanEncoder<int>>
void SZ_decompress_Interp(const Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == ALGO_INTERP);
    auto cmpDataPos = cmpData;
    auto sz = make_compressor_sz_generic<T, N>(
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2)),
        Encoder(), Lossless_zstd());
    sz->decompress(conf, cmpDataPos, cmpSize, decData);
}

template <class T, uint N>
double interp_compress_test(
    const std::vector<std::vector<T>> sampled_blocks, const Config conf, uchar *cmpData,
    size_t cmpCap) {  // test interp cmp on a set of sampled data blocks and return the compression ratio
    auto sz =
        make_decomposition_interpolation<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));

    std::vector<int> total_quant_bins;
    for (size_t k = 0; k < sampled_blocks.size(); k++) {
        auto cur_block = sampled_blocks[k];
        auto quant_bins = sz.compress(conf, cur_block.data());
        total_quant_bins.insert(total_quant_bins.end(), quant_bins.begin(),
                                quant_bins.end());  // merge the quant bins. Lossless them together
    }

    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();
    encoder.preprocess_encode(total_quant_bins, sz.get_out_range().second);
    size_t bufferSize =
        std::max<size_t>(1000, 1.2 * (sz.size_est() + encoder.size_est() + sizeof(T) * total_quant_bins.size()));

    // Owned, because the encoder and the lossless layer below can throw.
    std::unique_ptr<uchar[]> buffer_owner(new uchar[bufferSize]);
    uchar *const buffer = buffer_owner.get();
    uchar *buffer_pos = buffer;
    sz.save(buffer_pos);
    encoder.save(buffer_pos);

    // store the size of quant_inds is necessary as it is not always equal to conf.num
    write<uint64_t>(total_quant_bins.size(), buffer_pos);
    encoder.encode(total_quant_bins, buffer_pos);
    encoder.postprocess_encode();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    auto compression_ratio = conf.num * sampled_blocks.size() * sizeof(T) * 1.0 / cmpSize;
    return compression_ratio;
}

template <class T, uint N>
double lorenzo_compress_test(
    const std::vector<std::vector<T>> sampled_blocks, const Config &conf, uchar *cmpData,
    size_t cmpCap) {  // test lorenzo cmp on a set of sampled data blocks and return the compression ratio
    std::vector<int> total_quant_bins;
    // if ((N == 3 && !conf.regression2) || (N == 1 && !conf.regression && !conf.regression2)) {
    std::vector<std::shared_ptr<concepts::PredictorInterface<T, N>>> predictors;
    predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
    predictors.push_back(std::make_shared<LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
    auto sz = make_decomposition_blockwise<T, N>(conf, ComposedPredictor<T, N>(predictors),
                                                 LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2));
    // auto sz = make_decomposition_lorenzo_regression<T, N>(conf, LinearQuantizer<T>(conf.absErrorBound,
    // conf.quantbinCnt / 2));
    for (size_t k = 0; k < sampled_blocks.size(); k++) {
        auto cur_block = sampled_blocks[k];
        auto quant_bins = sz.compress(conf, cur_block.data());
        total_quant_bins.insert(total_quant_bins.end(), quant_bins.begin(),
                                quant_bins.end());  // merge the quant bins. Lossless them together
    }
    auto encoder = HuffmanEncoder<int>();
    auto lossless = Lossless_zstd();
    encoder.preprocess_encode(total_quant_bins, conf.quantbinCnt);
    size_t bufferSize = std::max<size_t>(1000, 1.2 * (encoder.size_est() + sizeof(T) * total_quant_bins.size()));

    // Owned, because the encoder and the lossless layer below can throw.
    std::unique_ptr<uchar[]> buffer_owner(new uchar[bufferSize]);
    uchar *const buffer = buffer_owner.get();
    uchar *buffer_pos = buffer;
    sz.save(buffer_pos);
    encoder.save(buffer_pos);

    // store the size of quant_inds is necessary as it is not always equal to conf.num
    write<uint64_t>(total_quant_bins.size(), buffer_pos);
    encoder.encode(total_quant_bins, buffer_pos);
    encoder.postprocess_encode();
    auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
    auto compression_ratio = conf.num * sampled_blocks.size() * sizeof(T) * 1.0 / cmpSize;
    return compression_ratio;
    // }
    // else{
    //     return 0.0;
    // }
}

template <class T, uint N>
size_t SZ_compress_Interp_lorenzo(Config &conf, T *data, uchar *cmpData, size_t cmpCap) {
    assert(conf.cmprAlgo == ALGO_INTERP_LORENZO);

    //        Timer timer(true);
    calAbsErrorBound(conf, data);

    if (conf.interpAnchorStride < 0) {  // set default anchor stride
        std::array<size_t, 4> anchor_strides = {4096, 128, 32, 16};
        conf.interpAnchorStride = anchor_strides[N - 1];
    }

    std::array<double, 4> sample_Rates = {0.005, 0.005, 0.005,
                                          0.005};  // default data sample rate. todo: add a config var to control
    auto sampleRate = sample_Rates[N - 1];
    std::array<size_t, 4> sampleBlock_Sizes = {4096, 128, 32,
                                               16};  // default sampled data block rate. Should better be no smaller
                                                     // than the anchor stride. todo: add a config var to control
    size_t sampleBlockSize = sampleBlock_Sizes[N - 1];
    size_t shortest_edge = conf.dims[0];
    for (size_t i = 0; i < N; i++) {
        shortest_edge = conf.dims[i] < shortest_edge ? conf.dims[i] : shortest_edge;
    }
    // Automatically adjust sampleblocksize.
    while (sampleBlockSize >= shortest_edge) sampleBlockSize /= 2;
    while (sampleBlockSize >= 16 && (pow(sampleBlockSize + 1, N) / conf.num) > 1.5 * sampleRate) sampleBlockSize /= 2;
    if (sampleBlockSize < 8) sampleBlockSize = 8;

    bool to_tune = pow(sampleBlockSize + 1, N) <= 0.05 * conf.num;  // to further revise
    for (auto &dim : conf.dims) {
        if (dim < sampleBlockSize) {
            to_tune = false;
            break;
        }
    }

    if (!to_tune) {  // if the sampled data would be too many (currently it is 5% of the input), skip the tuning
        conf.cmprAlgo = ALGO_INTERP;
        return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    }
    std::vector<std::vector<T>> sampled_blocks;
    size_t per_block_ele_num = pow(sampleBlockSize + 1, N);
    size_t sampling_num;
    std::vector<std::vector<size_t>> starts;
    auto profStride = sampleBlockSize / 4;  // larger is faster, smaller is better
    profiling_block<T, N>(data, conf.dims, starts, sampleBlockSize, conf.absErrorBound,
                          profStride);  // filter out the non-constant data blocks
    size_t num_filtered_blocks = starts.size();
    bool profiling = num_filtered_blocks * per_block_ele_num >= 0.5 * sampleRate * conf.num;  // temp. to refine
    // bool profiling = false;
    sampleBlocks<T, N>(data, conf.dims, sampleBlockSize, sampled_blocks, sampleRate, profiling,
                       starts);  // sample out same data blocks
    sampling_num = sampled_blocks.size() * per_block_ele_num;

    if (sampling_num == 0 || sampling_num >= conf.num * 0.2) {
        conf.cmprAlgo = ALGO_INTERP;
        return SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    }
    double best_lorenzo_ratio = 0, best_interp_ratio = 0, ratio;
    size_t bufferCap = conf.num * sizeof(T);
    // Owned, because the compressions below can throw.
    std::unique_ptr<uchar[]> buffer_owner(new uchar[bufferCap]);
    uchar *const buffer = buffer_owner.get();
    Config lorenzo_config = conf;

    {
        // tune interp
        conf.interpDirection = 0;
        conf.interpAlpha = 1.25;
        conf.interpBeta = 2.0;
        auto testConfig = conf;
        std::vector<size_t> dims(N, sampleBlockSize + 1);
        testConfig.setDims(dims.begin(), dims.end());
        for (auto &interp_op : {INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC}) {
            testConfig.interpAlgo = interp_op;
            ratio = interp_compress_test<T, N>(sampled_blocks, testConfig, buffer, bufferCap);
            if (ratio > best_interp_ratio) {
                best_interp_ratio = ratio;
                conf.interpAlgo = interp_op;
            }
        }

        testConfig.interpAlgo = conf.interpAlgo;
        testConfig.interpDirection = factorial(N) - 1;
        ratio = interp_compress_test<T, N>(sampled_blocks, testConfig, buffer, bufferCap);
        if (ratio > best_interp_ratio * 1.02) {
            best_interp_ratio = ratio;
            conf.interpDirection = testConfig.interpDirection;
        }
        testConfig.interpDirection = conf.interpDirection;
        // test more alpha-beta pairs for best compression ratio,
        auto alphalist = std::vector<double>{1.0, 1.5, 2.0};
        auto betalist = std::vector<double>{1.0, 2.5, 3.0};
        for (size_t i = 0; i < alphalist.size(); i++) {
            auto alpha = alphalist[i];
            auto beta = betalist[i];
            testConfig.interpAlpha = alpha;
            testConfig.interpBeta = beta;
            ratio = interp_compress_test<T, N>(sampled_blocks, testConfig, buffer, bufferCap);
            if (ratio > best_interp_ratio * 1.02) {
                best_interp_ratio = ratio;
                conf.interpAlpha = alpha;
                conf.interpBeta = beta;
            }
        }
    }
    {
        // only test lorenzo for 1D
        if (N == 1 && best_interp_ratio < 50) {
            std::vector<size_t> sample_dims(N, sampleBlockSize + 1);
            lorenzo_config.cmprAlgo = ALGO_LORENZO_REG;
            lorenzo_config.setDims(sample_dims.begin(), sample_dims.end());
            lorenzo_config.lorenzo = true;
            lorenzo_config.lorenzo2 = true;
            lorenzo_config.regression = false;
            lorenzo_config.regression2 = false;
            lorenzo_config.openmp = false;
            lorenzo_config.blockSize = 5;
            //        lorenzo_config.quantbinCnt = 65536 * 2;
            best_lorenzo_ratio = lorenzo_compress_test<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
            //            delete[]cmprData;
            //    printf("Lorenzo ratio = %.2f\n", ratio);
        }
    }

    bool useInterp = !(best_lorenzo_ratio >= best_interp_ratio * 1.1 && best_lorenzo_ratio < 50 &&
                       best_interp_ratio < 50);  // 1.1 is a fix coefficient. subject to revise
    size_t cmpSize = 0;
    if (useInterp) {
        conf.cmprAlgo = ALGO_INTERP;
        cmpSize = SZ_compress_Interp<T, N>(conf, data, cmpData, cmpCap);
    } else {
        // no need to tune lorenzo for 3D anymore
        // if (N == 3) {
        //     float pred_freq, mean_freq;
        //     T mean_guess;
        //     lorenzo_config.quantbinCnt = optimize_quant_invl_3d<T>(
        //         data, conf.dims[0], conf.dims[1], conf.dims[2], conf.absErrorBound, pred_freq, mean_freq,
        //         mean_guess);
        //     lorenzo_config.pred_dim = 2;
        //     ratio  =
        //         lorenzo_compress_test<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
        //     if (ratio > best_lorenzo_ratio * 1.02) {
        //         best_lorenzo_ratio = ratio;
        //     } else {
        //         lorenzo_config.pred_dim = 3;
        //     }
        // }

        if (conf.relErrorBound < 1.01e-6 && best_lorenzo_ratio > 5 && lorenzo_config.quantbinCnt != 16384) {
            auto quant_num = lorenzo_config.quantbinCnt;
            lorenzo_config.quantbinCnt = 16384;
            ratio = lorenzo_compress_test<T, N>(sampled_blocks, lorenzo_config, buffer, bufferCap);
            if (ratio > best_lorenzo_ratio * 1.02) {
                best_lorenzo_ratio = ratio;
            } else {
                lorenzo_config.quantbinCnt = quant_num;
            }
        }
        lorenzo_config.setDims(conf.dims.begin(), conf.dims.end());
        conf = lorenzo_config;
        //            double tuning_time = timer.stop();
        cmpSize = SZ_compress_LorenzoReg<T, N>(conf, data, cmpData, cmpCap);
    }

    return cmpSize;
}
}  // namespace SZ3
#endif
