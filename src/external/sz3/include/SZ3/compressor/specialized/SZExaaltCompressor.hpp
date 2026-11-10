#ifndef SZ3_EXAALT_COMPRESSSOR_HPP
#define SZ3_EXAALT_COMPRESSSOR_HPP

#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

/**
 * specialized compressor for solid-state molecular dynamics simulation data
 * it utilizes the data character of levels to boost compression ratio
 * see ICDE'22 MDZ paper for details
 */
namespace SZ3 {
template <class T, uint N, class Quantizer, class Encoder, class Lossless>
class SZExaaltCompressor : public SZ3::concepts::CompressorInterface<T> {
   public:
    SZExaaltCompressor(Quantizer quantizer, Encoder encoder, Lossless lossless, int timestep_op)
        : quantizer(quantizer), encoder(encoder), lossless(lossless), timestep_op(timestep_op) {
        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quatizer interface");
        static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                      "must implement the encoder interface");
        static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                      "must implement the lossless interface");
    }

    // compress given the error bound
    //        uchar *compress(T *data, size_t &compressed_size) {
    size_t compress(const Config &conf, T *data, uchar *cmpData, size_t cmpCap) override {
        assert(!(timestep_op > 0 && conf.dims.size() != 2) && "timestep prediction requires 2d dataset");

        std::vector<int> quant_inds(conf.num);
        std::vector<int> pred_inds(conf.num);
        quantizer.precompress_data();

        //            Timer timer(true);

        auto l0 = quantize_to_level(data[0]);
        pred_inds[0] = l0 + level_num;
        quant_inds[0] = quantizer.quantize_and_overwrite(data[0], level(l0));

        if (timestep_op == 0) {
            for (size_t i = 1; i < conf.num; i++) {
                auto l = quantize_to_level(data[i]);
                pred_inds[i] = l - l0 + level_num;
                quant_inds[i] = quantizer.quantize_and_overwrite(data[i], level(l));
                l0 = l;
            }
        } else {
            std::vector<int> levels(conf.dims[1]);
            levels[0] = l0;
            for (size_t i = 1; i < conf.dims[1]; i++) {
                levels[i] = quantize_to_level(data[i]);
                pred_inds[i] = levels[i] - levels[i - 1] + level_num;
                quant_inds[i] = quantizer.quantize_and_overwrite(data[i], level(levels[i]));
            }
            auto pred_idx = conf.dims[1];
            if (timestep_op == 1) {
                for (size_t i = 0; i < conf.dims[1]; i++) {
                    for (size_t t = 1; t < conf.dims[0]; t++) {
                        size_t idx = t * conf.dims[1] + i;
                        quant_inds[pred_idx++] =
                            quantizer.quantize_and_overwrite(data[idx], data[(t - 1) * conf.dims[1] + i]);
                    }
                }
                pred_inds.resize(conf.dims[1]);
            } else {
                for (size_t i = 0; i < conf.dims[1]; i++) {
                    l0 = levels[i];
                    for (size_t t = 1; t < conf.dims[0]; t++) {
                        size_t idx = t * conf.dims[1] + i;
                        auto l = quantize_to_level(data[idx]);
                        pred_inds[pred_idx] = l - l0 + level_num;
                        quant_inds[pred_idx++] = quantizer.quantize_and_overwrite(data[idx], level(l));
                        l0 = l;
                    }
                }
            }
            assert(pred_idx == conf.num);
        }

        //            timer.stop("Predition & Quantization");

        quantizer.postcompress_data();

        auto buffer = static_cast<uchar *>(malloc(4 * conf.num * sizeof(T)));
        uchar *buffer_pos = buffer;
        quantizer.save(buffer_pos);
        //            quantizer.print();

        
        if (quantizer.get_out_range().first != 0) {
            throw std::runtime_error("The output range of the quantizer must start from 0 for this compressor");
        }
        encoder.preprocess_encode(quant_inds, quantizer.get_out_range().second);
        encoder.save(buffer_pos);
        encoder.encode(quant_inds, buffer_pos);
        encoder.postprocess_encode();

        //            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << std::endl;
        //            std::cout << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;

        encoder.preprocess_encode(pred_inds, level_num * 2 + 1);
        encoder.save(buffer_pos);
        encoder.encode(pred_inds, buffer_pos);
        encoder.postprocess_encode();

        auto cmpSize = lossless.compress(buffer, buffer_pos - buffer, cmpData, cmpCap);
        free(buffer);
        return cmpSize;
    }

    //        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
    T *decompress(const Config &conf, uchar const *cmpData, size_t cmpSize, T *dec_data) override {
        uchar *buffer = nullptr;
        size_t bufferSize = lossless.decompress(cmpData, cmpSize, buffer, 0);
        // The parsing below walks the decompressed buffer, so bufferSize is its bound, not cmpSize.
        size_t remaining_length = bufferSize;
        uchar const *buffer_pos = buffer;

        quantizer.load(buffer_pos, remaining_length);
        encoder.load(buffer_pos, remaining_length);
        auto quant_inds = encoder.decode(buffer_pos, conf.num, remaining_length);
        encoder.postprocess_decode();

        encoder.load(buffer_pos, remaining_length);
        auto pred_inds_num = (timestep_op == 1) ? conf.dims[1] : conf.num;
        auto pred_inds = encoder.decode(buffer_pos, pred_inds_num, remaining_length);
        encoder.postprocess_decode();

        free(buffer);

        quantizer.predecompress_data();

        auto l = pred_inds[0] - level_num;
        dec_data[0] = quantizer.recover(level(l), quant_inds[0]);

        if (timestep_op == 0) {
            for (size_t i = 1; i < conf.num; i++) {
                l += pred_inds[i] - level_num;
                dec_data[i] = quantizer.recover(level(l), quant_inds[i]);
            }
        } else {
            std::vector<int> levels(conf.dims[1]);
            levels[0] = l;
            for (size_t i = 1; i < conf.dims[1]; i++) {
                l += pred_inds[i] - level_num;
                dec_data[i] = quantizer.recover(level(l), quant_inds[i]);
                levels[i] = l;
            }
            auto pred_idx = conf.dims[1];
            if (timestep_op == 1) {
                for (size_t i = 0; i < conf.dims[1]; i++) {
                    for (size_t t = 1; t < conf.dims[0]; t++) {
                        dec_data[t * conf.dims[1] + i] =
                            quantizer.recover(dec_data[(t - 1) * conf.dims[1] + i], quant_inds[pred_idx++]);
                    }
                }
            } else {
                for (size_t i = 0; i < conf.dims[1]; i++) {
                    l = levels[i];
                    for (size_t t = 1; t < conf.dims[0]; t++) {
                        size_t idx = t * conf.dims[1] + i;
                        l += pred_inds[pred_idx] - level_num;
                        dec_data[idx] = quantizer.recover(level(l), quant_inds[pred_idx++]);
                    }
                }
            }
            assert(pred_idx == conf.num);
        }

        quantizer.postdecompress_data();
        encoder.postprocess_decode();
        return dec_data;
    }

    void set_level(float level_start_, float level_offset_, int level_num_) {
        this->level_start = level_start_;
        this->level_offset = level_offset_;
        this->level_num = level_num_ + 200;
    }

    inline int quantize_to_level(T data) { return round((data - level_start) / level_offset); }

    inline T level(int l) { return level_start + nofma(l * level_offset); }

   private:
    Quantizer quantizer;
    Encoder encoder;
    Lossless lossless;
    float level_start{};
    float level_offset{};
    int level_num{};
    int timestep_op;
};

template <class T, uint N, class Quantizer, class Encoder, class Lossless>
std::shared_ptr<SZExaaltCompressor<T, N, Quantizer, Encoder, Lossless>> make_compressor_exaalt(Quantizer quantizer,
                                                                                               Encoder encoder,
                                                                                               Lossless lossless,
                                                                                               int timestep_op) {
    return std::make_shared<SZExaaltCompressor<T, N, Quantizer, Encoder, Lossless>>(quantizer, encoder, lossless,
                                                                                    timestep_op);
}
}  // namespace SZ3
#endif
