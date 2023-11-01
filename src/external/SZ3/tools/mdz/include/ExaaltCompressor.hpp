#ifndef _SZ_EXAALT_COMPRESSSOR_HPP
#define _SZ_EXAALT_COMPRESSSOR_HPP

#include <SZ3/predictor/Predictor.hpp>
#include <SZ3/predictor/LorenzoPredictor.hpp>
#include <SZ3/quantizer/Quantizer.hpp>
#include <SZ3/encoder/Encoder.hpp>
#include <SZ3/lossless/Lossless.hpp>
#include <SZ3/utils/Iterator.hpp>
#include <SZ3/utils/MemoryUtil.hpp>
#include <SZ3/utils/Config.hpp>
#include <SZ3/def.hpp>

namespace SZ3 {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZ_Exaalt_Compressor {
    public:


        SZ_Exaalt_Compressor(const Config &conf,
                             Quantizer quantizer,
                             Encoder encoder,
                             Lossless lossless,
                             int timestep_op) :
                quantizer(quantizer), encoder(encoder), lossless(lossless),
                num_elements(conf.num), timestep_op(timestep_op) {
            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
            assert(!(timestep_op > 0 && conf.dims.size() != 2) &&
                   "timestep prediction requires 2d dataset");
            global_dimensions[0] = conf.dims[0];
            global_dimensions[1] = conf.dims[1];
        }

        void set_level(float level_start_, float level_offset_, int level_num_) {
            this->level_start = level_start_;
            this->level_offset = level_offset_;
            this->level_num = level_num_ + 200;
        }

        inline int quantize_to_level(T data) {
            return round((data - level_start) / level_offset);
        }

        inline T level(int l) {
            return level_start + l * level_offset;
        }

        // compress given the error bound
        uchar *compress(T *data, size_t &compressed_size) {

            std::vector<int> quant_inds(num_elements);
            std::vector<int> pred_inds(num_elements);
            quantizer.precompress_data();

//            Timer timer(true);

            auto l0 = quantize_to_level(data[0]);
            pred_inds[0] = l0 + level_num;
            quant_inds[0] = quantizer.quantize_and_overwrite(data[0], level(l0));

            if (timestep_op == 0) {
                for (size_t i = 1; i < num_elements; i++) {
                    auto l = quantize_to_level(data[i]);
                    pred_inds[i] = l - l0 + level_num;
                    quant_inds[i] = quantizer.quantize_and_overwrite(data[i], level(l));
                    l0 = l;
                }
            } else {

                std::vector<int> levels(global_dimensions[1]);
                levels[0] = l0;
                for (size_t i = 1; i < global_dimensions[1]; i++) {
                    levels[i] = quantize_to_level(data[i]);
                    pred_inds[i] = levels[i] - levels[i - 1] + level_num;
                    quant_inds[i] = quantizer.quantize_and_overwrite(data[i], level(levels[i]));
                }
                auto pred_idx = global_dimensions[1];
                if (timestep_op == 1) {
                    for (size_t i = 0; i < global_dimensions[1]; i++) {
                        for (size_t t = 1; t < global_dimensions[0]; t++) {
                            size_t idx = t * global_dimensions[1] + i;
                            quant_inds[pred_idx++] = quantizer.quantize_and_overwrite(
                                    data[idx], data[(t - 1) * global_dimensions[1] + i]);
                        }
                    }
                    pred_inds.resize(global_dimensions[1]);
                } else {
                    for (size_t i = 0; i < global_dimensions[1]; i++) {
                        l0 = levels[i];
                        for (size_t t = 1; t < global_dimensions[0]; t++) {
                            size_t idx = t * global_dimensions[1] + i;
                            auto l = quantize_to_level(data[idx]);
                            pred_inds[pred_idx] = l - l0 + level_num;
                            quant_inds[pred_idx++] = quantizer.quantize_and_overwrite(
                                    data[idx], level(l));
                            l0 = l;
                        }
                    }
                }
                assert(pred_idx == num_elements);
            }

//            timer.stop("Predition & Quantization");

            quantizer.postcompress_data();

            uchar *compressed_data;
            compressed_data = new uchar[4 * num_elements * sizeof(T)];
            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            quantizer.save(compressed_data_pos);
//            quantizer.print();

            encoder.preprocess_encode(quant_inds, 2 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

//            std::cout << *std::min_element(pred_inds.begin(), pred_inds.end()) << std::endl;
//            std::cout << *std::max_element(pred_inds.begin(), pred_inds.end()) << std::endl;

            encoder.preprocess_encode(pred_inds, level_num * 2 + 1);
            encoder.save(compressed_data_pos);
            encoder.encode(pred_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data, compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }

        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            size_t remaining_length = length;
            auto compressed_data = lossless.decompress(lossless_compressed_data, remaining_length);
            uchar const *compressed_data_pos = compressed_data;
            read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
//                std::cout << d << " ";
            }
//            std::cout << std::endl;
            quantizer.load(compressed_data_pos, remaining_length);
            encoder.load(compressed_data_pos, remaining_length);
            auto quant_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();

            encoder.load(compressed_data_pos, remaining_length);
            auto pred_inds_num = (timestep_op == 1) ? global_dimensions[1] : num_elements;
            auto pred_inds = encoder.decode(compressed_data_pos, pred_inds_num);
            encoder.postprocess_decode();

            lossless.postdecompress_data(compressed_data);

            auto dec_data = new T[num_elements];

            quantizer.predecompress_data();

            auto l = pred_inds[0] - level_num;
            dec_data[0] = quantizer.recover(level(l), quant_inds[0]);

            if (timestep_op == 0) {
                for (size_t i = 1; i < num_elements; i++) {
                    l += pred_inds[i] - level_num;
                    dec_data[i] = quantizer.recover(level(l), quant_inds[i]);
                }
            } else {
                std::vector<int> levels(global_dimensions[1]);
                levels[0] = l;
                for (size_t i = 1; i < global_dimensions[1]; i++) {
                    l += pred_inds[i] - level_num;
                    dec_data[i] = quantizer.recover(level(l), quant_inds[i]);
                    levels[i] = l;
                }
                auto pred_idx = global_dimensions[1];
                if (timestep_op == 1) {
                    for (size_t i = 0; i < global_dimensions[1]; i++) {
                        for (size_t t = 1; t < global_dimensions[0]; t++) {
                            dec_data[t * global_dimensions[1] + i] =
                                    quantizer.recover(
                                            dec_data[(t - 1) * global_dimensions[1] + i],
                                            quant_inds[pred_idx++]);
                        }
                    }
                } else {
                    for (size_t i = 0; i < global_dimensions[1]; i++) {
                        l = levels[i];
                        for (size_t t = 1; t < global_dimensions[0]; t++) {
                            size_t idx = t * global_dimensions[1] + i;
                            l += pred_inds[pred_idx] - level_num;
                            dec_data[idx] = quantizer.recover(level(l), quant_inds[pred_idx++]);
                        }
                    }
                }
                assert(pred_idx == num_elements);
            }


            quantizer.postdecompress_data();
            encoder.postprocess_decode();
            return dec_data;
        }


    private:
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        float level_start;
        float level_offset;
        int level_num;
        int timestep_op;
    };


}
#endif

