#ifndef _SZ_BLOCK_INTERPOLATION_COMPRESSOR_HPP
#define _SZ_BLOCK_INTERPOLATION_COMPRESSOR_HPP

#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/Iterator.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Interpolators.hpp"
#include "SZ3/def.hpp"
#include <cstring>
#include <cmath>

namespace SZ3 {
    template<class T, uint N, class Quantizer, class Encoder, class Lossless>
    class SZBlockInterpolationCompressor {
    public:


        SZBlockInterpolationCompressor(Quantizer quantizer, Encoder encoder, Lossless lossless) :
                quantizer(quantizer), encoder(encoder), lossless(lossless) {

            static_assert(std::is_base_of<concepts::QuantizerInterface<T>, Quantizer>::value,
                          "must implement the quatizer interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");

        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {

            size_t remaining_length = cmpSize;

            uchar *buffer = lossless.decompress(cmpData, remaining_length);
            uchar const *buffer_pos = buffer;

            read(global_dimensions.data(), N, buffer_pos, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
            }
            read(block_size, buffer_pos, remaining_length);
            read(interpolator_id, buffer_pos, remaining_length);
            read(direction_sequence_id, buffer_pos, remaining_length);

            quantizer.load(buffer_pos, remaining_length);
            encoder.load(buffer_pos, remaining_length);
            quant_inds = encoder.decode(buffer_pos, num_elements);

            encoder.postprocess_decode();

            lossless.postdecompress_data(buffer);

            auto range = std::make_shared<multi_dimensional_range<T, N>>(decData,
                                                                             std::begin(global_dimensions),
                                                                             std::end(global_dimensions),
                                                                             block_size,
                                                                             0);

            quantizer.predecompress_data();

//            debug.resize(num_elements, 0);

            auto inter_begin = range->begin();
            auto inter_end = range->end();
            size_t block_idx = 0;

            for (auto block = inter_begin; block != inter_end; block++) {
                auto interp_end_idx = block.get_global_index();
                uint max_interp_level = 1;
                auto block_global_idx = block.get_global_index();
                for (int i = 0; i < N; i++) {
                    size_t block_dim = (block_global_idx[i] + block_size > global_dimensions[i]) ?
                                       global_dimensions[i] - block_global_idx[i] : block_size;
                    interp_end_idx[i] += block_dim - 1;
                    if (max_interp_level < ceil(log2(block_dim))) {
                        max_interp_level = (uint) ceil(log2(block_dim));
                    }
                }

                *block = quantizer.recover(0, quant_inds[quant_index++]);

                for (uint level = max_interp_level; level > 0 && level <= max_interp_level; level--) {
                    size_t stride_ip = 1U << (level - 1);
                    block_interpolation(decData, block.get_global_index(), interp_end_idx, PB_recover,
                                        interpolators[interpolator_id], direction_sequence_id, stride_ip);
                }
            }


//            assert(quant_index == num_elements);
            quantizer.postdecompress_data();

            return decData;
        }


        // compress given the error bound
        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {

            block_size = conf.blockSize;
            num_elements = conf.num;
            interpolator_id = conf.interpAlgo;
            direction_sequence_id = conf.interpDirection;


            std::copy_n(conf.dims.begin(), N, global_dimensions.begin());

            quant_inds.clear();
            auto range = std::make_shared<multi_dimensional_range<T, N>>(data,
                                                                             std::begin(global_dimensions),
                                                                             std::end(global_dimensions),
                                                                             block_size, 0);
            quantizer.precompress_data();
            for (auto block = range->begin(); block != range->end(); ++block) {

                auto block_global_idx = block.get_global_index();
                auto interp_end_idx = block.get_global_index();
                uint max_interp_level = 1;
                for (int i = 0; i < N; i++) {
                    size_t block_dim = (block_global_idx[i] + block_size > global_dimensions[i]) ?
                                       global_dimensions[i] - block_global_idx[i] : block_size;
                    interp_end_idx[i] += block_dim - 1;
                    if (max_interp_level < ceil(log2(block_dim))) {
                        max_interp_level = (uint) ceil(log2(block_dim));
                    }
                }
                quant_inds.push_back(quantizer.quantize_and_overwrite(*block, 0));

                for (uint level = max_interp_level; level > 0 && level <= max_interp_level; level--) {
                    uint stride_ip = 1U << (level - 1);
                    block_interpolation(data, block.get_global_index(), interp_end_idx, PB_predict_overwrite,
                                        interpolators[interpolator_id], direction_sequence_id, stride_ip);
                }
            }
            quantizer.postcompress_data();
//            predictor.print();

            encoder.preprocess_encode(quant_inds, quantizer.get_radius() * 2);
            size_t bufferSize = 1.2 * (quantizer.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());

            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;

            write(global_dimensions.data(), N, buffer_pos);
            write(block_size, buffer_pos);
            write(interpolator_id, buffer_pos);
            write(direction_sequence_id, buffer_pos);

            quantizer.save(buffer_pos);

            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();

            assert(buffer_pos - buffer < bufferSize);

            uchar *lossless_data = lossless.compress(buffer,
                                                     buffer_pos - buffer,
                                                     compressed_size);
            lossless.postcompress_data(buffer);

            return lossless_data;
        }

    private:

        enum PredictorBehavior {
            PB_predict_overwrite, PB_predict, PB_recover
        };

        size_t offset2(std::array<size_t, N> idx) {
            size_t offset = idx[0];
            for (int i = 1; i < N; i++) {
                offset = offset * global_dimensions[i] + idx[i];
            }
            return offset;
        }

        template<class... Args>
        size_t offset(Args &&... args) {
            std::array<size_t, N> idx{static_cast<size_t>(std::forward<Args>(args))...};
            return offset2(idx);
        }

        inline void quantize(T &d, T pred) {
            quant_inds.push_back(quantizer.quantize_and_overwrite(d, pred));
        }

        inline void recover(T &d, T pred) {
            d = quantizer.recover(pred, quant_inds[quant_index++]);
        };


        double block_interpolation_1d(T *data, size_t begin, size_t end, size_t stride,
                                      const std::string &interp_func,
                                      const PredictorBehavior pb) {
            size_t n = (end - begin) / stride + 1;
            if (n <= 1) {
                return 0;
            }
            double predict_error = 0;

            size_t stride3x = 3 * stride;
            size_t stride5x = 5 * stride;
            if (interp_func == "linear" || n < 5) {
                if (pb == PB_predict_overwrite) {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        quantize(*d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            quantize(*d, *(d - stride));
                        } else {
                            quantize(*d, interp_linear1(*(d - stride3x), *(d - stride)));
                        }
                    }
                } else {
                    for (size_t i = 1; i + 1 < n; i += 2) {
                        T *d = data + begin + i * stride;
                        recover(*d, interp_linear(*(d - stride), *(d + stride)));
                    }
                    if (n % 2 == 0) {
                        T *d = data + begin + (n - 1) * stride;
                        if (n < 4) {
                            recover(*d, *(d - stride));
                        } else {
                            recover(*d, interp_linear1(*(d - stride3x), *(d - stride)));
                        }
                    }
                }
            } else {
                if (pb == PB_predict_overwrite) {

                    T *d = data + begin + stride;
                    quantize(*d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (size_t i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        quantize(*d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    if (n % 2 == 0) {
                        d = data + begin + (n - 3) * stride;
                        quantize(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                        d += 2 * stride;
                        quantize(*d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    } else {
                        d = data + begin + (n - 2) * stride;
                        quantize(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    }
                } else {
                    T *d = data + begin + stride;
                    recover(*d, interp_quad_1(*(d - stride), *(d + stride), *(d + stride3x)));

                    for (size_t i = 3; i + 3 < n; i += 2) {
                        d = data + begin + i * stride;
                        recover(*d, interp_cubic(*(d - stride3x), *(d - stride), *(d + stride), *(d + stride3x)));
                    }
                    if (n % 2 == 0) {
                        d = data + begin + (n - 3) * stride;
                        recover(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                        d += 2 * stride;
                        recover(*d, interp_quad_3(*(d - stride5x), *(d - stride3x), *(d - stride)));
                    } else {
                        d = data + begin + (n - 2) * stride;
                        recover(*d, interp_quad_2(*(d - stride3x), *(d - stride), *(d + stride)));
                    }

                }
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 1, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            return block_interpolation_1d(data, offset2(begin), offset2(end), stride_ip, interp_func, pb);
        }

        template<uint NN = N>
        typename std::enable_if<NN == 2, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            double predict_error = 0;
            size_t stride_ip2 = stride_ip * 2;
            if (direction == 0) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                    predict_error += block_interpolation_1d(data, offset(begin[0], j), offset(end[0], j),
                                                            stride_ip * global_dimensions[1],
                                                            interp_func, pb);
                }
                for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                    predict_error += block_interpolation_1d(data, offset(i, begin[1]), offset(i, end[1]), stride_ip, interp_func,
                                                            pb);
                }
            } else {
                for (size_t i = begin[0]; i <= end[0]; i += stride_ip2) {
                    predict_error += block_interpolation_1d(data, offset(i, begin[1]), offset(i, end[1]), stride_ip, interp_func,
                                                            pb);
                }
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                    predict_error += block_interpolation_1d(data, offset(begin[0], j), offset(end[0], j),
                                                            stride_ip * global_dimensions[1],
                                                            interp_func, pb);
                }
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 4, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            double predict_error = 0;
            size_t stride_ip2 = stride_ip * 2;
            for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                    for (size_t t = (begin[3] ? begin[3] + stride_ip2 : 0); t <= end[3]; t += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(begin[0], j, k, t), offset(end[0], j, k, t),
                                                                stride_ip * global_dimensions[1] * global_dimensions[2] *
                                                                global_dimensions[3],
                                                                interp_func, pb);
                    }
                }
            }
            for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                    for (size_t t = (begin[3] ? begin[3] + stride_ip2 : 0); t <= end[3]; t += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(i, begin[1], k, t), offset(i, end[1], k, t),
                                                                stride_ip * global_dimensions[2] * global_dimensions[3],
                                                                interp_func, pb);
                    }
                }
            }
            for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                    for (size_t t = (begin[3] ? begin[3] + stride_ip2 : 0); t <= end[3]; t += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(i, j, begin[2], t), offset(i, j, end[2], t),
                                                                stride_ip * global_dimensions[3],
                                                                interp_func, pb);
                    }
                }
            }
            for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                    for (size_t k = begin[2]; k <= end[2]; k += stride_ip) {
                        predict_error += block_interpolation_1d(data, offset(i, j, k, begin[3]), offset(i, j, k, end[3]),
                                                                stride_ip, interp_func, pb);
                    }
                }
            }
            return predict_error;
        }

        template<uint NN = N>
        typename std::enable_if<NN == 3, double>::type
        block_interpolation(T *data, std::array<size_t, N> begin, std::array<size_t, N> end, const PredictorBehavior pb,
                            const std::string &interp_func, const int direction, uint stride_ip = 1) {
            double predict_error = 0;
            size_t stride_ip2 = stride_ip * 2;

            if (direction == 0 || direction == 1) {
                for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                    for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(begin[0], j, k), offset(end[0], j, k),
                                                                stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                interp_func,
                                                                pb);
                    }
                }
                if (direction == 0) {
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t k = begin[2]; k <= end[2]; k += stride_ip2) {
                            predict_error += block_interpolation_1d(data, offset(i, begin[1], k), offset(i, end[1], k),
                                                                    stride_ip * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t j = begin[1]; j <= end[1]; j += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(i, j, begin[2]), offset(i, j, end[2]), stride_ip,
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                } else {
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t j = begin[1]; j <= end[1]; j += stride_ip2) {
                            predict_error += block_interpolation_1d(data, offset(i, j, begin[2]), offset(i, j, end[2]), stride_ip,
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                    for (size_t i = begin[0]; i <= end[0]; i += stride_ip) {
                        for (size_t k = begin[2]; k <= end[2]; k += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(i, begin[1], k), offset(i, end[1], k),
                                                                    stride_ip * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                }

            } else if (direction == 2 || direction == 3) {
                for (size_t k = begin[0]; k <= end[0]; k += stride_ip2) {
                    for (size_t j = begin[2]; j <= end[2]; j += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(k, begin[1], j), offset(k, end[1], j),
                                                                stride_ip * global_dimensions[2],
                                                                interp_func,
                                                                pb);
                    }
                }
                if (direction == 2) {
                    for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                        for (size_t j = begin[2]; j <= end[2]; j += stride_ip2) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], i, j), offset(end[0], i, j),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                    for (size_t k = begin[0]; k <= end[0]; k += stride_ip) {
                        for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(k, i, begin[2]), offset(k, i, end[2]),
                                                                    stride_ip,
                                                                    interp_func, pb);
                        }
                    }
                } else {
                    for (size_t k = begin[0]; k <= end[0]; k += stride_ip2) {
                        for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(k, i, begin[2]), offset(k, i, end[2]),
                                                                    stride_ip,
                                                                    interp_func, pb);
                        }
                    }
                    for (size_t i = begin[1]; i <= end[1]; i += stride_ip) {
                        for (size_t j = begin[2]; j <= end[2]; j += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], i, j), offset(end[0], i, j),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func,
                                                                    pb);
                        }
                    }
                }

            } else if (direction == 4 || direction == 5) {

                for (size_t j = begin[0]; j <= end[0]; j += stride_ip2) {
                    for (size_t k = begin[1]; k <= end[1]; k += stride_ip2) {
                        predict_error += block_interpolation_1d(data, offset(j, k, begin[2]), offset(j, k, end[2]),
                                                                stride_ip, interp_func, pb);
                    }
                }
                if (direction == 4) {
                    for (size_t k = begin[1]; k <= end[1]; k += stride_ip2) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], k, i), offset(end[0], k, i),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                    for (size_t j = begin[0]; j <= end[0]; j += stride_ip) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(j, begin[1], i), offset(j, end[1], i),
                                                                    stride_ip * global_dimensions[2], interp_func,
                                                                    pb);
                        }
                    }
                } else {
                    for (size_t j = begin[0]; j <= end[0]; j += stride_ip2) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(j, begin[1], i), offset(j, end[1], i),
                                                                    stride_ip * global_dimensions[2], interp_func,
                                                                    pb);
                        }
                    }
                    for (size_t k = begin[1]; k <= end[1]; k += stride_ip) {
                        for (size_t i = begin[2]; i <= end[2]; i += stride_ip) {
                            predict_error += block_interpolation_1d(data, offset(begin[0], k, i), offset(end[0], k, i),
                                                                    stride_ip * global_dimensions[1] * global_dimensions[2],
                                                                    interp_func, pb);
                        }
                    }
                }
            }
            return predict_error;
        }

        int interpolator_id;
        int direction_sequence_id;
        std::vector<std::string> interpolators = {"linear", "cubic"};
        std::vector<int> quant_inds;
        size_t quant_index = 0; // for decompress
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        uint block_size;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

};


#endif

