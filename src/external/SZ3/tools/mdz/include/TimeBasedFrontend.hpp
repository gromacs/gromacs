#ifndef SZ_TIMEBASED_FRONTEND
#define SZ_TIMEBASED_FRONTEND

#include <SZ3/frontend/Frontend.hpp>
#include <SZ3/def.hpp>
#include <SZ3/predictor/Predictor.hpp>
#include <SZ3/predictor/LorenzoPredictor.hpp>
#include <SZ3/quantizer/Quantizer.hpp>
#include <SZ3/utils/Iterator.hpp>
#include <SZ3/utils/Config.hpp>
#include <SZ3/utils/MemoryUtil.hpp>

namespace SZ3 {


    template<class T, uint N, class Predictor, class Quantizer>
    class TimeBasedFrontend : public concepts::FrontendInterface<T, N> {
    public:

        TimeBasedFrontend(const Config &conf, Predictor predictor, Quantizer quantizer, T *data_ts0) :
                fallback_predictor(LorenzoPredictor<T, N - 1, 1>(conf.absErrorBound)),
                predictor(predictor),
                quantizer(quantizer),
                block_size(conf.blockSize),
                stride(conf.stride),
                num_elements(conf.num),
                data_ts0(data_ts0) {
            assert((conf.dims.size() == 2) &&
                   "timestep prediction requires 2d dataset");
            global_dimensions[0] = conf.dims[0];
            global_dimensions[1] = conf.dims[1];
        }

        ~TimeBasedFrontend() = default;

        std::vector<int> compress(T *data) {
            std::vector<int> quant_inds(num_elements);
            size_t quant_count = 0;
            if (data_ts0 != nullptr) {
                for (size_t j = 0; j < global_dimensions[1]; j++) {
                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(data[j], data_ts0[j]);
                }
            } else {
                std::array<size_t, N - 1> global_dims;
                for (int i = 0; i < N - 1; i++) {
                    global_dims[i] = global_dimensions[i + 1];
                };

                auto inter_block_range = std::make_shared<SZ3::multi_dimensional_range<T, N - 1>>(
                        data, std::begin(global_dims), std::end(global_dims), stride, 0);
                auto intra_block_range = std::make_shared<SZ3::multi_dimensional_range<T, N - 1>>(
                        data, std::begin(global_dims), std::end(global_dims), 1, 0);

//                std::array<size_t, N - 1> intra_block_dims;
                predictor.precompress_data(inter_block_range->begin());
                quantizer.precompress_data();
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    intra_block_range->update_block_range(block, block_size);

//                    // std::cout << *block << " " << lp.predict(block) << std::endl;
//                    for (int i = 0; i < intra_block_dims.size(); i++) {
//                        size_t cur_index = block.get_local_index(i);
//                        size_t dims = inter_block_range->get_dimensions(i);
//                        intra_block_dims[i] = (cur_index == dims - 1 &&
//                                               global_dims[i] - cur_index * stride < block_size) ?
//                                              global_dims[i] - cur_index * stride : block_size;
//                    }
//
//                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
//                    intra_block_range->set_offsets(block.get_offset());
//                    intra_block_range->set_starting_position(block.get_local_index());
                    concepts::PredictorInterface<T, N - 1> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    predictor_withfallback->precompress_block_commit();
//                    quantizer.precompress_block();
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                                *element, predictor_withfallback->predict(element));
                    }
                }
                predictor.postcompress_data(inter_block_range->begin());
            }
            for (size_t j = 0; j < global_dimensions[1]; j++) {
                for (size_t i = 1; i < global_dimensions[0]; i++) {
                    size_t idx = i * global_dimensions[1] + j;
                    size_t idx_prev = (i - 1) * global_dimensions[1] + j;
                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(data[idx], data[idx_prev]);
                }
            }
            assert(quant_count == num_elements);
            quantizer.postcompress_data();
            return quant_inds;
        }

        T *decompress(std::vector<int> &quant_inds, T *dec_data) {

            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::array<size_t, N - 1> intra_block_dims;
//            auto dec_data = new T[num_elements];

            if (data_ts0 != nullptr) {
                for (size_t j = 0; j < global_dimensions[1]; j++) {
                    dec_data[j] = quantizer.recover(data_ts0[j], *(quant_inds_pos++));
                }
            } else {
                std::array<size_t, N - 1> global_dims;
                for (int i = 0; i < N - 1; i++) {
                    global_dims[i] = global_dimensions[i + 1];
                };
                auto inter_block_range = std::make_shared<SZ3::multi_dimensional_range<T, N - 1>>(
                        dec_data, std::begin(global_dims), std::end(global_dims), stride, 0);

                auto intra_block_range = std::make_shared<SZ3::multi_dimensional_range<T, N - 1>>(
                        dec_data, std::begin(global_dims), std::end(global_dims), 1, 0);

                predictor.predecompress_data(inter_block_range->begin());
                quantizer.predecompress_data();

                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    intra_block_range->update_block_range(block, block_size);
//                    for (int i = 0; i < intra_block_dims.size(); i++) {
//                        size_t cur_index = block.get_local_index(i);
//                        size_t dims = inter_block_range->get_dimensions(i);
//                        intra_block_dims[i] = (cur_index == dims - 1) ? global_dims[i] - cur_index * block_size
//                                                                      : block_size;
//                    }
//                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
//                    intra_block_range->set_offsets(block.get_offset());
//                    intra_block_range->set_starting_position(block.get_local_index());

                    concepts::PredictorInterface<T, N - 1> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    auto intra_begin = intra_block_range->begin();
                    auto intra_end = intra_block_range->end();
                    for (auto element = intra_begin; element != intra_end; ++element) {
                        *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++));
                    }
                }
                predictor.postdecompress_data(inter_block_range->begin());
            }

            for (size_t j = 0; j < global_dimensions[1]; j++) {
                for (size_t i = 1; i < global_dimensions[0]; i++) {
                    size_t idx = i * global_dimensions[1] + j;
                    size_t idx_prev = (i - 1) * global_dimensions[1] + j;
                    dec_data[idx] = quantizer.recover(dec_data[idx_prev], *(quant_inds_pos++));
                }
            }

            quantizer.postdecompress_data();
            return dec_data;
        }

        void save(uchar *&c) {
            write(global_dimensions.data(), N, c);
            write(block_size, c);

            predictor.save(c);
            quantizer.save(c);
        }

        void load(const uchar *&c, size_t &remaining_length) {
            read(global_dimensions.data(), N, c, remaining_length);
            num_elements = 1;
            for (const auto &d: global_dimensions) {
                num_elements *= d;
//                std::cout << d << " ";
            }
//            std::cout << std::endl;
            read(block_size, c, remaining_length);
            stride = block_size;
            predictor.load(c, remaining_length);
            quantizer.load(c, remaining_length);
        }

        void print() {
//            predictor.print();
//            quantizer.print();
        }

        void clear() {
            predictor.clear();
            fallback_predictor.clear();
            quantizer.clear();
        }

        int get_radius() const { return quantizer.get_radius(); }

        size_t get_num_elements() const { return num_elements; };

        size_t size_est() { return 0; };
    private:
        Predictor predictor;
        LorenzoPredictor<T, N - 1, 1> fallback_predictor;
        Quantizer quantizer;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
        T *data_ts0 = nullptr;
    };

    template<class T, uint N, class Predictor, class Quantizer>
    TimeBasedFrontend<T, N, Predictor, Quantizer>
    make_sz3_timebased_frontend(const Config &conf, Predictor predictor, Quantizer quantizer, T *data_ts0) {
        return TimeBasedFrontend<T, N, Predictor, Quantizer>(conf, predictor, quantizer, data_ts0);
    }
}

#endif
