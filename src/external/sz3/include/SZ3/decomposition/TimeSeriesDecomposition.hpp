#ifndef SZ3_TIME_SERIES_DECOMPOSITION_HPP
#define SZ3_TIME_SERIES_DECOMPOSITION_HPP

#include <cassert>
#include <limits>
#include <stdexcept>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/BlockwiseIterator.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {
template <class T, uint N, class Predictor, class Quantizer>
class TimeSeriesDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    using Block_iter = typename block_data<T, N - 1>::block_iterator;

    TimeSeriesDecomposition(const Config& conf, Predictor predictor, Quantizer quantizer, T* data_ts0)
        : predictor(predictor),
          fallback_predictor(LorenzoPredictor<T, N - 1, 1>(conf.absErrorBound)),
          quantizer(quantizer),
          num_elements(conf.num),
          data_ts0(data_ts0) {
        static_assert(std::is_base_of<concepts::PredictorInterface<T, N - 1>, Predictor>::value,
                      "must implement the predictor interface");
        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quantizer interface");
        assert((conf.dims.size() == 2) && "timestep prediction requires 2d dataset");
    }

    ~TimeSeriesDecomposition() = default;

    std::vector<int> compress(const Config& conf, T* data) override {
        std::vector<int> quant_inds(num_elements);
        size_t quant_count = 0;
        // The timestep loop below predicts from the reconstruction of timestep 0.
        const T* ts0_recon = data;
        std::shared_ptr<block_data<T, N - 1>> data_with_padding;
        if (data_ts0 != nullptr) {
            for (size_t j = 0; j < conf.dims[1]; j++) {
                quant_inds[quant_count++] = quantizer.quantize_and_overwrite(data[j], data_ts0[j]);
            }
        } else {
            std::vector<size_t> spatial_dims(N - 1);
            for (uint i = 0; i < N - 1; i++) {
                spatial_dims[i] = conf.dims[i + 1];
            };

            data_with_padding =
                std::make_shared<block_data<T, N - 1>>(data, spatial_dims, predictor.get_padding(), true);
            auto block = data_with_padding->block_iter(conf.blockSize);
            do {
                concepts::PredictorInterface<T, N - 1>* predictor_withfallback = &predictor;
                if (!predictor.precompress(block)) {
                    predictor_withfallback = &fallback_predictor;
                }
                predictor_withfallback->precompress_block_commit();
                Block_iter::foreach (block, [&](T* c, const std::array<size_t, N - 1>& index) {
                    T pred = predictor_withfallback->predict(block, c, index);
                    quant_inds[quant_count++] = quantizer.quantize_and_overwrite(*c, pred);
                });
            } while (block.next());

            ts0_recon = data_with_padding->values();
        }

        for (size_t j = 0; j < conf.dims[1]; j++) {
            T prev = ts0_recon[j];
            for (size_t i = 1; i < conf.dims[0]; i++) {
                size_t idx = i * conf.dims[1] + j;
                quant_inds[quant_count++] = quantizer.quantize_and_overwrite(data[idx], prev);
                prev = data[idx];  // quantize_and_overwrite left the reconstruction here
            }
        }
        assert(quant_count == num_elements);
        quantizer.postcompress_data();
        return quant_inds;
    }

    T* decompress(const Config& conf, std::vector<int>& quant_inds, T* dec_data) override {
        if (quant_inds.size() < num_elements) {
            throw std::out_of_range("SZ3 time series: fewer bins than the grid consumes");
        }
        int const* quant_inds_pos = quant_inds.data();
        // std::array<size_t, N - 1> intra_block_dims;
        //            auto dec_data = new T[num_elements];

        if (data_ts0 != nullptr) {
            for (size_t j = 0; j < conf.dims[1]; j++) {
                dec_data[j] = quantizer.recover(data_ts0[j], *(quant_inds_pos++));
            }
        } else {
            std::vector<size_t> spatial_dims(N - 1);
            for (uint i = 0; i < N - 1; i++) {
                spatial_dims[i] = conf.dims[i + 1];
            };

            auto data_with_padding =
                std::make_shared<block_data<T, N - 1>>(dec_data, spatial_dims, predictor.get_padding(), false);
            auto block = data_with_padding->block_iter(conf.blockSize);
            do {
                concepts::PredictorInterface<T, N - 1>* predictor_withfallback = &predictor;
                if (!predictor.predecompress(block)) {
                    predictor_withfallback = &fallback_predictor;
                }
                Block_iter::foreach (block, [&](T* c, const std::array<size_t, N - 1>& index) {
                    T pred = predictor_withfallback->predict(block, c, index);
                    *c = quantizer.recover(pred, *(quant_inds_pos++));
                });
            } while (block.next());
        }

        for (size_t j = 0; j < conf.dims[1]; j++) {
            for (size_t i = 1; i < conf.dims[0]; i++) {
                size_t idx = i * conf.dims[1] + j;
                size_t idx_prev = (i - 1) * conf.dims[1] + j;
                dec_data[idx] = quantizer.recover(dec_data[idx_prev], *(quant_inds_pos++));
            }
        }

        quantizer.postdecompress_data();
        return dec_data;
    }

    size_t size_est() override { return fallback_predictor.size_est() + predictor.size_est() + quantizer.size_est(); }

    void save(uchar*& c) override {
        fallback_predictor.save(c);
        predictor.save(c);
        quantizer.save(c);
    }

    void load(const uchar*& c, size_t& remaining_length) override {
        fallback_predictor.load(c, remaining_length);
        predictor.load(c, remaining_length);
        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    Predictor predictor;
    LorenzoPredictor<T, N - 1, 1> fallback_predictor;
    Quantizer quantizer;
    size_t num_elements;
    T* data_ts0 = nullptr;
};

template <class T, uint N, class Predictor, class Quantizer>
TimeSeriesDecomposition<T, N, Predictor, Quantizer> make_decomposition_timeseries(const Config& conf,
                                                                                  Predictor predictor,
                                                                                  Quantizer quantizer, T* data_ts0) {
    return TimeSeriesDecomposition<T, N, Predictor, Quantizer>(conf, predictor, quantizer, data_ts0);
}
}  // namespace SZ3

#endif