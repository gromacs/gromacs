#ifndef SZ3_BLOCKWISE_DECOMPOSITION_HPP
#define SZ3_BLOCKWISE_DECOMPOSITION_HPP

#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/predictor/LorenzoPredictor.hpp"
#include "SZ3/predictor/Predictor.hpp"
#include "SZ3/quantizer/LinearQuantizer.hpp"
#include "SZ3/utils/BlockwiseIterator.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {
template <class T, uint N, class Predictor, class Quantizer>
class BlockwiseDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    using Block_iter = typename block_data<T, N>::block_iterator;

    BlockwiseDecomposition(const Config &conf, Predictor predictor_, Quantizer quantizer_)
        : predictor(predictor_), quantizer(quantizer_), fallback_predictor(conf.absErrorBound) {
        static_assert(std::is_base_of<concepts::PredictorInterface<T, N>, Predictor>::value,
                      "must implement the Predictor interface");
        // A block the predictor declines falls back to Lorenzo, which reads neighbours the predictor may not pad for.
        padding = std::max(predictor.get_padding(), fallback_predictor.get_padding());
    }

    std::vector<int> compress(const Config &conf, T *data) override {
        auto data_with_padding = std::make_shared<block_data<T, N>>(data, conf.dims, padding, true);
        auto block = data_with_padding->block_iter(conf.blockSize);
        std::vector<int> quant_inds;
        quant_inds.reserve(conf.num);
        do {
            concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
            if (!predictor.precompress(block)) {
                predictor_withfallback = &fallback_predictor;
            }
            predictor_withfallback->precompress_block_commit();
            Block_iter::foreach (block, [&](T *c, const std::array<size_t, N> &index) {
                T pred = predictor_withfallback->predict(block, c, index);
                quant_inds.push_back(quantizer.quantize_and_overwrite(*c, pred));
            });

        } while (block.next());
        return quant_inds;
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        if (quant_inds.size() < conf.num) {
            throw std::out_of_range("SZ3 blockwise: fewer bins than the grid consumes");
        }
        int *quant_inds_pos = &quant_inds[0];

        auto data_with_padding = std::make_shared<block_data<T, N>>(dec_data, conf.dims, padding, false);
        auto block = data_with_padding->block_iter(conf.blockSize);
        do {
            concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
            if (!predictor.predecompress(block)) {
                predictor_withfallback = &fallback_predictor;
            }
            Block_iter::foreach (block, [&](T *c, const std::array<size_t, N> &index) {
                T pred = predictor_withfallback->predict(block, c, index);
                *c = quantizer.recover(pred, *(quant_inds_pos++));
            });

        } while (block.next());

        return dec_data;
    }

    size_t size_est() override { return fallback_predictor.size_est() + predictor.size_est() + quantizer.size_est(); }

    void save(uchar *&c) override {
        fallback_predictor.save(c);
        predictor.save(c);
        quantizer.save(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        fallback_predictor.load(c, remaining_length);
        predictor.load(c, remaining_length);
        quantizer.load(c, remaining_length);
    }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    Predictor predictor;
    Quantizer quantizer;
    LorenzoPredictor<T, N, 1> fallback_predictor;
    size_t padding;
};

template <class T, uint N, class Predictor, class Quantizer>
BlockwiseDecomposition<T, N, Predictor, Quantizer> make_decomposition_blockwise(const Config &conf, Predictor predictor,
                                                                                Quantizer quantizer) {
    return BlockwiseDecomposition<T, N, Predictor, Quantizer>(conf, predictor, quantizer);
}

}  // namespace SZ3
#endif
