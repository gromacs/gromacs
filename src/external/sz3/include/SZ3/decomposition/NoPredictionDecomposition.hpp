#ifndef SZ3_NO_PREDICTION_DECOMPOSITION_HPP
#define SZ3_NO_PREDICTION_DECOMPOSITION_HPP

#include <stdexcept>

#include "Decomposition.hpp"
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ3 {
template <class T, uint N, class Quantizer>
class NoPredictionDecomposition : public concepts::DecompositionInterface<T, int, N> {
   public:
    NoPredictionDecomposition(const Config & /*conf*/, Quantizer quantizer_) : quantizer(quantizer_) {
        static_assert(std::is_base_of<concepts::QuantizerInterface<T, int>, Quantizer>::value,
                      "must implement the quantizer interface");
    }

    T *decompress(const Config &conf, std::vector<int> &quant_inds, T *dec_data) override {
        if (quant_inds.size() < conf.num) {
            throw std::out_of_range("SZ3 no-prediction: fewer bins than the data has elements");
        }
        for (size_t i = 0; i < conf.num; i++) {
            dec_data[i] = quantizer.recover(0, quant_inds[i]);
        }
        quantizer.postdecompress_data();
        return dec_data;
    }

    std::vector<int> compress(const Config &conf, T *data) override {
        std::vector<int> quant_inds(conf.num);
        for (size_t i = 0; i < conf.num; i++) {
            quant_inds[i] = quantizer.quantize_and_overwrite(data[i], 0);
        }
        quantizer.postcompress_data();
        return quant_inds;
    }

    size_t size_est() override { return quantizer.size_est() + 64; }

    void save(uchar *&c) override { quantizer.save(c); }

    void load(const uchar *&c, size_t &remaining_length) override { quantizer.load(c, remaining_length); }

    std::pair<int, int> get_out_range() override { return quantizer.get_out_range(); }

   private:
    Quantizer quantizer;
};

template <class T, uint N, class Quantizer>
NoPredictionDecomposition<T, N, Quantizer> make_decomposition_noprediction(const Config &conf, Quantizer quantizer) {
    return NoPredictionDecomposition<T, N, Quantizer>(conf, quantizer);
}

}  // namespace SZ3

#endif
