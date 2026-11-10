#ifndef SZ3_LORENZO_PREDICTOR_HPP
#define SZ3_LORENZO_PREDICTOR_HPP

#include <iostream>
#include <limits>
#include <type_traits>

#include "SZ3/predictor/Predictor.hpp"

namespace SZ3 {

// N-dimension L-layer lorenzo predictor
template <class T, uint N, uint L>
class LorenzoPredictor : public concepts::PredictorInterface<T, N> {
    static_assert(std::is_floating_point<T>::value, "LorenzoPredictor only takes a floating-point T");

   public:
    static const uint8_t predictor_id = 0b00000001;
    using block_iter = typename block_data<T, N>::block_iterator;

    LorenzoPredictor() { this->noise = 0; }

    LorenzoPredictor(double eb) {
        this->noise = 0;
        if (L == 1) {
            if (N == 1) {
                this->noise = 0.5 * eb;
            } else if (N == 2) {
                this->noise = 0.81 * eb;
            } else if (N == 3) {
                this->noise = 1.22 * eb;
            } else if (N == 4) {
                this->noise = 1.79 * eb;
            }
        } else if (L == 2) {
            if (N == 1) {
                this->noise = 1.08 * eb;
            } else if (N == 2) {
                this->noise = 2.76 * eb;
            } else if (N == 3) {
                this->noise = 6.8 * eb;
            }
        }
    }

    bool precompress(const block_iter &) override { return true; }

    void precompress_block_commit() noexcept override {}

    bool predecompress(const block_iter &) override { return true; }

    void save(uchar *&) override {}

    void load(const uchar *&, size_t &) override {}

    void print() const override {
        std::cout << L << "-Layer " << N << "D Lorenzo predictor, noise = " << noise << "\n";
    }

    size_t get_padding() override { return 2; }

    T estimate_error(const block_iter &block, T *d, const std::array<size_t, N> &index) override {
        return fabs(*d - predict(block, d, index)) + this->noise;
    }

    T predict(const block_iter &block, T *d, const std::array<size_t, N> & /*index*/) override {
        auto ds = block.get_dim_strides();
        if constexpr (N == 1 && L == 1) {
            return prev1(d, 1);
        } else if constexpr (N == 2 && L == 1) {
            return prev2(d, ds, 0, 1) + prev2(d, ds, 1, 0) - prev2(d, ds, 1, 1);
        } else if constexpr (N == 3 && L == 1) {
            return prev3(d, ds, 0, 0, 1) + prev3(d, ds, 0, 1, 0) + prev3(d, ds, 1, 0, 0) - prev3(d, ds, 0, 1, 1) -
                   prev3(d, ds, 1, 0, 1) - prev3(d, ds, 1, 1, 0) + prev3(d, ds, 1, 1, 1);
        } else if constexpr (N == 4 && L == 1) {
            return prev4(d, ds, 0, 0, 0, 1) + prev4(d, ds, 0, 0, 1, 0) - prev4(d, ds, 0, 0, 1, 1) +
                   prev4(d, ds, 0, 1, 0, 0) - prev4(d, ds, 0, 1, 0, 1) - prev4(d, ds, 0, 1, 1, 0) +
                   prev4(d, ds, 0, 1, 1, 1) + prev4(d, ds, 1, 0, 0, 0) - prev4(d, ds, 1, 0, 0, 1) -
                   prev4(d, ds, 1, 0, 1, 0) + prev4(d, ds, 1, 0, 1, 1) - prev4(d, ds, 1, 1, 0, 0) +
                   prev4(d, ds, 1, 1, 0, 1) + prev4(d, ds, 1, 1, 1, 0) - prev4(d, ds, 1, 1, 1, 1);
        } else if constexpr (N == 1 && L == 2) {
            return 2 * prev1(d, 1) - prev1(d, 2);
        } else if constexpr (N == 2 && L == 2) {
            return 2 * prev2(d, ds, 0, 1) - prev2(d, ds, 0, 2) + 2 * prev2(d, ds, 1, 0) - 4 * prev2(d, ds, 1, 1) +
                   2 * prev2(d, ds, 1, 2) - prev2(d, ds, 2, 0) + 2 * prev2(d, ds, 2, 1) - prev2(d, ds, 2, 2);
        } else if constexpr (N == 3 && L == 2) {
            return 2 * prev3(d, ds, 0, 0, 1) - prev3(d, ds, 0, 0, 2) + 2 * prev3(d, ds, 0, 1, 0) -
                   4 * prev3(d, ds, 0, 1, 1) + 2 * prev3(d, ds, 0, 1, 2) - prev3(d, ds, 0, 2, 0) +
                   2 * prev3(d, ds, 0, 2, 1) - prev3(d, ds, 0, 2, 2) + 2 * prev3(d, ds, 1, 0, 0) -
                   4 * prev3(d, ds, 1, 0, 1) + 2 * prev3(d, ds, 1, 0, 2) - 4 * prev3(d, ds, 1, 1, 0) +
                   8 * prev3(d, ds, 1, 1, 1) - 4 * prev3(d, ds, 1, 1, 2) + 2 * prev3(d, ds, 1, 2, 0) -
                   4 * prev3(d, ds, 1, 2, 1) + 2 * prev3(d, ds, 1, 2, 2) - prev3(d, ds, 2, 0, 0) +
                   2 * prev3(d, ds, 2, 0, 1) - prev3(d, ds, 2, 0, 2) + 2 * prev3(d, ds, 2, 1, 0) -
                   4 * prev3(d, ds, 2, 1, 1) + 2 * prev3(d, ds, 2, 1, 2) - prev3(d, ds, 2, 2, 0) +
                   2 * prev3(d, ds, 2, 2, 1) - prev3(d, ds, 2, 2, 2);
        } else {
            static_assert(N <= 4 && L <= 2, "Unsupported dimension or layer configuration");
            return T(0);
        }
    }

   protected:
    T noise = 0;

   private:
    // Helper functions for Lorenzo prediction. The neighbour offsets are unsigned, so `d - offset` keeps
    // the step negative instead of wrapping it into an out-of-bounds pointer.
    T prev1(T *d, size_t i) { return *(d - i); }
    T prev2(T *d, const std::array<size_t, N> &ds, size_t j, size_t i) { return *(d - (j * ds[0] + i)); }
    T prev3(T *d, const std::array<size_t, N> &ds, size_t k, size_t j, size_t i) {
        return *(d - (k * ds[1] + j * ds[0] + i));
    }
    T prev4(T *d, const std::array<size_t, N> &ds, size_t t, size_t k, size_t j, size_t i) {
        return *(d - (t * ds[2] + k * ds[1] + j * ds[0] + i));
    }
};
}  // namespace SZ3
#endif
