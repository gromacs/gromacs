#ifndef SZ3_BYPASS_ENCODER_HPP
#define SZ3_BYPASS_ENCODER_HPP

#include <cassert>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "Encoder.hpp"
#include "SZ3/def.hpp"

namespace SZ3 {

template <class T>
class BypassEncoder : public concepts::EncoderInterface<T> {
   public:
    void preprocess_encode(const std::vector<T> &bins, int stateNum) override {
    }

    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        memcpy(bytes, &bins[0], sizeof(T) * bins.size());
        bytes += sizeof(T) * bins.size();
        return 0;
    }

    void postprocess_encode() override {}

    void preprocess_decode() override {}

    std::vector<T> decode(const uchar *&bytes, size_t targetLength, size_t &remaining_length) override {
        if (targetLength > remaining_length / sizeof(T)) {
            throw std::out_of_range("SZ3 bypass encoder: more bins requested than the buffer holds");
        }
        std::vector<T> bins(targetLength);
        memcpy(bins.data(), bytes, sizeof(T) * targetLength);
        bytes += sizeof(T) * targetLength;
        remaining_length -= sizeof(T) * targetLength;
        return bins;
    }

    void postprocess_decode() override {}

    void save(uchar *&c) override {}

    void load(const uchar *&c, size_t &remaining_length) override {}
};
}  // namespace SZ3
#endif
