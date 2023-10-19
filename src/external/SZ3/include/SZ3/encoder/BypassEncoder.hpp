#ifndef _SZ_BYPASS_ENCODER_HPP
#define _SZ_BYPASS_ENCODER_HPP

#include "Encoder.hpp"
#include "SZ3/def.hpp"
#include <vector>

namespace SZ {

    template<class T>
    class BypassEncoder : public concepts::EncoderInterface<T> {
    public:

        ~BypassEncoder() = default;

        void preprocess_encode(const std::vector<T> &bins, int stateNum) {
            assert(stateNum <= 256 && "stateNum should be no more than 256.");
        };

        size_t encode(const std::vector<T> &bins, uchar *&bytes) {
            for (auto &bin: bins) {
                *bytes++ = uchar(bin);
            }
            return 0;
        };

        void postprocess_encode() {};

        void preprocess_decode() {};

        std::vector<T> decode(const uchar *&bytes, size_t targetLength) {
            std::vector<T> bins(targetLength);
            for (auto &bin: bins) {
                bin = *bytes++;
            }
            return bins;
        };

        void postprocess_decode() {};

        uint save(uchar *&c) {
            return 0;
        };

        void load(const uchar *&c, size_t &remaining_length) {};

    };
}
#endif
