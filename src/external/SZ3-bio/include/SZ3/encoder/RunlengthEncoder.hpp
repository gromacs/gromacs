#ifndef _SZ_RUNLENGTH_ENCODER_HPP
#define _SZ_RUNLENGTH_ENCODER_HPP

#include "Encoder.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/def.hpp"
#include <vector>

namespace SZ3 {

    template<class T>
    class RunlengthEncoder : public concepts::EncoderInterface<T> {
    public:

        ~RunlengthEncoder() = default;

        void preprocess_encode(const std::vector<T> &bins, int stateNum) {
        };

        size_t encode(const std::vector<T> &bins, uchar *&bytes) {
            int max = 0;
            size_t s = 0;
            for (size_t i = 1; i < bins.size(); i++) {
                if (bins[i] != bins[i - 1]) {
                    write(bins[i - 1], bytes);
                    write(int(i - s), bytes);
                    if (int(i - s) > max) {
                        max = int(i - s);
                    }
                    s = i;
                }
            }
            write(bins[bins.size() - 1], bytes);
            write(int(bins.size() - s), bytes);
            return 0;
        };

        void postprocess_encode() {};

        void preprocess_decode() {};

        std::vector<T> decode(const uchar *&bytes, size_t targetLength) {
            std::vector<T> bins(targetLength, 0);
            T value;
            int cnt;
            for (size_t i = 0; i < bins.size();) {
                read(value, bytes);
                read(cnt, bytes);
                for (size_t j = i; j < i + cnt; j++) {
                    bins[j] = value;
                }
                i += cnt;
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
