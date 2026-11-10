#ifndef SZ3_SEGMENTED_ENCODER_HPP
#define SZ3_SEGMENTED_ENCODER_HPP

#include <limits>
#include <stdexcept>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

namespace SZ3 {

/**
 * Bins that are a sequence of segments [n, n bins]: each segment gets its own Encoder, so a decomposition that
 * concatenates streams of different statistics gets a code per stream. Decoding returns the bins as they were.
 *
 * save() writes the segment count as a uint32, then for each segment its size as a LEB128 varint (at most that of
 * INT_MAX) and, unless it is empty, its encoder's save(); encode() writes the encoders' encode() in segment order.
 * A stream with more segments than max_segments is refused on load(), before an encoder is built for each.
 */
template <class Encoder>
class SegmentedEncoder : public concepts::EncoderInterface<int> {
   public:
    explicit SegmentedEncoder(uint32_t max_segments) : max_segments_(max_segments) {}

    void preprocess_encode(const std::vector<int> &bins, int stateNum) override {
        segments_.clear();
        encoders_.clear();
        for (size_t k = 0; k < bins.size();) {
            const size_t n = size_t(bins[k]);
            if (bins[k] < 0 || n > bins.size() - k - 1)
                throw std::invalid_argument(
                    "SZ3 SegmentedEncoder: bad input format, a segment length is negative or runs past the end");
            segments_.emplace_back(bins.begin() + k + 1, bins.begin() + k + 1 + n);
            if (n) {
                encoders_.emplace_back();
                encoders_.back().preprocess_encode(segments_.back(), stateNum);
            }
            k += n + 1;
        }
    }

    size_t encode(const std::vector<int> & /*bins*/, uchar *&bytes) override {
        uchar *const start = bytes;
        for (size_t s = 0, e = 0; s < segments_.size(); s++)
            if (!segments_[s].empty()) encoders_[e++].encode(segments_[s], bytes);
        return bytes - start;
    }

    void postprocess_encode() override {
        for (auto &encoder : encoders_) encoder.postprocess_encode();
        segments_.clear();
    }

    void preprocess_decode() override {
        for (auto &encoder : encoders_) encoder.preprocess_decode();
    }

    std::vector<int> decode(const uchar *&bytes, size_t targetLength, size_t &remaining_length) override {
        // a count and the bins of each segment; checked before anything is allocated
        uint64_t total = segment_sizes_.size();
        for (const uint64_t n : segment_sizes_) {
            if (total > targetLength || n > targetLength - total)
                throw std::out_of_range(
                    "SZ3 SegmentedEncoder: corrupt stream, the segment lengths do not add up to the bin count");
            total += n;
        }
        if (total != targetLength)
            throw std::out_of_range(
                "SZ3 SegmentedEncoder: corrupt stream, the segment lengths do not add up to the bin count");
        std::vector<int> out;
        out.reserve(targetLength);
        for (size_t s = 0, e = 0; s < segment_sizes_.size(); s++) {
            out.push_back(int(segment_sizes_[s]));
            if (!segment_sizes_[s]) continue;
            const auto v = encoders_[e++].decode(bytes, segment_sizes_[s], remaining_length);
            out.insert(out.end(), v.begin(), v.end());
        }
        return out;
    }

    void postprocess_decode() override {
        for (auto &encoder : encoders_) encoder.postprocess_decode();
    }

    void save(uchar *&c) override {
        write(uint32_t(segments_.size()), c);
        for (size_t s = 0, e = 0; s < segments_.size(); s++) {
            write_varint(segments_[s].size(), c);
            if (!segments_[s].empty()) encoders_[e++].save(c);
        }
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        uint32_t count = 0;
        read(count, c, remaining_length);
        if (count > remaining_length)
            throw std::out_of_range("SZ3 SegmentedEncoder: corrupt stream, more segments than bytes");
        if (count > max_segments_)
            throw std::out_of_range("SZ3 SegmentedEncoder: corrupt stream, more segments than the data has");
        // an encoder per non-empty segment as it is read, so memory follows the stream rather than the count it claims
        segment_sizes_.clear();
        encoders_.clear();
        for (uint32_t s = 0; s < count; s++) {
            const uint64_t size = read_varint(c, remaining_length);
            if (size > uint64_t(std::numeric_limits<int>::max()))
                throw std::out_of_range("SZ3 SegmentedEncoder: corrupt stream, a segment length beyond INT_MAX");
            segment_sizes_.push_back(size);
            if (size) {
                encoders_.emplace_back();
                encoders_.back().load(c, remaining_length);
            }
        }
    }

    size_t size_est() override {
        size_t bound = sizeof(uint32_t) + 5 * segments_.size();  // 5: a varint of 31 bits
        for (auto &encoder : encoders_) bound += encoder.size_est();
        return bound;
    }

   private:
    std::vector<Encoder> encoders_;  // one per non-empty segment
    std::vector<std::vector<int>> segments_;
    std::vector<uint64_t> segment_sizes_;
    uint32_t max_segments_;
};

}  // namespace SZ3
#endif
