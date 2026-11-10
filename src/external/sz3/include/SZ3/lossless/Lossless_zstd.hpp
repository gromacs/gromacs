//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ3_LOSSLESS_ZSTD_HPP
#define SZ3_LOSSLESS_ZSTD_HPP

#include <memory>
#include <stdexcept>

#include "SZ3/def.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/MemoryUtil.hpp"

// Do not include <zstd.h> here: these headers are installed, so every consumer would need it too.
#ifdef SZ3_USE_ZSTD_HEADER
#include <zstd.h>
#else
#include <stddef.h>
extern "C" {
size_t ZSTD_compress(void *dst, size_t dstCapacity, const void *src, size_t srcSize, int compressionLevel);
size_t ZSTD_decompress(void *dst, size_t dstCapacity, const void *src, size_t compressedSize);
size_t ZSTD_compressBound(size_t srcSize);
unsigned ZSTD_isError(size_t code);
}
#endif

namespace SZ3 {
class Lossless_zstd : public concepts::LosslessInterface {
   public:
    Lossless_zstd() = default;

    Lossless_zstd(int comp_level) : compression_level(comp_level) {}

    /** ZSTD_compressBound is declared only in this header, so callers elsewhere go through this. */
    static size_t compress_bound(size_t srcLen) { return ZSTD_compressBound(srcLen); }

    /**
     * compress data with lossless compressors
     * @param src  data to be compressed
     * @param srcLen length (in bytes) of the data to be compressed
     * @param dst compressed data
     * @param dstCap capacity (in bytes) for storing the compressed data
     * @return length (in bytes) of the data compressed
     */
    size_t compress(const uchar *src, size_t srcLen, uchar *dst, size_t dstCap) override {
        if (dstCap < sizeof(uint64_t)) throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        write<uint64_t>(srcLen, dst);
        dstCap -= sizeof(uint64_t);  // reserve space for srcLen
        size_t dstLen = ZSTD_compress(dst, dstCap, src, srcLen, compression_level);
        if (ZSTD_isError(dstLen)) {
            throw std::length_error(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
        }
        return dstLen + sizeof(uint64_t);
    }

    /**
     * reverse of compress(), decompress the data with lossless compressors
     * @param src data to be decompressed
     * @param srcLen length (in bytes) of that data
     * @param dst buffer to decompress into; when null on entry the callee allocates it with malloc()
     *            and the caller frees it
     * @param dstCap the capacity of dst, ignored when dst is null
     * @return length (in bytes) of the data decompressed
     */
    size_t decompress(const uchar *src, size_t srcLen, uchar *&dst, size_t dstCap) override {
        // The stream is a decompressed-size field followed by the zstd frame, all untrusted.
        if (srcLen < sizeof(uint64_t)) {
            throw std::out_of_range("SZ3 lossless: compressed data is smaller than the size header");
        }
        uint64_t dstLen = 0;
        read(dstLen, src);

        // malloc, because the caller frees what it gets back with free().
        std::unique_ptr<uchar, void (*)(void *)> owner(nullptr, &free);
        if (dst == nullptr) {
            owner.reset(static_cast<uchar *>(malloc(dstLen)));
            if (owner == nullptr) {
                throw std::runtime_error("SZ3 lossless: can not allocate the decompression buffer");
            }
        } else if (dstLen > dstCap) {
            throw std::out_of_range("SZ3 lossless: declared decompressed size exceeds the allowed capacity");
        }
        uchar *out = (dst != nullptr) ? dst : owner.get();

        // A short frame would leave the tail of the output uninitialized for the caller to read.
        size_t res = ZSTD_decompress(out, dstLen, src, srcLen - sizeof(uint64_t));
        if (ZSTD_isError(res) || res != dstLen) {
            throw std::runtime_error("SZ3 lossless: stream does not decompress to the size it declares");
        }

        dst = out;
        owner.release();
        return res;
    }

   private:
    int compression_level = 3;  // default setting of level is 3
};
}  // namespace SZ3
#endif  // SZ_LOSSLESS_ZSTD_HPP
