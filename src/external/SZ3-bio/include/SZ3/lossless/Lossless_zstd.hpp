//
// Created by Kai Zhao on 4/21/20.
//

#ifndef SZ_LOSSLESS_ZSTD_HPP
#define SZ_LOSSLESS_ZSTD_HPP

#include "zstd.h"
#include "SZ3/def.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/lossless/Lossless.hpp"

namespace SZ3 {
    class Lossless_zstd : public concepts::LosslessInterface {

    public:
        Lossless_zstd() = default;

        Lossless_zstd(int comp_level) : compression_level(comp_level) {};

        uchar *compress(uchar *data, size_t dataLength, size_t &outSize) {
            size_t estimatedCompressedSize = std::max(size_t(dataLength * 1.2), size_t(400));
            uchar *compressBytes = new uchar[estimatedCompressedSize];
            uchar *compressBytesPos = compressBytes;
            write(dataLength, compressBytesPos);

            outSize = ZSTD_compress(compressBytesPos, estimatedCompressedSize, data, dataLength,
                                    compression_level);
            outSize += sizeof(size_t);
            return compressBytes;
        }

        void postcompress_data(uchar *data) {
            delete[] data;
        }

        uchar *decompress(const uchar *data, size_t &compressedSize) {
            const uchar *dataPos = data;
            size_t dataLength = 0;
            read(dataLength, dataPos, compressedSize);

            uchar *oriData = new uchar[dataLength];
            ZSTD_decompress(oriData, dataLength, dataPos, compressedSize);
            compressedSize = dataLength;
            return oriData;
        }

        void postdecompress_data(uchar *data) {
            delete[] data;
        }

    private:
        int compression_level = 3;  //default setting of level is 3
    };
}
#endif //SZ_LOSSLESS_ZSTD_HPP
