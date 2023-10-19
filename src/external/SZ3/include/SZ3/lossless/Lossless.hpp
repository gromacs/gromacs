//
// Created by Kai Zhao on 6/12/20.
//

#ifndef SZ_LOSSLESS_HPP
#define SZ_LOSSLESS_HPP


namespace SZ {
    namespace concepts {

        class LosslessInterface {
        public:
            virtual void postcompress_data(uchar *data) = 0;

            virtual void postdecompress_data(uchar *data) = 0;

            virtual uchar *compress(uchar *data, size_t dataLength, size_t &outSize) = 0;

            virtual uchar *decompress(const uchar *data, size_t& compressedSize) = 0;
        };
    }
}

#endif //SZ_LOSSLESS_HPP
