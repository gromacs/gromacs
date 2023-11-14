//
// Created by Kai Zhao on 6/12/20.
//

#ifndef SZ_LOSSLESS_HPP
#define SZ_LOSSLESS_HPP


namespace SZ3 {
    namespace concepts {

        class LosslessInterface {
        public:
            virtual void postcompress_data(uchar *data) = 0;

            virtual void postdecompress_data(uchar *data) = 0;

            /**
             * compress data with lossless compressors
             * @param data data to be compressed
             * @param dataLength length (in bytes) of the data to be compressed
             * @param outSize compressed size (in bytes)
             * @return compressed data
             */
            virtual uchar *compress(uchar *data, size_t dataLength, size_t &outSize) = 0;

            /**
             * reverse of compress(), decompress the data with lossless compressors
             * @param data data to be decompressed
             * @param compressedSize length (in bytes) of the data to be decompressed (as input) or the data decompressed (as output).
             * @return decompressed data
             */
            virtual uchar *decompress(const uchar *data, size_t& compressedSize) = 0;
        };
    }
}

#endif //SZ_LOSSLESS_HPP
