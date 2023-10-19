#ifndef SZ_COMPRESSOR_HPP
#define SZ_COMPRESSOR_HPP

#include "SZ3/def.hpp"
#include "SZ3/utils/Config.hpp"

namespace SZ {
    namespace concepts {
        template<class T>
        class CompressorInterface {
        public:
            virtual T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) = 0;

            virtual T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) = 0;

            virtual uchar *compress(const Config &conf, T *data, size_t &compressed_size) = 0;
        };
    }
}
#endif
