#ifndef SZ_Truncate_COMPRESSOR_HPP
#define SZ_Truncate_COMPRESSOR_HPP

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/frontend/Frontend.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/def.hpp"
#include <cstring>

namespace SZ3 {
    template<class T, uint N, class Lossless>
    class SZTruncateCompressor : public concepts::CompressorInterface<T> {
    public:


        SZTruncateCompressor(const Config &conf, Lossless lossless, int byteLens) :
                lossless(lossless), conf(conf), byteLen(byteLens) {
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {

            auto compressed_data = new uchar[conf.num * sizeof(T)];
            auto compressed_data_pos = (uchar *) compressed_data;

            Timer timer(true);
            truncateArray(data, conf.num, byteLen, compressed_data_pos);
            timer.stop("Prediction & Quantization");

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     (uchar *) compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t& cmpSize, T *decData) {
            size_t remaining_length = cmpSize;

            auto compressed_data = lossless.decompress(cmpData, remaining_length);
            auto compressed_data_pos = (uchar *) compressed_data;

            Timer timer(true);
//            auto dec_data = new T[conf.num];
            truncateArrayRecover(compressed_data_pos, conf.num, byteLen, decData);

            lossless.postdecompress_data(compressed_data);
            timer.stop("Prediction & Recover");
            return decData;
        }



    private:
        Lossless lossless;
        Config conf;
        int byteLen = 2;
    };

    template<class T, uint N, class Lossless>
    SZTruncateCompressor<T, N, Lossless>
    make_sz_truncate_compressor(const Config &conf, Lossless lossless, int byteLens) {
        return SZTruncateCompressor<T, N, Lossless>(conf, lossless, byteLens);
    }
}
#endif
