#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Statistic.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/api/impl/SZInterp.hpp"
#include "SZ3/api/impl/SZLorenzoReg.hpp"
#include <cmath>

namespace SZ3 {
    template<class T, uint N>
    char *SZ_compress_dispatcher(Config &conf, T *data, size_t &outSize) {

        assert(N == conf.N);
        calAbsErrorBound(conf, data);

        char *cmpData;
        if (conf.absErrorBound == 0) {
            auto zstd = Lossless_zstd();
            cmpData = (char *) zstd.compress((uchar *) data, conf.num * sizeof(T), outSize);
        } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
            cmpData = (char *) SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
        } else if (conf.cmprAlgo == ALGO_INTERP) {
            cmpData = (char *) SZ_compress_Interp<T, N>(conf, data, outSize);
        } else if (conf.cmprAlgo == ALGO_INTERP_LORENZO) {
            cmpData = (char *) SZ_compress_Interp_lorenzo<T, N>(conf, data, outSize);
        }
        return cmpData;
    }


    template<class T, uint N>
    void SZ_decompress_dispatcher(Config &conf, char *cmpData, size_t cmpSize, T *decData) {
        if (conf.absErrorBound == 0) {
            auto zstd = Lossless_zstd();
            auto zstdOut = zstd.decompress((uchar *) cmpData, cmpSize);
            memcpy(decData, zstdOut, conf.num * sizeof(T));
        } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
            SZ_decompress_LorenzoReg<T, N>(conf, cmpData, cmpSize, decData);
        } else if (conf.cmprAlgo == ALGO_INTERP) {
            SZ_decompress_Interp<T, N>(conf, cmpData, cmpSize, decData);
        } else {
            printf("SZ_decompress_dispatcher, Method not supported\n");
            exit(0);
        }

    }
}
#endif