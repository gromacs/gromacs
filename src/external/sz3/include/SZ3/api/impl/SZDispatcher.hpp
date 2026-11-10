#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

#include <cmath>
#include <memory>
#include <stdexcept>

#include "SZ3/api/impl/SZAlgoBioMD.hpp"
#include "SZ3/api/impl/SZAlgoInterp.hpp"
#include "SZ3/api/impl/SZAlgoLorenzoReg.hpp"
#include "SZ3/api/impl/SZAlgoNopred.hpp"
#include "SZ3/deprecated/HuffmanEncoderV1.hpp"
#include "SZ3/lossless/Lossless_zstd.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Statistic.hpp"

namespace SZ3 {
template <class T, uint N>
size_t SZ_compress_dispatcher(Config &conf, const T *data, uchar *cmpData, size_t cmpCap) {
    assert(N == conf.N);
    if (conf.cmprAlgo == ALGO_BIOMDXTC) {
        throw std::invalid_argument("SZ3: ALGO_BIOMDXTC was removed in 3.4.0; use ALGO_BIOMD");
    }
    calAbsErrorBound(conf, data);
    size_t cmpSize = 0;

    // a bound that is 0, negative or not finite (a relative one over data with Inf among them): lossless only mode
    if (!(conf.absErrorBound > 0) || !std::isfinite(conf.absErrorBound)) {
        conf.cmprAlgo = ALGO_LOSSLESS;
    }

    // do lossy compression
    bool isCmpCapSufficient = true;
    if (conf.cmprAlgo != ALGO_LOSSLESS) {
        try {
            // BIOMD and LORENZO_REG only read their input: LORENZO_REG's blockwise decomposition works on a copy padded
            // by at least 2, its fallback Lorenzo's padding
            if (conf.cmprAlgo == ALGO_BIOMD) {
                return SZ_compress_bioMD<T, N>(conf, data, cmpData, cmpCap);
            } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
                cmpSize = SZ_compress_LorenzoReg<T, N>(conf, const_cast<T *>(data), cmpData, cmpCap);
            } else {
                // the others write their reconstruction into the data they are given, so they get a copy
                std::vector<T> dataCopy(data, data + conf.num);
                if (conf.cmprAlgo == ALGO_INTERP) {
                    cmpSize = SZ_compress_Interp<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
                } else if (conf.cmprAlgo == ALGO_INTERP_LORENZO) {
                    cmpSize = SZ_compress_Interp_lorenzo<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
                } else if (conf.cmprAlgo == ALGO_NOPRED) {
                    cmpSize = SZ_compress_nopred<T, N>(conf, dataCopy.data(), cmpData, cmpCap);
                } else {
                    throw std::invalid_argument("Unknown compression algorithm");
                }
            }

        } catch (std::length_error &e) {
            if (std::string(e.what()) == SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH) {
                isCmpCapSufficient = false;
                // printf("SZ is downgraded to lossless mode because the buffer for compressed data is not large enough.\n");
            } else {
                throw;
            }
        }
    }

    // do lossless only compression if 1) cmpr algorithm is lossless or 2) compressed buffer not large enough for lossy
    if (conf.cmprAlgo == ALGO_LOSSLESS || !isCmpCapSufficient) {
        conf.cmprAlgo = ALGO_LOSSLESS;
        auto zstd = Lossless_zstd();
        return zstd.compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), cmpData, cmpCap);
    }

    // if lossy compression ratio < 3, test if lossless only mode has a better ratio than lossy
    if (conf.num * sizeof(T) / 1.0 / cmpSize < 3) {
        auto zstd = Lossless_zstd();
        // No bigger than the lossy stream: zstd gives up as soon as its output could not be smaller.
        std::unique_ptr<uchar[]> zstdCmpData(new uchar[cmpSize]);
        try {
            size_t zstdCmpSize =
                zstd.compress(reinterpret_cast<const uchar *>(data), conf.num * sizeof(T), zstdCmpData.get(), cmpSize);
            if (zstdCmpSize < cmpSize) {
                conf.cmprAlgo = ALGO_LOSSLESS;
                memcpy(cmpData, zstdCmpData.get(), zstdCmpSize);
                cmpSize = zstdCmpSize;
            }
        } catch (const std::length_error &) {
        }
    }
    return cmpSize;
}

template <class T, uint N>
void SZ_decompress_dispatcher(Config &conf, const uchar *cmpData, size_t cmpSize, T *decData) {
    {
        // backward compatibility for v3.3.2; its Huffman streams are HuffmanEncoderV1's.
        if (conf.sz3DataVer < versionInt("3.4.0")) {
            if (conf.cmprAlgo == ALGO_LORENZO_REG)
                return SZ_decompress_LorenzoReg<T, N, HuffmanEncoderV1<int>>(conf, cmpData, cmpSize, decData);
            if (conf.cmprAlgo == ALGO_INTERP)
                return SZ_decompress_Interp<T, N, HuffmanEncoderV1<int>>(conf, cmpData, cmpSize, decData);
            if (conf.cmprAlgo != ALGO_LOSSLESS)
                throw std::invalid_argument("SZ3: " + enum_to_string(static_cast<ALGO>(conf.cmprAlgo), ALGO_MAP) +
                                            " data of version " + versionStr(conf.sz3DataVer) +
                                            " is not supported; use SZ3 v" + versionStr(conf.sz3DataVer) +
                                            " to decompress it");
        }
    }
    if (conf.cmprAlgo == ALGO_LOSSLESS) {
        auto zstd = Lossless_zstd();
        auto decDataPos = reinterpret_cast<uchar *>(decData);
        size_t decDataSize = zstd.decompress(cmpData, cmpSize, decDataPos, conf.num * sizeof(T));
        if (decDataSize != conf.num * sizeof(T)) {
            throw std::runtime_error("Decompressed data size does not match the original data size");
        }
    } else if (conf.cmprAlgo == ALGO_LORENZO_REG) {
        SZ_decompress_LorenzoReg<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_INTERP) {
        SZ_decompress_Interp<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_NOPRED) {
        SZ_decompress_nopred<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == ALGO_BIOMD) {
        SZ_decompress_bioMD<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        throw std::invalid_argument("Unknown compression algorithm");
    }
}
}  // namespace SZ3
#endif
