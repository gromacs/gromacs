/**
 * @file sz.hpp
 * @brief SZ3 compression and decompression API.
 *
 * This header provides the main API functions for compressing and decompressing data using SZ3.
 *
 * Compressed Data Format of SZ3:
 * The compressed data is always stored in little-endian order.
 * The compressed data consists of three main sections:
 * 1. Header (16 bytes): Contains metadata about the compressed data.
 *    - Magic Number (4 bytes): Identifies the data as SZ3-compressed.
 *    - Version (4 bytes): Indicates the version of the SZ3 format.
 *    - Compressed Size (8 bytes): Specifies the size of the compressed payload.
 * 2. Compressed Payload: The actual compressed data.
 * 3. Configuration: Stores the compression configuration used.
 *
 * The layout can be visualized as follows:
 * [ Header (16 bytes) | Compressed Payload | Configuration ]
 */

#ifndef SZ3_SZ_HPP
#define SZ3_SZ_HPP

#include <memory>
#include <stdexcept>

#include "SZ3/api/impl/SZImpl.hpp"
#include "SZ3/version.hpp"

/**
 * Compresses the input data using the provided configuration and stores the result in a pre-allocated buffer.
 * @tparam T The data type of the source data.
 * @param config The compression configuration.
 * @param data Pointer to the source data array.
 * @param cmpData Pointer to the pre-allocated buffer for compressed data.
 * @param cmpCap The size of the pre-allocated buffer in bytes.
 * @return The size of the compressed data in bytes.
 * @example
 * SZ3::Config conf(100, 200, 300); // 300 is the fastest dimension
 * conf.errorBoundMode = SZ3::EB_ABS; // Refer to def.hpp for supported error bound modes
 * conf.absErrorBound = 1E-3; // Absolute error bound of 1e-3
 * size_t outSize = SZ_compress(conf, data, outBuff, outBuffCap);
 */
template <class T>
size_t SZ_compress(const SZ3::Config& config, const T* data, char* cmpData, size_t cmpCap) {
    using namespace SZ3;
    Config conf(config);

    if (cmpCap < SZ_compress_size_bound<T>(conf)) {
        throw std::invalid_argument(SZ3_ERROR_COMP_BUFFER_NOT_LARGE_ENOUGH);
    }

    auto cmpDataPos = reinterpret_cast<uchar*>(cmpData);

    // save 16 bytes header
    write(static_cast<uint32_t>(SZ3_MAGIC_NUMBER), cmpDataPos);  // magic number (4 bytes)
    write(versionInt(SZ3_DATA_VER), cmpDataPos);                 // data version (4 bytes)
    auto sizeHeaderPos = cmpDataPos;
    cmpDataPos += 8; // reserve space for cmp data size (8 bytes)

    // begin compression
    auto cmpDataCap = cmpCap - 16 - conf.size_est() * 2;
    uint64_t cmpDataSize = 0;
    if (conf.N == 1) {
        cmpDataSize = SZ_compress_impl<T, 1>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 2) {
        cmpDataSize = SZ_compress_impl<T, 2>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 3) {
        cmpDataSize = SZ_compress_impl<T, 3>(conf, data, cmpDataPos, cmpDataCap);
    } else if (conf.N == 4) {
        cmpDataSize = SZ_compress_impl<T, 4>(conf, data, cmpDataPos, cmpDataCap);
    } else {
        throw std::invalid_argument("Data dimension higher than 4 is not supported.");
    }

    // save compressed data size back in header
    write(cmpDataSize, sizeHeaderPos);

    // save config
    cmpDataPos += cmpDataSize;
    auto confSize = conf.save(cmpDataPos);

    return 16 + cmpDataSize + confSize;
}

/**
 * Compresses the input data using the provided configuration and returns a newly allocated buffer containing the compressed data.
 * @tparam T The data type of the source data.
 * @param config The compression configuration.
 * @param data Pointer to the source data array.
 * @param cmpSize Output parameter set to the size of the compressed data in bytes.
 * @return Pointer to the newly allocated buffer containing the compressed data. The caller is responsible for deleting this buffer using 'delete[]'.
 * @note This function allocates memory for the compressed data. Ensure to free it when no longer needed.
 */
template <class T>
char* SZ_compress(const SZ3::Config& config, const T* data, size_t& cmpSize) {
    using namespace SZ3;

    size_t bufferLen = SZ_compress_size_bound<T>(config);
    std::unique_ptr<char[]> buffer(new char[bufferLen]);
    cmpSize = SZ_compress(config, data, buffer.get(), bufferLen);

    return buffer.release();
}

/**
 * Decompresses the compressed data into a pre-allocated buffer using the configuration loaded from the compressed data.
 * @tparam T The data type of the decompressed data.
 * @param config Configuration placeholder that will be overwritten with the compression configuration from the compressed data.
 * @param cmpData Pointer to the compressed data.
 * @param cmpSize The size of the compressed data in bytes.
 * @param decData Reference to a pointer for the pre-allocated buffer for decompressed data. If null, a new buffer is allocated.
 * @example
 * auto decData = new float[100 * 200 * 300];
 * SZ3::Config conf;
 * SZ_decompress(conf, cmpData, cmpSize, decData);
 */
template <class T>
void SZ_decompress(SZ3::Config& config, const char* cmpData, size_t cmpSize, T*& decData) {
    using namespace SZ3;

    auto cmpDataPos = reinterpret_cast<const uchar*>(cmpData);

    // Header layout: magic number (4) + data version (4) + compressed payload size (8) = 16 bytes.
    if (cmpSize < 16) {
        throw std::out_of_range("SZ3: compressed data is smaller than the header");
    }

    read(config.sz3MagicNumber, cmpDataPos);
    if (config.sz3MagicNumber != SZ3_MAGIC_NUMBER) {
        throw std::invalid_argument("magic number mismatch, the input data is not compressed by SZ3");
    }

    read(config.sz3DataVer, cmpDataPos);
    if (config.sz3DataVer < versionInt(SZ3_DATA_VER_OLDEST) || config.sz3DataVer > versionInt(SZ3_DATA_VER)) {
        // Do not print here. In the HDF5 filter stdout is the user's output file.
        std::stringstream ss;
        ss << "SZ3 " << SZ3_VER << " reads data versions " << SZ3_DATA_VER_OLDEST << " to " << SZ3_DATA_VER
           << ", but this data is version " << versionStr(config.sz3DataVer) << ". ";
        if (config.sz3DataVer > versionInt(SZ3_DATA_VER))
            ss << "Upgrade SZ3 to decompress it.";
        else
            ss << "Use SZ3 v" << versionStr(config.sz3DataVer) << " to decompress it.";
        throw std::invalid_argument(ss.str());
    }

    uint64_t cmpDataSize = 0;
    read(cmpDataSize, cmpDataPos);

    // The compressed payload is followed by the serialized config; both must fit in the remaining bytes.
    if (cmpDataSize > cmpSize - 16) {
        throw std::out_of_range("SZ3: compressed payload size exceeds the buffer");
    }
    auto cmpConfPos = cmpDataPos + cmpDataSize;
    size_t confRemaining = cmpSize - 16 - cmpDataSize;
    config.load(cmpConfPos, confRemaining);

    if (decData == nullptr) {
        decData = new T[config.num];
    }
    if (config.N == 1) {
        SZ_decompress_impl<T, 1>(config, cmpDataPos, cmpDataSize, decData);
    } else if (config.N == 2) {
        SZ_decompress_impl<T, 2>(config, cmpDataPos, cmpDataSize, decData);
    } else if (config.N == 3) {
        SZ_decompress_impl<T, 3>(config, cmpDataPos, cmpDataSize, decData);
    } else if (config.N == 4) {
        SZ_decompress_impl<T, 4>(config, cmpDataPos, cmpDataSize, decData);
    } else {
        throw std::invalid_argument("Data dimension higher than 4 is not supported.");
    }
}

/**
 * Decompresses the compressed data into a buffer it allocates itself, using the configuration loaded from the
 * compressed data. Pass a buffer to the four-argument overload instead to decompress into memory you already hold.
 * @tparam T The data type of the decompressed data.
 * @param config Configuration placeholder that will be overwritten with the compression configuration from the compressed data.
 * @param cmpData Pointer to the compressed data.
 * @param cmpSize The size of the compressed data in bytes.
 * @return A new[] buffer of config.num elements, which the caller owns and must delete[]. There is no null return to
 * test: the four-argument overload throws, so a value arriving here is always a buffer.
 * @example
 * SZ3::Config conf;
 * float* decData = SZ_decompress<float>(conf, cmpData, cmpSize);
 * // conf now describes the data; conf.num is the element count of decData
 * delete[] decData;
 */
template <class T>
T* SZ_decompress(SZ3::Config& config, const char* cmpData, size_t cmpSize) {
    using namespace SZ3;
    T* decData = nullptr;
    try {
        SZ_decompress<T>(config, cmpData, cmpSize, decData);
    } catch (...) {
        // The four-argument overload allocates through the reference and can throw afterwards, and
        // the caller never sees the pointer when it does.
        delete[] decData;
        throw;
    }
    return decData;
}

#endif
