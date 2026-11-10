#include "H5Z_SZ3.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>

// The message is a printf format, so anything that is not a literal goes through "%s".
#define H5Z_SZ_PUSH_AND_GOTO(MAJ, MIN, RET, ...)                                                  \
    do {                                                                                          \
        H5Epush(H5E_DEFAULT, __FILE__, _funcname_, __LINE__, H5E_ERR_CLS, MAJ, MIN, __VA_ARGS__); \
        return RET;                                                                               \
    } while (0)

// The filter's own callbacks. HDF5 reaches them through H5Z_SZ3 below, so they stay out of the
// public header, where every consumer of the filter would see them declared and never defined.
extern "C" {
static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id);

static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes,
                             size_t* buf_size, void** buf);
}

const H5Z_class2_t H5Z_SZ3[1] = {{
    H5Z_CLASS_T_VERS,                                       /* H5Z_class_t version */
    H5Z_FILTER_SZ3,                                         /* Filter id number */
    1,                                                      /* encoder_present flag (set to true) */
    1,                                                      /* decoder_present flag (set to true) */
    // Written into every dataset that uses the filter, and quoted back by HDF5 when the filter
    // is missing, so it has to say where to get it.
    "H5Z-SZ3-" SZ3_VER " (data format " SZ3_DATA_VER "); see "
    "https://github.com/szcompressor/SZ3/tree/master/tools/H5Z-SZ3",
    NULL,                                                   /* The "can apply" callback */
    H5Z_sz3_set_local,                                      /* The "set local" callback */
    static_cast<H5Z_func_t>(H5Z_filter_sz3),                /* The actual filter function */
}};

HDF5SZ3_EXPORT H5PL_type_t H5PLget_plugin_type(void) { return H5PL_TYPE_FILTER; }

HDF5SZ3_EXPORT const void* H5PLget_plugin_info(void) { return H5Z_SZ3; }

// set_SZ3_conf_to_H5() and get_sz3_conf_from_cdvalues() keep cd_values little-endian, even on big-endian hosts.
static void get_sz3_conf_from_cdvalues(const unsigned int* cd, size_t cd_nelmts, SZ3::Config& conf) {
    std::vector<unsigned char> buffer(cd_nelmts * sizeof(unsigned int));
    auto out = buffer.data();
    SZ3::write(cd, cd_nelmts, out);
    const unsigned char* bytes = buffer.data();
    size_t len = buffer.size();
    {
        // backward compatibility for v3.2.0 to v3.3.0; they store the magic number and the data version first.
        if (cd_nelmts > 1 && cd[0] == SZ3_MAGIC_NUMBER) {
            throw std::invalid_argument("SZ3 HDF5 filter: data is in SZ3 data format v" + versionStr(cd[1]) +
                                        ", this build reads v" SZ3_DATA_VER);
        }
    }
    {
        // backward compatibility for v3.3.2; it stores no version, so the first byte is the Config's length, never 0.
        if (bytes[0] != 0) {
            conf.load(bytes, len);
            return;
        }
    }
    // Another data version may lay out the Config differently.
    if (cd[0] < versionInt(SZ3_DATA_VER_OLDEST) || cd[0] > versionInt(SZ3_DATA_VER))
        throw std::invalid_argument("SZ3 HDF5 filter: data is in SZ3 data format v" + versionStr(cd[0]) +
                                    ", this build reads v" SZ3_DATA_VER_OLDEST " to v" SZ3_DATA_VER);
    bytes += sizeof(unsigned int);
    len -= sizeof(unsigned int);
    conf.load(bytes, len);
}

// Do not use H5Zfilter_avail() here: it answers for the library, not for this property list.
static bool sz3_filter_on_plist(hid_t propertyList) {
    const int nfilters = H5Pget_nfilters(propertyList);
    for (int i = 0; i < nfilters; i++) {
        unsigned int flags = 0;
        size_t cd_nelmts = 0;
        unsigned int cd_values[1] = {0};
        if (H5Z_FILTER_SZ3 ==
            H5Pget_filter2(propertyList, static_cast<unsigned>(i), &flags, &cd_nelmts, cd_values, 0, NULL, NULL)) {
            return true;
        }
    }
    return false;
}

herr_t set_SZ3_conf_to_H5(hid_t propertyList, SZ3::Config& conf) {
    static char const* _funcname_ = "set_SZ3_conf_to_H5";

    std::vector<unsigned int> cd_values(1 + conf.size_est(), 0);
    cd_values[0] = versionInt(SZ3_DATA_VER);
    auto pos = reinterpret_cast<unsigned char*>(cd_values.data() + 1);
    cd_values.resize(1 + (conf.save(pos) + sizeof(unsigned int) - 1) / sizeof(unsigned int));
#if SZ3_BIG_ENDIAN
    for (size_t i = 1; i < cd_values.size(); i++) cd_values[i] = SZ3::byteswap(cd_values[i]);
#endif
    size_t cd_nelmts = cd_values.size();

    if (sz3_filter_on_plist(propertyList)) {
        // Keep the flags the filter was set with: H5Z_FLAG_OPTIONAL lets HDF5 store a chunk raw when SZ3 fails.
        unsigned int flags = H5Z_FLAG_MANDATORY;
        size_t no_values = 0;
        if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, &flags, &no_values, NULL, 0, NULL, NULL) ||
            0 > H5Pmodify_filter(propertyList, H5Z_FILTER_SZ3, flags, cd_nelmts, cd_values.data())) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
        }
    } else {
        // filter not set, set filter. Notice that calling H5Pset_filter twice with the same filter id
        // will cause unexpected errors for decompression
        if (0 > H5Pset_filter(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data())) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
        }
    }

    return 1;
}

herr_t H5Pset_sz3(hid_t propertyList, int cmprAlgo, int errorBoundMode, double absErrorBound, double relErrorBound,
                  double psnrErrorBound, double l2normErrorBound) {
    static char const* _funcname_ = "H5Pset_sz3";
    static_assert(H5Z_SZ3_EB_ABS == SZ3::EB_ABS && H5Z_SZ3_EB_REL == SZ3::EB_REL && H5Z_SZ3_EB_PSNR == SZ3::EB_PSNR &&
                      H5Z_SZ3_EB_L2NORM == SZ3::EB_L2NORM && H5Z_SZ3_EB_ABS_AND_REL == SZ3::EB_ABS_AND_REL &&
                      H5Z_SZ3_EB_ABS_OR_REL == SZ3::EB_ABS_OR_REL,
                  "H5Z_SZ3.hpp error-bound modes must match SZ3::EB");
    static_assert(H5Z_SZ3_ALGO_LORENZO_REG == SZ3::ALGO_LORENZO_REG &&
                      H5Z_SZ3_ALGO_INTERP_LORENZO == SZ3::ALGO_INTERP_LORENZO &&
                      H5Z_SZ3_ALGO_INTERP == SZ3::ALGO_INTERP && H5Z_SZ3_ALGO_NOPRED == SZ3::ALGO_NOPRED &&
                      H5Z_SZ3_ALGO_LOSSLESS == SZ3::ALGO_LOSSLESS && H5Z_SZ3_ALGO_BIOMD == SZ3::ALGO_BIOMD &&
                      H5Z_SZ3_ALGO_BIOMDXTC == SZ3::ALGO_BIOMDXTC,
                  "H5Z_SZ3.hpp algorithms must match SZ3::ALGO");
    if (cmprAlgo < H5Z_SZ3_ALGO_LORENZO_REG || cmprAlgo > H5Z_SZ3_ALGO_BIOMDXTC) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "unknown SZ3 algorithm %d", cmprAlgo);
    }
    if (cmprAlgo == H5Z_SZ3_ALGO_BIOMDXTC) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1,
                             "SZ3 ALGO_BIOMDXTC was removed in 3.4.0; use H5Z_SZ3_ALGO_BIOMD");
    }
    if (errorBoundMode < H5Z_SZ3_EB_ABS || errorBoundMode > H5Z_SZ3_EB_ABS_OR_REL) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "unknown SZ3 error-bound mode %d", errorBoundMode);
    }
    for (double bound : {absErrorBound, relErrorBound, psnrErrorBound, l2normErrorBound}) {
        if (!std::isfinite(bound) || bound < 0) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "SZ3 error bounds must be finite and not negative");
        }
    }
    SZ3::Config conf;
    conf.cmprAlgo = static_cast<uint8_t>(cmprAlgo);
    conf.errorBoundMode = static_cast<uint8_t>(errorBoundMode);
    conf.absErrorBound = absErrorBound;
    conf.relErrorBound = relErrorBound;
    conf.psnrErrorBound = psnrErrorBound;
    conf.l2normErrorBound = l2normErrorBound;
    return set_SZ3_conf_to_H5(propertyList, conf) > 0 ? 1 : -1;
}

herr_t get_SZ3_conf_from_H5(hid_t propertyList, SZ3::Config& conf) {
    static char const* _funcname_ = "get_SZ3_conf_from_H5";

    // H5Pget_filter_by_id fails when the list carries no SZ3 filter, which is not an error to report.
    if (!sz3_filter_on_plist(propertyList)) {
        return 0;
    }
    // Ask for the count before the values: a stored config carries the dataset's dimensions, so its
    // length is not the caller's to guess, and HDF5 reports the real count even when it filled less.
    size_t cd_nelmts = 0;
    unsigned int probe[1] = {0};
    if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, &cd_nelmts, probe, 0, NULL, NULL)) {
        return -1;
    }
    if (cd_nelmts > 0) {
        std::vector<unsigned int> cd_values(cd_nelmts, 0);
        if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, &cd_nelmts, cd_values.data(), 0,
                                    NULL, NULL)) {
            return -1;
        }
        try {
            SZ3::Config loaded;
            get_sz3_conf_from_cdvalues(cd_values.data(), cd_values.size(), loaded);
            conf = loaded;
        } catch (const std::exception& e) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "%s", e.what());
        }
    }
    return 1;
}

static herr_t H5Z_sz3_set_local_impl(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    // printf("start H5Z_sz3_set_local\n");

    // printf("start in H5Z_sz3_set_local, dcpl_id = %d\n", dcpl_id);
    static char const* _funcname_ = "H5Z_sz3_set_local";

    // herr_t ret = H5Zregister(H5Z_SZ3);

    SZ3::Config conf;
    if (get_SZ3_conf_from_H5(dcpl_id, conf) < 0) return -1;

    int ndims;
    hsize_t dims_all[H5S_MAX_RANK];
    if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims_all, NULL)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
    std::vector<size_t> dims(dims_all, dims_all + ndims);
    // The filter uses a chunk's bytes as host values, so the datatype has to be identical to a host
    // type; H5Tequal also compares byte order, precision, offset and padding.
    const std::pair<hid_t, uint8_t> host_types[] = {{H5T_NATIVE_FLOAT, SZ_FLOAT}, {H5T_NATIVE_DOUBLE, SZ_DOUBLE},
                                                    {H5T_NATIVE_INT8, SZ_INT8},   {H5T_NATIVE_UINT8, SZ_UINT8},
                                                    {H5T_NATIVE_INT16, SZ_INT16}, {H5T_NATIVE_UINT16, SZ_UINT16},
                                                    {H5T_NATIVE_INT32, SZ_INT32}, {H5T_NATIVE_UINT32, SZ_UINT32},
                                                    {H5T_NATIVE_INT64, SZ_INT64}, {H5T_NATIVE_UINT64, SZ_UINT64}};
    const auto host_type = std::find_if(std::begin(host_types), std::end(host_types),
                                        [&](const auto& t) { return H5Tequal(type_id, t.first) > 0; });
    if (host_type == std::end(host_types))
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, -1,
                             "datatype must be the host's float, double, or 1, 2, 4 or 8-byte integer");
    conf.dataType = host_type->second;
    // update conf with dims
    conf.setDims(std::begin(dims), std::end(dims));
    if (conf.N > 4)
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "SZ3 compresses at most 4 dimensions longer than 1, not %d",
                             static_cast<int>(conf.N));
    if (conf.cmprAlgo == SZ3::ALGO_BIOMD && (conf.N > 3 || conf.dims[conf.N - 1] != 3)) {
        std::string chunk;
        for (size_t d : dims) chunk += (chunk.empty() ? "" : ", ") + std::to_string(d);
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1,
                             "SZ3 ALGO_BIOMD needs chunks whose last dimension is 3, (frames, atoms, 3); these are (%s)",
                             chunk.c_str());
    }
    //  need to update magic number and data version,
    //  as the config may be from cd_values passed by users
    conf.sz3MagicNumber = SZ3_MAGIC_NUMBER;
    conf.sz3DataVer = versionInt(SZ3_DATA_VER);

    if (set_SZ3_conf_to_H5(dcpl_id, conf) <= 0) return -1;
    return 1;
}

template <typename T>
void process_data(SZ3::Config& conf, void** buf, size_t* buf_size, size_t nbytes, bool is_decompress) {
    if (is_decompress) {
        // HDF5 frees what this returns, so it has to come from malloc. On null SZ_decompress would
        // allocate with new[] instead, and that pairing is undefined.
        std::unique_ptr<T, decltype(&free)> processedData(static_cast<T*>(malloc(conf.num * sizeof(T))), &free);
        if (!processedData) throw std::bad_alloc();
        T* decData = processedData.get();
        SZ_decompress(conf, static_cast<char*>(*buf), nbytes, decData);
        free(*buf);
        *buf = processedData.release();
        *buf_size = conf.num * sizeof(T);
    } else {
        if (nbytes != conf.num * sizeof(T))
            throw std::invalid_argument("SZ3 HDF5 filter: a filter before SZ3 changed the chunk's size");
        // The bound assumes the payload fits in the raw size, so leave headroom on top of it for
        // algorithms whose output can reach or exceed that.
        size_t cmpCap = std::max(SZ3::SZ_compress_size_bound<T>(conf), sizeof(T) * conf.num * 2);
        std::unique_ptr<char, decltype(&free)> cmpData(static_cast<char*>(malloc(cmpCap)), &free);
        if (!cmpData) throw std::bad_alloc();
        *buf_size = SZ_compress(conf, static_cast<T*>(*buf), cmpData.get(), cmpCap);
        free(*buf);
        *buf = cmpData.release();
    }
}

/**
 * https://docs.hdfgroup.org/hdf5/v1_14/_f_i_l_t_e_r.html
 * The flags, cd_nelmts, and cd_values are the same as for the H5Pset_filter() function with the additional flag
 * H5Z_FLAG_REVERSE which is set when the filter is called as part of the input pipeline. The input buffer is pointed to
 * by *buf and has a total size of *buf_size bytes but only nbytes are valid data. The filter should perform the
 * transformation in place if possible and return the number of valid bytes or zero for failure. If the transformation
 * cannot be done in place then the filter should allocate a new buffer with malloc() and assign it to *buf, assigning
 * the allocated size of that buffer to *buf_size. The old buffer should be freed by calling free().
 */
static size_t H5Z_filter_sz3_impl(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[],
                                  size_t nbytes, size_t* buf_size, void** buf) {
    // printf("start H5Z_filter_sz3\n");

    if (cd_nelmts == 0) // this is special data such as string, which should not be treated as values.
        return nbytes;

    bool is_decompress = flags & H5Z_FLAG_REVERSE;
    SZ3::Config conf;

    if (is_decompress) {
        // The chunk carries its own Config, which gives the element count and type to decompress.
        uint32_t magic = 0, dataVer = 0;
        uint64_t cmpDataSize = 0;
        auto pos = static_cast<const unsigned char*>(*buf);
        if (nbytes >= 16) {
            SZ3::read(magic, pos);
            SZ3::read(dataVer, pos);
            SZ3::read(cmpDataSize, pos);
        }
        if (magic != SZ3_MAGIC_NUMBER) {
            {
                // backward compatibility for v3.3.2; it stores chunks of fewer than 20 elements raw, with no header,
                // so in fewer than 20 * 8 bytes.
                if (nbytes < 20 * sizeof(double)) return nbytes;
            }
            throw std::invalid_argument("SZ3 HDF5 filter: chunk is not in a data format this SZ3 reads, or is damaged");
        }
        if (dataVer < versionInt(SZ3_DATA_VER_OLDEST) || dataVer > versionInt(SZ3_DATA_VER))
            throw std::invalid_argument("SZ3 HDF5 filter: data is in SZ3 data format v" + versionStr(dataVer) +
                                        ", this build reads v" SZ3_DATA_VER_OLDEST " to v" SZ3_DATA_VER +
                                        (dataVer > versionInt(SZ3_DATA_VER)
                                             ? "; upgrade SZ3 to read it"
                                             : "; use SZ3 v" + versionStr(dataVer) + " to read it"));
        if (cmpDataSize > nbytes - 16) throw std::invalid_argument("SZ3 HDF5 filter: chunk is truncated");
        size_t remaining = nbytes - 16 - cmpDataSize;
        pos += cmpDataSize;
        conf.load(pos, remaining);
    } else {
        get_sz3_conf_from_cdvalues(cd_values, cd_nelmts, conf);
    }

    switch (conf.dataType) {
        case SZ_FLOAT:
            process_data<float>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_DOUBLE:
            process_data<double>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT8:
            process_data<int8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT8:
            process_data<uint8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT16:
            process_data<int16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT16:
            process_data<uint16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT32:
            process_data<int32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT32:
            process_data<uint32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT64:
            process_data<int64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT64:
            process_data<uint64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        default:
            throw std::invalid_argument("SZ3 HDF5 filter: unknown datatype in cd_values");
    }
    return *buf_size;
}

// HDF5 calls these from C, where an exception that escapes ends the application. Zero is how a
// filter reports failure; set_local uses a negative return.
static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes,
                             size_t* buf_size, void** buf) {
    static char const* _funcname_ = "H5Z_filter_sz3";
    try {
        return H5Z_filter_sz3_impl(flags, cd_nelmts, cd_values, nbytes, buf_size, buf);
    } catch (const std::exception& e) {
#if H5_VERSION_GE(1, 14, 5)
        // 1.14.5 and 1.14.6 pause the error stack around filters, which drops what is pushed below. Print only
        // then: stdout would land in h5dump's output.
        bool paused = false;
        if (H5Eis_paused(H5E_DEFAULT, &paused) >= 0 && paused) fprintf(stderr, "%s\n", e.what());
#endif
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, 0, "%s", e.what());
    } catch (...) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, 0, "unknown error");
    }
}

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    static char const* _funcname_ = "H5Z_sz3_set_local";
    try {
        return H5Z_sz3_set_local_impl(dcpl_id, type_id, chunk_space_id);
    } catch (const std::exception& e) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, -1, "%s", e.what());
    } catch (...) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, -1, "unknown error");
    }
}
