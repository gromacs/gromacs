//
// Created by arham23 on 2/8/22.
//

#include <memory>
#include <fstream>
#include <iterator>
#include "H5PLextern.h"
#include "H5Z_SZ3.hpp"

hid_t H5Z_SZ_ERRCLASS = -1;

//h5repack -f UD=32024,0 /home/arham23/Software/SZ3/test/testfloat_8_8_128.dat.h5 tf_8_8_128.dat.sz.h5

//filter definition
const H5Z_class2_t H5Z_SZ3[1] = {{
                                     H5Z_CLASS_T_VERS,              /* H5Z_class_t version */
                                     (H5Z_filter_t) H5Z_FILTER_SZ3, /* Filter id number */
                                     1,              /* encoder_present flag (set to true) */
                                     1,              /* decoder_present flag (set to true) */
                                     "SZ3 compressor/decompressor for floating-point data.", /* Filter name for debugging */
                                     NULL,                          /* The "can apply" callback */
                                     H5Z_sz3_set_local,                          /* The "set local" callback */
                                     (H5Z_func_t) H5Z_filter_sz3,   /* The actual filter function */
                                 }
};

H5PL_type_t H5PLget_plugin_type(void) {
    return H5PL_TYPE_FILTER;
}

const void *H5PLget_plugin_info(void) {
    return H5Z_SZ3;
}

herr_t set_SZ3_conf_to_H5(const hid_t propertyList, SZ3::Config & conf) {
    static char const *_funcname_ = "set_SZ3_conf_to_H5";
    
    //save conf into cd_values
    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts, 0);
    auto buffer = (unsigned char *) (cd_values.data());
    
    conf.save(buffer);
    
    /* update cd_values for the filter */
    if (0 > H5Pmodify_filter(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data()))
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
    
    return (herr_t) 1;
}

herr_t get_SZ3_conf_from_H5(const hid_t propertyList, SZ3::Config & conf) {
    static char const *_funcname_ = "get_SZ3_conf_from_H5";
    
    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts, 0);
    
    //read cd_values from HDF5
    //note that cd_nelmts must be non-zero, otherwise, cd_values cannot be filled.
    if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, &cd_nelmts, cd_values.data(), 0, NULL, NULL))
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current SZ cd_values");
    
    //load cd_values into config
    if (cd_nelmts != 0) {
        auto buffer = (const unsigned char *) (cd_values.data());
        conf.load(buffer);
    }
    return (herr_t) 1;
}

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    
    printf("start H5Z_sz3_set_local\n");
    
    //printf("start in H5Z_sz3_set_local, dcpl_id = %d\n", dcpl_id);
    static char const *_funcname_ = "H5Z_sz3_set_local";
    
    // herr_t ret = H5Zregister(H5Z_SZ3);
    
    SZ3::Config conf;
    get_SZ3_conf_from_H5(dcpl_id, conf);
    
    //read datatype and dims from HDF5
    H5T_class_t dclass;
    if (0 > (dclass = H5Tget_class(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");
    
    size_t dsize;
    if (0 == (dsize = H5Tget_size(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");
    
    int ndims;
    hsize_t dims_all[H5S_MAX_RANK];
    if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims_all, 0)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
    std::vector<size_t> dims(dims_all, dims_all + ndims);
    
    //update conf with datatype
    conf.dataType = SZ_FLOAT;
    if (dclass == H5T_FLOAT)
        conf.dataType = dsize == 4 ? SZ_FLOAT : SZ_DOUBLE;
    else if (dclass == H5T_INTEGER) {
        H5T_sign_t dsign;
        if (0 > (dsign = H5Tget_sign(type_id)))
            H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");
        if (dsign == H5T_SGN_NONE) //unsigned
        {
            switch (dsize) {
                case 1:conf.dataType = SZ_UINT8;
                    break;
                case 2:conf.dataType = SZ_UINT16;
                    break;
                case 4:conf.dataType = SZ_UINT32;
                    break;
                case 8:conf.dataType = SZ_UINT64;
                    break;
            }
        } else {
            switch (dsize) {
                case 1:conf.dataType = SZ_INT8;
                    break;
                case 2:conf.dataType = SZ_INT16;
                    break;
                case 4:conf.dataType = SZ_INT32;
                    break;
                case 8:conf.dataType = SZ_INT64;
                    break;
            }
        }
    } else {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
    }
    
    //update conf with dims
    conf.setDims(std::begin(dims), std::end(dims));
    
    set_SZ3_conf_to_H5(dcpl_id, conf);
    
    return (herr_t) 1;
}

template<typename T>
void process_data(SZ3::Config &conf, void **buf, size_t *buf_size, size_t nbytes, bool is_decompress) {
    if (is_decompress) {
        T *processedData = (T *) malloc(conf.num * sizeof(T));
        SZ_decompress(conf, (char *) *buf, nbytes, processedData);
        free(*buf);
        *buf = processedData;
        *buf_size = conf.num * sizeof(T);
    } else {
        size_t outSize = 0;
        char *processedData = SZ_compress(conf, (T *) *buf, outSize);
        free(*buf);
        *buf = processedData;
        *buf_size = outSize;
    }
}

/**
 * https://docs.hdfgroup.org/hdf5/v1_14/_f_i_l_t_e_r.html
 * The flags, cd_nelmts, and cd_values are the same as for the H5Pset_filter() function with the additional flag H5Z_FLAG_REVERSE which is set when the filter is called as part of the input pipeline.
 * The input buffer is pointed to by *buf and has a total size of *buf_size bytes but only nbytes are valid data.
 * The filter should perform the transformation in place if possible and return the number of valid bytes or zero for failure.
 * If the transformation cannot be done in place then the filter should allocate a new buffer with malloc() and assign it to *buf, assigning the allocated size of that buffer to *buf_size.
 * The old buffer should be freed by calling free().
 */
static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf) {
    printf("start H5Z_filter_sz3\n");
    
    if (cd_nelmts == 0) //this is special data such as string, which should not be treated as values.
        return nbytes;
    
    SZ3::Config conf;
    
    auto buffer = (const unsigned char *) (cd_values);
    conf.load(buffer);
//    conf.print();
    
    if (conf.num < 20)
        return nbytes;
    
    bool is_decompress = flags & H5Z_FLAG_REVERSE;
    switch (conf.dataType) {
        case SZ_FLOAT: process_data<float>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_DOUBLE: process_data<double>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT8: process_data<int8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT8: process_data<uint8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT16: process_data<int16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT16: process_data<uint16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT32: process_data<int32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT32: process_data<uint32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT64: process_data<int64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT64: process_data<uint64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        default: std::cerr << (is_decompress ? "Decompression" : "Compression") << " Error: Unknown Datatype" << std::endl;
            std::exit(EXIT_FAILURE);
    }
    return *buf_size;
}