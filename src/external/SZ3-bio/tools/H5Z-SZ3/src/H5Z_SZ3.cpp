//
// Created by arham23 on 2/8/22.
//

#include <memory>
#include "H5Z_SZ3.hpp"
#include <fstream>
#include "H5PLextern.h"
#include "SZ3/api/sz.hpp"
#include "SZ3/utils/ByteUtil.hpp"



int sysEndianType = LITTLE_ENDIAN_SYSTEM;
int dataEndianType = LITTLE_ENDIAN_DATA;
hid_t H5Z_SZ_ERRCLASS = -1;

using namespace SZ3;

//h5repack -f UD=32024,0 /home/arham23/Software/SZ3/test/testfloat_8_8_128.dat.h5 tf_8_8_128.dat.sz.h5

//load from "sz3.config" in local directory if 1 else use default values or cd values

#define SZ3_CONFIG_PATH "SZ3_CONFIG_PATH"
SZ3::Config sz3_conf;
bool sz3_conf_loaded = false;
int MAX_CHUNK_SIZE = INT_MAX;

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

/*FILTER FUNCTIONS*/

/**
 * to be used in compression, and to be called outside H5Z_filter_sz().
 * */

void SZ_refreshDimForCdArray(int dataType, size_t old_cd_nelmts, unsigned int *old_cd_values, size_t *new_cd_nelmts, unsigned int **new_cd_values,
                             size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    unsigned char bytes[8] = {0};
    *new_cd_values = (unsigned int *) malloc(sizeof(unsigned int) * 16);
    memset(*new_cd_values, 0, sizeof(unsigned int) * 16);

    //correct dimension if needed
    size_t _r[5];
    filterDimension(r5, r4, r3, r2, r1, _r);
    size_t _r5 = _r[4];
    size_t _r4 = _r[3];
    size_t _r3 = _r[2];
    size_t _r2 = _r[1];
    size_t _r1 = _r[0];

    int i = 0;
    int oldDim = computeDimension(r5, r4, r3, r2, r1);
    int newDim = computeDimension(_r5, _r4, _r3, _r2, _r1);
    (*new_cd_values)[0] = newDim;
    (*new_cd_values)[1] = dataType;


    switch (newDim) {
        case 1:
            longToBytes_bigEndian(bytes, (uint64_t) r1);
            (*new_cd_values)[2] = bytesToInt_bigEndian(bytes);
            (*new_cd_values)[3] = bytesToInt_bigEndian(&bytes[4]);
            if (old_cd_nelmts == 0)
                *new_cd_nelmts = 4;
            else {
                (*new_cd_values)[4] = old_cd_values[0];
                (*new_cd_values)[5] = old_cd_values[1];
                (*new_cd_values)[6] = old_cd_values[2];
                (*new_cd_values)[7] = old_cd_values[3];
                (*new_cd_values)[8] = old_cd_values[4];
                (*new_cd_values)[9] = old_cd_values[5];
                (*new_cd_values)[10] = old_cd_values[6];
                (*new_cd_values)[11] = old_cd_values[7];
                (*new_cd_values)[12] = old_cd_values[8];
                *new_cd_nelmts = 13;
            }
            break;
        case 2:
            (*new_cd_values)[2] = (unsigned int) _r2;
            (*new_cd_values)[3] = (unsigned int) _r1;
            if (old_cd_nelmts == 0)
                *new_cd_nelmts = 4;
            else {
                (*new_cd_values)[4] = old_cd_values[0];
                (*new_cd_values)[5] = old_cd_values[1];
                (*new_cd_values)[6] = old_cd_values[2];
                (*new_cd_values)[7] = old_cd_values[3];
                (*new_cd_values)[8] = old_cd_values[4];
                (*new_cd_values)[9] = old_cd_values[5];
                (*new_cd_values)[10] = old_cd_values[6];
                (*new_cd_values)[11] = old_cd_values[7];
                (*new_cd_values)[12] = old_cd_values[8];
                *new_cd_nelmts = 13;
            }
            break;
        case 3:
            (*new_cd_values)[2] = (unsigned int) _r3;
            (*new_cd_values)[3] = (unsigned int) _r2;
            (*new_cd_values)[4] = (unsigned int) _r1;
            if (old_cd_nelmts == 0)
                *new_cd_nelmts = 5;
            else {
                (*new_cd_values)[5] = old_cd_values[0];
                (*new_cd_values)[6] = old_cd_values[1];
                (*new_cd_values)[7] = old_cd_values[2];
                (*new_cd_values)[8] = old_cd_values[3];
                (*new_cd_values)[9] = old_cd_values[4];
                (*new_cd_values)[10] = old_cd_values[5];
                (*new_cd_values)[11] = old_cd_values[6];
                (*new_cd_values)[12] = old_cd_values[7];
                (*new_cd_values)[13] = old_cd_values[8];
                *new_cd_nelmts = 14;
            }
            break;
        case 4:
            (*new_cd_values)[2] = (unsigned int) _r4;
            (*new_cd_values)[3] = (unsigned int) _r3;
            (*new_cd_values)[4] = (unsigned int) _r2;
            (*new_cd_values)[5] = (unsigned int) _r1;
            if (old_cd_nelmts == 0)
                *new_cd_nelmts = 6;
            else {
                (*new_cd_values)[6] = old_cd_values[0];
                (*new_cd_values)[7] = old_cd_values[1];
                (*new_cd_values)[8] = old_cd_values[2];
                (*new_cd_values)[9] = old_cd_values[3];
                (*new_cd_values)[10] = old_cd_values[4];
                (*new_cd_values)[11] = old_cd_values[5];
                (*new_cd_values)[12] = old_cd_values[6];
                (*new_cd_values)[13] = old_cd_values[7];
                (*new_cd_values)[14] = old_cd_values[8];
                *new_cd_nelmts = 15;
                break;
            }
        default:
            (*new_cd_values)[2] = (unsigned int) _r5;
            (*new_cd_values)[3] = (unsigned int) _r4;
            (*new_cd_values)[4] = (unsigned int) _r3;
            (*new_cd_values)[5] = (unsigned int) _r2;
            (*new_cd_values)[6] = (unsigned int) _r1;
            if (old_cd_nelmts == 0)
                *new_cd_nelmts = 7;
            else {
                (*new_cd_values)[7] = old_cd_values[0];
                (*new_cd_values)[8] = old_cd_values[1];
                (*new_cd_values)[9] = old_cd_values[2];
                (*new_cd_values)[10] = old_cd_values[3];
                (*new_cd_values)[11] = old_cd_values[4];
                (*new_cd_values)[12] = old_cd_values[5];
                (*new_cd_values)[13] = old_cd_values[6];
                (*new_cd_values)[14] = old_cd_values[7];
                (*new_cd_values)[15] = old_cd_values[8];
                *new_cd_nelmts = 16;
            }
    }
}


void
SZ_errConfigToCdArray(size_t *cd_nelmts, unsigned int **cd_values, int error_bound_mode, double abs_error, double rel_error, double l2normErrorBound,
                      double psnr) {
    *cd_values = (unsigned int *) malloc(sizeof(unsigned int) * 9);
    int k = 0;
    (*cd_values)[k++] = error_bound_mode;
    unsigned char b[8];
    doubleToBytes(b, abs_error);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b + 4);
    doubleToBytes(b, rel_error);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b + 4);
    doubleToBytes(b, l2normErrorBound);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b + 4);
    doubleToBytes(b, psnr);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b);
    (*cd_values)[k++] = bytesToInt32_bigEndian(b + 4);
    *cd_nelmts = k;
}

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {

    //printf("get into H5Z_sz3_set_local\n");
    detectSysEndianType();

    //printf("start in H5Z_sz3_set_local, dcpl_id = %d\n", dcpl_id);
    static char const *_funcname_ = "H5Z_sz3_set_local";
    size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0, dsize;

    int i, ndims, ndims_used = 0;
    hsize_t dims[H5S_MAX_RANK], dims_used[5] = {0, 0, 0, 0, 0};
    herr_t retval = 0;
    H5T_class_t dclass;
    H5T_sign_t dsign;
    unsigned int flags = 0;
    size_t mem_cd_nelmts = 9, cd_nelmts = 0;
    unsigned int mem_cd_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    //H5Z_FILTER_SZ
    //note that mem_cd_nelmts must be non-zero, otherwise, mem_cd_values cannot be filled.
    if (0 > H5Pget_filter_by_id(dcpl_id, H5Z_FILTER_SZ3, &flags, &mem_cd_nelmts, mem_cd_values, 0, NULL, NULL))
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current SZ cd_values");

    //set default value for error bound
    sz3_conf.errorBoundMode = EB_ABS;
    sz3_conf.absErrorBound = 1e-3;
    if (!sz3_conf_loaded) {
        if (const char *conf_file = std::getenv(SZ3_CONFIG_PATH)) {
            sz3_conf.loadcfg(conf_file);
            sz3_conf_loaded = true;
        }
    }
//    if (mem_cd_nelmts == 0) //this means that the error information is missing from the cd_values
//    {
//        //printf("mem_cd_nelmets is 0, so let's try using sz3.config to load error configuration....\n");
//        std::ifstream f(CONFIG_PATH);
//        if (f.good()) {
//            printf("sz3.config found!\n");
//            sz3_conf_loaded = 1;
//        } else
//            printf("sz3.config not found, using default parameters\n");
//        f.close();
//    } else //this means that the error information is included in the cd_values
//    {
//        sz3_conf_loaded = 0;
//        //printf("mem_cd_nelmets is non-zero, so let's use the parameters set through cd_values.....\n");
//    }
    herr_t ret = H5Zregister(H5Z_SZ3);

    int dataType = SZ_FLOAT;

    //printf("DC\n");
    if (0 > (dclass = H5Tget_class(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

    //printf("DS\n");
    if (0 == (dsize = H5Tget_size(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

    //printf("ND\n");
    if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims, 0)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");

    for (i = 0; i < ndims; i++)
        dims_used[i] = dims[i];


    //printf("NDIM: %i\n", ndims);
    //printf("N_USE: %i\n", ndims_used);
    //printf("DCLASS: %i\n", dclass);
    //printf("DSIZE: %zu\n", dsize);

    //for(i = 0; i < ndims_used; i++){
    //	printf("DIMS[%i] : %zu\n", i, dims_used[i]);
    //}
    //printf("\nDCEQ\n");

    if (dclass == H5T_FLOAT)
        dataType = dsize == 4 ? SZ_FLOAT : SZ_DOUBLE;
    else if (dclass == H5T_INTEGER) {
        if (0 > (dsign = H5Tget_sign(type_id)))
            H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");
        if (dsign == H5T_SGN_NONE) //unsigned
        {
            switch (dsize) {
                case 1:
                    dataType = SZ_UINT8;
                    break;
                case 2:
                    dataType = SZ_UINT16;
                    break;
                case 4:
                    dataType = SZ_UINT32;
                    break;
                case 8:
                    dataType = SZ_UINT64;
                    break;
            }
        } else {
            switch (dsize) {
                case 1:
                    dataType = SZ_INT8;
                    break;
                case 2:
                    dataType = SZ_INT16;
                    break;
                case 4:
                    dataType = SZ_INT32;
                    break;
                case 8:
                    dataType = SZ_INT64;
                    break;
            }
        }
    } else {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
    }

    unsigned int *cd_values = NULL;
    if (mem_cd_nelmts != 0 && mem_cd_nelmts != 9) {
        H5Epush(H5E_DEFAULT, __FILE__, "H5Z_sz3_set_local", __LINE__, H5E_ERR_CLS, H5E_ARGS, H5E_BADVALUE,
                "Wrong number of cd_values: The new version has 9 integer elements in cd_values. Please check 'test/print_h5repack_args' to get the correct cd_values.");
        H5Eprint(H5E_DEFAULT, stderr);
        return -1;
    }
    SZ_refreshDimForCdArray(dataType, mem_cd_nelmts, mem_cd_values, &cd_nelmts, &cd_values, dims_used[4], dims_used[3], dims_used[2], dims_used[1],
                            dims_used[0]);

    /* Now, update cd_values for the filter */
    if (0 > H5Pmodify_filter(dcpl_id, H5Z_FILTER_SZ3, flags, cd_nelmts, cd_values))
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");

    free(cd_values);

    retval = 1;
    done:
    return retval;
}


static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf) {
    //printf("get into H5Z_filter_sz3\n");
    size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
    int dimSize = 0, dataType = 0;

    if (cd_nelmts == 0) //this is special data such as string, which should not be treated as values.
        return nbytes;

    int withErrInfo = checkCDValuesWithErrors(cd_nelmts, cd_values);
    int error_mode = 0;
//    int cmp_algo = 1;
//    int interp_algo = 1;
    double abs_error = 0, rel_error = 0, l2norm_error = 0, psnr = 0;
    if (withErrInfo)
        SZ_cdArrayToMetaDataErr(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1, &error_mode, &abs_error, &rel_error,
                                &l2norm_error, &psnr);
    else
        SZ_cdArrayToMetaData(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1);

    /*int i=0;
    for(i=0;i<cd_nelmts;i++)
    	printf("cd_values[%d]=%u\n", i, cd_values[i]);
    printf("dimSize=%d, r1=%u, r2=%u, r3=%u, r4=%u, r5=%u\n", dimSize, r1, r2, r3, r4, r5);*/
    auto buff = (float *) *buf;
//    printf("dimSize=%d, r1=%u, r2=%u, r3=%u, r4=%u, r5=%u nbytes %lu %.3f %.3f %.3f \n", dimSize, r1, r2, r3, r4, r5, nbytes, buff[10000],
//           buff[20000], buff[30000]);

    size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);
    if (nbEle < 20)
        return nbytes;

    if (flags & H5Z_FLAG_REVERSE) {

        /* decompress data */
        SZ3::Config conf;

        switch (dataType) {
            case SZ_FLOAT: //FLOAT
            {
                float *f_decompressedData = new float[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, f_decompressedData);
                free(*buf);
                *buf = f_decompressedData;
                *buf_size = nbEle * sizeof(float);
                break;
            }

            case SZ_DOUBLE: //DOUBLE
            {
                double *d_decompressedData = new double[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, d_decompressedData);
                free(*buf);
                *buf = d_decompressedData;
                *buf_size = nbEle * sizeof(double);
                break;
            }

            case SZ_INT8: //INT 8
            {
                int8_t *c_decompressedData = new int8_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, c_decompressedData);
                free(*buf);
                *buf = c_decompressedData;
                *buf_size = nbEle * sizeof(int8_t);
                break;
            }

            case SZ_UINT8: //UINT 8
            {
                uint8_t *uc_decompressedData = new uint8_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, uc_decompressedData);
                free(*buf);
                *buf = uc_decompressedData;
                *buf_size = nbEle * sizeof(uint8_t);
                break;
            }

            case SZ_INT16: //INT 16
            {
                int16_t *s_decompressedData = new int16_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, s_decompressedData);
                free(*buf);
                *buf = s_decompressedData;
                *buf_size = nbEle * sizeof(int16_t);
                break;
            }

            case SZ_UINT16: //UINT 16
            {
                uint16_t *us_decompressedData = new uint16_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, us_decompressedData);
                free(*buf);
                *buf = us_decompressedData;
                *buf_size = nbEle * sizeof(uint16_t);
                break;
            }

            case SZ_INT32: //INT 32
            {
                int32_t *i_decompressedData = new int32_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, i_decompressedData);
                free(*buf);
                *buf = i_decompressedData;
                *buf_size = nbEle * sizeof(int32_t);
                break;
            }

            case SZ_UINT32: //UINT 32
            {
                uint32_t *ui_decompressedData = new uint32_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, ui_decompressedData);
                free(*buf);
                *buf = ui_decompressedData;
                *buf_size = nbEle * sizeof(uint32_t);
                break;
            }

            case SZ_INT64: //INT 64
            {
                int64_t *l_decompressedData = new int64_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, l_decompressedData);
                free(*buf);
                *buf = l_decompressedData;
                *buf_size = nbEle * sizeof(int64_t);
                break;
            }

            case SZ_UINT64: //UINT 64
            {
                uint64_t *ul_decompressedData = new uint64_t[nbEle];
                SZ_decompress(conf, (char *) *buf, nbytes, ul_decompressedData);
                free(*buf);
                *buf = ul_decompressedData;
                *buf_size = nbEle * sizeof(uint64_t);
                break;
            }

            default: {
                printf("Decompression Error: Unknown Datatype");
                exit(0);
            }
        }


    } else {
        /*compress data*/
        //based on # dimensions, get relevant dimensions and load config object with them
        if (dimSize <= 0) {
            printf("Error: Number of Dimensions is <= 0");
            exit(0);
        }
        //printf("\nDIMS_CMP:\n");
        //printf("r1 %u r2 %u r3 %u r4 %u r5 %u\n", r1,r2,r3,r4,r5);

        SZ3::Config conf(sz3_conf);
        std::vector<size_t> dims;
        if (r2 == 0) {
            dims = {r1};
        } else if (r3 == 0) {
            dims = {r2, r1};
        } else if (r4 == 0) {
            dims = {r3, r2, r1};
        } else if (r5 == 0) {
            dims = {r4, r3, r2, r1};
        } else {
            dims = {r5, r4, r3, r2, r1};
        }
        conf.setDims(dims.begin(), dims.end());



        //if config file found and no user defined params, read the config file
        if (withErrInfo) {
            if (error_mode < 0 || error_mode > 5) {
                printf("Invalid error mode: %i, error mode should be in [0,5]", error_mode);
                exit(0);
            }
            conf.errorBoundMode = error_mode;
            conf.absErrorBound = abs_error;
            conf.relErrorBound = rel_error;
            conf.l2normErrorBound = l2norm_error;
            conf.psnrErrorBound = psnr;
        }
        //printf("PARAMS: mode|%i, abs_eb|%f, rel_eb|%f, l2_eb|%f, psnr_eb|%f\n", error_mode, abs_error, rel_error, l2norm_error, psnr);

        size_t outSize = 0;
        char *compressedData = NULL;


        switch (dataType) {
            case SZ_FLOAT: //FLOAT
            {
                compressedData = SZ_compress(conf, (float *) *buf, outSize);
                break;
            }

            case SZ_DOUBLE: //DOUBLE
            {
                compressedData = SZ_compress(conf, (double *) *buf, outSize);
                break;
            }

            case SZ_INT8: //INT 8
            {
                compressedData = SZ_compress(conf, (int8_t *) *buf, outSize);
                break;
            }

            case SZ_UINT8: //UINT 8
            {
                compressedData = SZ_compress(conf, (uint8_t *) *buf, outSize);
                break;
            }

            case SZ_INT16: //INT 16
            {
                compressedData = SZ_compress(conf, (int16_t *) *buf, outSize);
                break;
            }

            case SZ_UINT16: //UINT 16
            {
                compressedData = SZ_compress(conf, (uint16_t *) *buf, outSize);
                break;
            }

            case SZ_INT32: //INT 32
            {
                compressedData = SZ_compress(conf, (int32_t *) *buf, outSize);
                break;
            }

            case SZ_UINT32: //UINT 32
            {
                compressedData = SZ_compress(conf, (uint32_t *) *buf, outSize);
                break;
            }

            case SZ_INT64: //INT 64
            {
                compressedData = SZ_compress(conf, (int64_t *) *buf, outSize);
                break;
            }

            case SZ_UINT64: //UINT 64
            {
                compressedData = SZ_compress(conf, (uint64_t *) *buf, outSize);
                break;
            }

            default: {
                printf("Compression Error: Unknown Datatype");
                exit(0);
            }
        }

        //printf("\nOS: %u \n", outSize);
        free(*buf);
        *buf = compressedData;
        *buf_size = outSize;


    }

    return *buf_size;
}

/*HELPER FUNCTIONS*/
//use to convert HDF5 cd_array to SZ params inside filter
void
SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int *dimSize, int *dataType, size_t *r5, size_t *r4, size_t *r3, size_t *r2,
                     size_t *r1) {
    assert(cd_nelmts >= 4);
    unsigned char bytes[8];
    *dimSize = cd_values[0];
    *dataType = cd_values[1];

    switch (*dimSize) {
        case 1:
            SZ3::int32ToBytes_bigEndian(bytes, cd_values[2]);
            SZ3::int32ToBytes_bigEndian(&bytes[4], cd_values[3]);
            if (sizeof(size_t) == 4)
                *r1 = (unsigned int) SZ3::bytesToInt64_bigEndian(bytes);
            else
                *r1 = (uint64_t) SZ3::bytesToInt64_bigEndian(bytes);
            *r2 = *r3 = *r4 = *r5 = 0;
            break;
        case 2:
            *r3 = *r4 = *r5 = 0;
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        case 3:
            *r4 = *r5 = 0;
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        case 4:
            *r5 = 0;
            *r4 = cd_values[5];
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        default:
            *r5 = cd_values[6];
            *r4 = cd_values[5];
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
    }
}

void
SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int *dimSize, int *dataType, size_t *r5, size_t *r4, size_t *r3, size_t *r2,
                        size_t *r1, int *error_bound_mode, double *abs_error, double *rel_error, double *l2norm_error, double *psnr) {
    //get dimension, datatype metadata from cd_values
    SZ_cdArrayToMetaData(cd_nelmts, cd_values, dimSize, dataType, r5, r4, r3, r2, r1);
    //read in error bound value information
    int dim = *dimSize;
    int k = dim == 1 ? 4 : dim + 2;
    unsigned char b[8];
    int32ToBytes_bigEndian(b, cd_values[k++]);
    *error_bound_mode = bytesToInt32_bigEndian(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    int32ToBytes_bigEndian(b + 4, cd_values[k++]);
    *abs_error = bytesToDouble(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    int32ToBytes_bigEndian(b + 4, cd_values[k++]);
    *rel_error = bytesToDouble(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    int32ToBytes_bigEndian(b + 4, cd_values[k++]);
    *l2norm_error = bytesToDouble(b);
    int32ToBytes_bigEndian(b, cd_values[k++]);
    int32ToBytes_bigEndian(b + 4, cd_values[k++]);
    *psnr = bytesToDouble(b);
}

void SZ_copymetaDataToCdArray(size_t *cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    unsigned char bytes[8] = {0};
    uint64_t size;
    int dim = computeDimension(r5, r4, r3, r2, r1);
    cd_values[0] = dim;
    cd_values[1] = dataType;    //0: FLOAT ; 1: DOUBLE ; 2,3,4,....: INTEGER....

    switch (dim) {
        case 1:
            size = (uint64_t) r1;
            SZ3::int64ToBytes_bigEndian(bytes, size);
            cd_values[2] = SZ3::bytesToInt32_bigEndian(bytes);
            cd_values[3] = SZ3::bytesToInt32_bigEndian(&bytes[4]);
            *cd_nelmts = 4;
            break;
        case 2:
            cd_values[2] = (unsigned int) r2;
            cd_values[3] = (unsigned int) r1;
            *cd_nelmts = 4;
            break;
        case 3:
            cd_values[2] = (unsigned int) r3;
            cd_values[3] = (unsigned int) r2;
            cd_values[4] = (unsigned int) r1;
            *cd_nelmts = 5;
            break;
        case 4:
            cd_values[2] = (unsigned int) r4;
            cd_values[3] = (unsigned int) r3;
            cd_values[4] = (unsigned int) r2;
            cd_values[5] = (unsigned int) r1;
            *cd_nelmts = 6;
            break;
        default:
            cd_values[2] = (unsigned int) r5;
            cd_values[3] = (unsigned int) r4;
            cd_values[4] = (unsigned int) r3;
            cd_values[5] = (unsigned int) r2;
            cd_values[6] = (unsigned int) r1;
            *cd_nelmts = 7;
    }
}

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]) {
    int result = 0; //0 means no-error-information-in-cd_values; 1 means cd_values contains error information
    int dimSize = cd_values[0];
    //printf("nc_nelmts = %d\n", cd_nelmts);
    switch (dimSize) {
        case 1:
            if (cd_nelmts > 4)
                result = 1;
            break;
        case 2:
            if (cd_nelmts > 4)
                result = 1;
            break;
        case 3:
            if (cd_nelmts > 5)
                result = 1;
            break;
        case 4:
            if (cd_nelmts > 6)
                result = 1;
            break;
        case 5:
            if (cd_nelmts > 7)
                result = 1;
            break;
    }
    return result;
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    size_t dataLength;
    if (r1 == 0) {
        dataLength = 0;
    } else if (r2 == 0) {
        dataLength = r1;
    } else if (r3 == 0) {
        dataLength = r1 * r2;
    } else if (r4 == 0) {
        dataLength = r1 * r2 * r3;
    } else if (r5 == 0) {
        dataLength = r1 * r2 * r3 * r4;
    } else {
        dataLength = r1 * r2 * r3 * r4 * r5;
    }
    return dataLength;
}

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    int dimension;
    if (r1 == 0) {
        dimension = 0;
    } else if (r2 == 0) {
        dimension = 1;
    } else if (r3 == 0) {
        dimension = 2;
    } else if (r4 == 0) {
        dimension = 3;
    } else if (r5 == 0) {
        dimension = 4;
    } else {
        dimension = 5;
    }
    return dimension;
}

void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1) {
    switch (dim) {
        case 1:
            dims[0] = r1;
            if (nbEle <= MAX_CHUNK_SIZE) //2^32-1
                chunk[0] = r1;
            else
                chunk[0] = 2147483648;//2^31
            break;
        case 2:
            dims[0] = r2;
            dims[1] = r1;
            if (nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r2;
                chunk[1] = r1;
            } else {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        case 3:
            dims[0] = r3;
            dims[1] = r2;
            dims[2] = r1;
            if (nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r3;
                chunk[1] = r2;
                chunk[2] = r1;
            } else {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        case 4:
            dims[0] = r4;
            dims[1] = r3;
            dims[2] = r2;
            dims[3] = r1;
            if (nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r4;
                chunk[1] = r3;
                chunk[2] = r2;
                chunk[3] = r1;
            } else {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        default:
            dims[0] = r5;
            dims[1] = r4;
            dims[2] = r3;
            dims[3] = r2;
            dims[4] = r1;
            if (nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r5;
                chunk[1] = r4;
                chunk[2] = r3;
                chunk[3] = r2;
                chunk[4] = r1;
            } else {
                printf("Error: size is too big!\n");
                exit(0);
            }
    }
}

//detect sys endian type
inline void detectSysEndianType() {
    //get sys endian type
    int x_temp = 1;
    char *y_temp = (char *) &x_temp;

    if (*y_temp == 1)
        sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        sysEndianType = BIG_ENDIAN_SYSTEM;
}

inline void symTransform_8bytes(unsigned char data[8]) {
    unsigned char tmp = data[0];
    data[0] = data[7];
    data[7] = tmp;

    tmp = data[1];
    data[1] = data[6];
    data[6] = tmp;

    tmp = data[2];
    data[2] = data[5];
    data[5] = tmp;

    tmp = data[3];
    data[3] = data[4];
    data[4] = tmp;
}

//the byte to input is in the big-endian format
inline double bytesToDouble(unsigned char *bytes) {
    ldouble buf;
    memcpy(buf.byte, bytes, 8);
    if (sysEndianType == LITTLE_ENDIAN_SYSTEM)
        symTransform_8bytes(buf.byte);
    return buf.value;
}

inline void doubleToBytes(unsigned char *b, double num) {
    ldouble buf;
    buf.value = num;
    memcpy(b, buf.byte, 8);
    if (sysEndianType == LITTLE_ENDIAN_SYSTEM)
        symTransform_8bytes(b);
}


/**
 * @brief		check dimension and correct it if needed
 * @return 	0 (didn't change dimension)
 * 					1 (dimension is changed)
 * 					2 (dimension is problematic)
 **/
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *correctedDimension) {
    int dimensionCorrected = 0;
    int dim = computeDimension(r5, r4, r3, r2, r1);
    correctedDimension[0] = r1;
    correctedDimension[1] = r2;
    correctedDimension[2] = r3;
    correctedDimension[3] = r4;
    correctedDimension[4] = r5;
    size_t *c = correctedDimension;
    if (dim == 1) {
        if (r1 < 1)
            return 2;
    } else if (dim == 2) {
        if (r2 == 1) {
            c[1] = 0;
            dimensionCorrected = 1;
        }
        if (r1 == 1) //remove this dimension
        {
            c[0] = c[1];
            c[1] = c[2];
            dimensionCorrected = 1;
        }
    } else if (dim == 3) {
        if (r3 == 1) {
            c[2] = 0;
            dimensionCorrected = 1;
        }
        if (r2 == 1) {
            c[1] = c[2];
            c[2] = c[3];
            dimensionCorrected = 1;
        }
        if (r1 == 1) {
            c[0] = c[1];
            c[1] = c[2];
            c[2] = c[3];
            dimensionCorrected = 1;
        }
    } else if (dim == 4) {
        if (r4 == 1) {
            c[3] = 0;
            dimensionCorrected = 1;
        }
        if (r3 == 1) {
            c[2] = c[3];
            c[3] = c[4];
            dimensionCorrected = 1;
        }
        if (r2 == 1) {
            c[1] = c[2];
            c[2] = c[3];
            c[3] = c[4];
            dimensionCorrected = 1;
        }
        if (r1 == 1) {
            c[0] = c[1];
            c[1] = c[2];
            c[2] = c[3];
            c[3] = c[4];
            dimensionCorrected = 1;
        }
    } else if (dim == 5) {
        if (r5 == 1) {
            c[4] = 0;
            dimensionCorrected = 1;
        }
        if (r4 == 1) {
            c[3] = c[4];
            c[4] = 0;
            dimensionCorrected = 1;
        }
        if (r3 == 1) {
            c[2] = c[3];
            c[3] = c[4];
            c[4] = 0;
            dimensionCorrected = 1;
        }
        if (r2 == 1) {
            c[1] = c[2];
            c[2] = c[3];
            c[3] = c[4];
            c[4] = 0;
            dimensionCorrected = 1;
        }
        if (r1 == 1) {
            c[0] = c[1];
            c[1] = c[2];
            c[2] = c[3];
            c[3] = c[4];
            c[4] = 0;
            dimensionCorrected = 1;
        }
    }

    return dimensionCorrected;

}

inline void longToBytes_bigEndian(unsigned char *b, uint64_t num) {
    b[0] = (unsigned char) (num >> 56);
    b[1] = (unsigned char) (num >> 48);
    b[2] = (unsigned char) (num >> 40);
    b[3] = (unsigned char) (num >> 32);
    b[4] = (unsigned char) (num >> 24);
    b[5] = (unsigned char) (num >> 16);
    b[6] = (unsigned char) (num >> 8);
    b[7] = (unsigned char) (num);
//	if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_8bytes(*b);
}

inline int bytesToInt_bigEndian(unsigned char *bytes) {
    int temp = 0;
    int res = 0;

    res <<= 8;
    temp = bytes[0] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[1] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[2] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[3] & 0xff;
    res |= temp;

    return res;
}

/**
 * @endianType: refers to the endian_type of unsigned char* b.
 * */
inline int64_t bytesToLong_bigEndian(unsigned char *b) {
    int64_t temp = 0;
    int64_t res = 0;

    res <<= 8;
    temp = b[0] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[1] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[2] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[3] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[4] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[5] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[6] & 0xff;
    res |= temp;

    res <<= 8;
    temp = b[7] & 0xff;
    res |= temp;

    return res;
}
