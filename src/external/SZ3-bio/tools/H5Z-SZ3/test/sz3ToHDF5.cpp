/**
 *  @file szToHDF5.c
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief This is an example of using compression interface (HDF5)
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "hdf5.h"
#include "H5Z_SZ3.hpp"
#include "SZ3/utils/FileUtil.hpp"

int sysEndianType = LITTLE_ENDIAN_SYSTEM;
int dataEndianType = LITTLE_ENDIAN_DATA;

#define DATASET "testdata_compressed"


using namespace SZ3;

int main(int argc, char *argv[]) {

    //(void) helper fn to detect system endian type
    //detectSysEndianType();
    //by default sysEndianType and dataEndianType are little endian, can set them manually here
    //dataEndianType = BIG_ENDIAN_DATA;

    size_t r5 = 0, r4 = 0, r3 = 0, r2 = 0, r1 = 0;
    int cmp_algo, interp_algo; //select compression and interpolation for SZ3
    char outDir[640], oriFilePath[640], outputFilePath[640];
    size_t cd_nelmts=0, nbEle;
    unsigned int *cd_values = NULL;
    //unsigned int cd_values[7];

    herr_t status;
    htri_t avail;
    unsigned filter_config;

    hid_t sid, idsid, cpid, fid;

    if (argc < 3) {
        printf("Test case: sz3ToHDF5 [dataType] [srcFilePath] [dimension sizes...]\n");
        printf("Example1 : sz3ToHDF5 -f testdata/x86/testfloat_8_8_128.dat 8 8 128\n");
        printf("Example 2: sz3ToHDF5 -i32 testdata/x86/testint32_8x8x8.dat 8 8 8\n");
        exit(0);
    }

    //printf("config file = %s\n", argv[2]);

    int dataType = 0;
    if (strcmp(argv[1], "-f") == 0)
        dataType = SZ_FLOAT;
    else if (strcmp(argv[1], "-d") == 0)
        dataType = SZ_DOUBLE;
    else if (strcmp(argv[1], "-i8") == 0)
        dataType = SZ_INT8;
    else if (strcmp(argv[1], "-u8") == 0)
        dataType = SZ_UINT8;
    else if (strcmp(argv[1], "-i16") == 0)
        dataType = SZ_INT16;
    else if (strcmp(argv[1], "-u16") == 0)
        dataType = SZ_UINT16;
    else if (strcmp(argv[1], "-i32") == 0)
        dataType = SZ_INT32;
    else if (strcmp(argv[1], "-u32") == 0)
        dataType = SZ_UINT32;
    else if (strcmp(argv[1], "-i64") == 0)
        dataType = SZ_INT64;
    else if (strcmp(argv[1], "-u64") == 0)
        dataType = SZ_UINT64;
    else {
        printf("Error: unknown data type in sz3ToHDF5.c!\n");
        exit(0);
    }

    printf("DTYPE: %i", dataType);
    snprintf(oriFilePath, 640, "%s", argv[2]);
    if (argc >= 4) {
        r1 = atoi(argv[3]); //8
    }
    if (argc >= 5) {
        r2 = atoi(argv[4]); //8
    }
    if (argc >= 6) {
        r3 = atoi(argv[5]); //128
    }
    if (argc >= 7) {
        r4 = atoi(argv[6]);
    }
    if (argc >= 8) {
        r5 = atoi(argv[7]);
    }

    //read in compression and interp algo
    //for testing set these here as defaults in config
    cmp_algo = 1;
    interp_algo = 1;

    //printf("cfgFile=%s\n", cfgFile);
    snprintf(outputFilePath, 640, "%s.sz3.h5", oriFilePath);

//	printf("argv[1]=%s, dataType=%d\n", argv[1], dataType);
    nbEle = computeDataLength(r5, r4, r3, r2, r1);

//	printf("nbEle=%u\n", nbEle);

    //Create cd_values
    printf("Dimension sizes: n5=%u, n4=%u, n3=%u, n2=%u, n1=%u\n", r5, r4, r3, r2, r1);
    int mode = 0; //0: ABS, 1: REL, ...
    SZ_errConfigToCdArray(&cd_nelmts, &cd_values, mode, 0.001, 0.001, 0,
                          0); //SZ_FLOAT or SZ_DOUBLE or SZ_INT 100x500x500 : 0, 0, 100, 500, 500, ABS, REL (0.01, 0.01*(max-min), PW_REL (0.01, 5, 6, 7, 8, 9 --> 5*0.01, 6*0.01, ...), PSNR (mean squared error)).
    //load_conffile_flag = 0;
    //										 REL
    //SZ_metaDataErrToCdArray(&cd_nelmts, &cd_values, dataType, r5, r4, r3, r2, r1, 1, 0.01, 0.01, 0, 0); //SZ_FLOAT or SZ_DOUBLE or SZ_INT 100x500x500 : 0, 0, 100, 500, 500, ABS, REL (0.01, 0.01*(max-min), PW_REL (0.01, 5, 6, 7, 8, 9 --> 5*0.01, 6*0.01, ...), PSNR (mean squared error)).
    /*cd_nelmts = 5;
    cd_values[0] = 3;
    cd_values[1] = 0;
    cd_values[2] = 128;
    cd_values[3] = 8;
    cd_values[4] = 8;
    cd_values[5] = 0;
    cd_values[6] = 0;*/

    int i = 0;
//	for(i=0;i<cd_nelmts;i++)
//		printf("cd_values[%d]=%u\n", i, cd_values[i]);

    //compute dimension
    int dim = computeDimension(r5, r4, r3, r2, r1);

    hsize_t dims[5] = {0, 0, 0, 0, 0}, chunk[5] = {0, 0, 0, 0, 0};
    init_dims_chunk(dim, dims, chunk, nbEle, r5, r4, r3, r2, r1);

    /* create HDF5 file */
    if (0 > (fid = H5Fcreate(outputFilePath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT))) {
        printf("Error in H5Pcreate");
        exit(0);
    }

    /*Create dataspace. Setting maximum size */
    if (0 > (sid = H5Screate_simple(dim, dims, NULL))) {
        printf("Error in H5Screate_simple");
        exit(0);
    }

    /* setup dataset creation properties */
    if (0 > (cpid = H5Pcreate(H5P_DATASET_CREATE))) {
        printf("Error in H5Pcreate");
        exit(0);
    }

    /* Add the SZ compression filter and set the chunk size */
    if (0 > H5Pset_filter(cpid, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values)) {
        printf("Error in H5Pcreate");
        exit(0);
    }
    avail = H5Zfilter_avail(H5Z_FILTER_SZ3);
    if (avail) {
        status = H5Zget_filter_info(H5Z_FILTER_SZ3, &filter_config);

        if (filter_config & H5Z_FILTER_CONFIG_ENCODE_ENABLED)
            printf("sz filter is available for encoding and decoding.\n");
    }
    if (0 > H5Pset_chunk(cpid, dim, chunk)) {
        printf("Error in H5Pcreate");
        exit(0);
    }

    //Initialize the configuration for SZ
    //You can also use the global variable conf_params to set the configuration for sz without cfgFile.
    //Example of setting an absolute error bound:
    //			sz_params* params = H5Z_SZ_Init_Default();
    //			params->errorBoundMode = ABS;
    //			params->absErrBound = 1E-4;

    //H5Z_SZ_Init(cfgFile);

    printf("....Writing SZ compressed data.............\n");

    if (dataType == SZ_FLOAT) {
        float *data = new float[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%f ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F32LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F32BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        };
    } else if (dataType == SZ_DOUBLE) {
        double *data = new double[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%f ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F64LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_IEEE_F64BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_IEEE_F64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        };
    } else if (dataType == SZ_INT8) {
        int8_t *data = new int8_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I8LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I8BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_UINT8) {
        uint8_t *data = new uint8_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U8LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U8BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U8BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_INT16) {

        int16_t *data = new int16_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I16LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I16BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_UINT16) {
        uint16_t *data = new uint16_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U16LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U16BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U16BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_INT32) {
        //printf("%i \t %i\n", sizeof(int), sizeof(int32_t));
        int32_t *data = new int32_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I32LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I32BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_UINT32) {
        uint32_t *data = new uint32_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%d ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U32LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U32BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_INT64) {
        int64_t *data = new int64_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%ld ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I64LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_I64BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_I64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else if (dataType == SZ_UINT64) {
        uint64_t *data = new uint64_t[nbEle];
        readfile(oriFilePath, nbEle, data);

        printf("original data = ");
        for (i = 0; i < 20; i++)
            printf("%ld ", data[i]);
        printf("....\n");

        if (dataEndianType == LITTLE_ENDIAN_DATA) {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U64LE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        } else //BIG_ENDIAN_DATA
        {
            if (0 > (idsid = H5Dcreate(fid, DATASET, H5T_STD_U64BE, sid, H5P_DEFAULT, cpid, H5P_DEFAULT))) {
                printf("Error in H5Dcreate");
                exit(0);
            }
            if (0 > H5Dwrite(idsid, H5T_STD_U64BE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) {
                printf("Error in H5Dwrite");
                exit(0);
            }
        }
        delete[] data;
        if (0 > H5Dclose(idsid)) {
            printf("Error in H5Dclose");
            exit(0);
        }
    } else {
        printf("Error: unknown data type in sz3ToHDF5.cpp!\n");
        exit(0);
    }

    /*Close and release resources*/
    if (0 > H5Sclose(sid)) {
        printf("Error in H5Sclose");
        exit(0);
    }
    if (0 > H5Pclose(cpid)) {
        printf("Error in H5Pclose");
        exit(0);
    }
    if (0 > H5Fclose(fid)) {
        printf("Error in H5Fclose");
        exit(0);
    }
    free(cd_values);
    printf("Output hdf5 file: %s\n", outputFilePath);
    herr_t ret = H5Zunregister(H5Z_FILTER_SZ3);
    if (ret < 0) return -1;
    H5close();
    return 0;
}
