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
#include <SZ3/utils/ByteUtil.hpp>

#include "hdf5.h"
#include "H5Cpp.h"
#include "H5Z_SZ3.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"

#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1
#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1
int sysEndianType = LITTLE_ENDIAN_SYSTEM;
int dataEndianType = LITTLE_ENDIAN_DATA;

#define DATASET "testdata_compressed"

int MAX_CHUNK_SIZE = INT_MAX;

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

template<typename T>
void process_data(const SZ3::Config &conf,
                  const char *oriFilePath,
                  int dataEndianType,
                  hid_t fid,
                  hid_t sid,
                  hid_t cpid,
                  const char *datasetName,
                  hid_t h5TypeLE,
                  hid_t h5TypeBE) {
    T *data = new T[conf.num];
    SZ3::readfile(oriFilePath, conf.num, data);
    
    std::cout << "original data = ";
    for (int i = 0; i < 20; ++i)
        std::cout << data[i] << " ";
    std::cout << "....\n";
    
    hid_t dataset;
    hid_t h5Type = (dataEndianType == LITTLE_ENDIAN_DATA) ? h5TypeLE : h5TypeBE;
    
    if ((dataset = H5Dcreate(fid, datasetName, h5Type, sid, H5P_DEFAULT, cpid, H5P_DEFAULT)) < 0) {
        std::cerr << "Error in H5Dcreate\n";
        delete[] data;
        exit(EXIT_FAILURE);
    }
    if (H5Dwrite(dataset, h5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
        std::cerr << "Error in H5Dwrite\n";
        delete[] data;
        exit(EXIT_FAILURE);
    }
    
    delete[] data;
    if (H5Dclose(dataset) < 0) {
        std::cerr << "Error in H5Dclose\n";
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[]) {
    
    char oriFilePath[640], outputFilePath[640];
    
    if (argc < 3) {
        printf("Test case: sz3ToHDF5 [dataType] [srcFilePath] [dimension sizes...]\n");
        printf("Example1 : sz3ToHDF5 -f testdata/x86/testfloat_8_8_128.dat 8 8 128\n");
        printf("Example 2: sz3ToHDF5 -i32 testdata/x86/testint32_8x8x8.dat 8 8 8\n");
        exit(0);
    }
    
    std::map<std::string, int> dataTypeMap = {
        {"-f", SZ_FLOAT}, {"-d", SZ_DOUBLE}, {"-i8", SZ_INT8}, {"-u8", SZ_UINT8},
        {"-i16", SZ_INT16}, {"-u16", SZ_UINT16}, {"-i32", SZ_INT32}, {"-u32", SZ_UINT32},
        {"-i64", SZ_INT64}, {"-u64", SZ_UINT64}
    };
    
    int dataType = 0;
    auto it = dataTypeMap.find(argv[1]);
    if (it != dataTypeMap.end()) {
        dataType = it->second;
    } else {
        std::cerr << "Error: unknown data type in sz3ToHDF5.c!\n";
        return 0;
    }
    
    snprintf(oriFilePath, 640, "%s", argv[2]);
    
    std::vector<int> dimensions;
    for (int i = 3; i < argc && i < 8; ++i) {
        dimensions.push_back(std::atoi(argv[i]));
    }
    std::reverse(dimensions.begin(), dimensions.end()); // slowest to fastest
    
    snprintf(outputFilePath, 640, "%s.sz3.h5", oriFilePath);
    
    hid_t sid, cpid, fid;
    /* create HDF5 file */
    if (0 > (fid = H5Fcreate(outputFilePath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT))) {
        printf("Error in H5Pcreate");
        exit(0);
    }
    
    
    //set up SZ configuration
    SZ3::Config conf;
    // use config.loadcfg("path_to_sz3_conf") to load configuration from a file if needed
    // setting up data related attributes (data type, dims, etc.) is not necessary, as they will be updated in H5Z_sz3_set_local automatically
    conf.setDims(dimensions.begin(), dimensions.end());
    // Set compression related attributes here
    conf.cmprAlgo = SZ3::ALGO_BIOMD;
    
    //save conf to cd_values
    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts);
    auto buffer = (unsigned char *) (cd_values.data());
    conf.save(buffer);
    // conf.print();
    
    std::vector<hsize_t> hdims(conf.dims.begin(), conf.dims.end());
    /*Create dataspace. Setting maximum size */
    if (0 > (sid = H5Screate_simple(conf.N, hdims.data(), NULL))) {
        printf("Error in H5Screate_simple");
        exit(0);
    }
    /* setup dataset creation properties */
    if (0 > (cpid = H5Pcreate(H5P_DATASET_CREATE))) {
        printf("Error in H5Pcreate");
        exit(0);
    }
    
    /* Add the SZ compression filter */
    if (0 > H5Pset_filter(cpid, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data())) {
        printf("Error in H5Pcreate");
        exit(0);
    }
    herr_t ret = H5Zregister(H5PLget_plugin_info());
    if (ret < 0) {
        printf("Error in H5Zregister");
        exit(0);
    }
    htri_t avail = H5Zfilter_avail(H5Z_FILTER_SZ3);
    if (avail) {
        unsigned filter_config;
        auto status = H5Zget_filter_info(H5Z_FILTER_SZ3, &filter_config);
        
        if (filter_config & H5Z_FILTER_CONFIG_ENCODE_ENABLED)
            printf("sz filter is available for encoding and decoding.\n");
    }
    
    /*  set the chunk size*/
    std::vector<hsize_t> hchunk(hdims);
//    hchunk[0] = 10;
    if (0 > H5Pset_chunk(cpid, conf.N, hchunk.data())) {
        printf("Error in H5Pcreate");
        exit(0);
    }
    
    {//This is an example to get/set SZ configuration from HDF5 file
        SZ3::Config conf1;
        get_SZ3_conf_from_H5(cpid, conf1);
//        conf1.absErrorBound = 1;
        set_SZ3_conf_to_H5(cpid, conf1);
    }
    printf("....Writing SZ compressed data.............\n");
    
    switch (dataType) {
        case SZ_FLOAT:process_data<float>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_IEEE_F32LE, H5T_IEEE_F32BE);
            break;
        case SZ_DOUBLE:process_data<double>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_IEEE_F64LE, H5T_IEEE_F64BE);
            break;
        case SZ_INT8:process_data<int8_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_I8LE, H5T_STD_I8BE);
            break;
        case SZ_UINT8:process_data<uint8_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_U8LE, H5T_STD_U8BE);
            break;
        case SZ_INT16:process_data<int16_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_I16LE, H5T_STD_I16BE);
            break;
        case SZ_UINT16:process_data<uint16_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_U16LE, H5T_STD_U16BE);
            break;
        case SZ_INT32:process_data<int32_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_I32LE, H5T_STD_I32BE);
            break;
        case SZ_UINT32:process_data<uint32_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_U32LE, H5T_STD_U32BE);
            break;
        case SZ_INT64:process_data<int64_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_I64LE, H5T_STD_I64BE);
            break;
        case SZ_UINT64:process_data<uint64_t>(conf, oriFilePath, dataEndianType, fid, sid, cpid, DATASET, H5T_STD_U64LE, H5T_STD_U64BE);
            break;
        default:std::cerr << "Error: Unknown data type\n";
            exit(EXIT_FAILURE);
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
    printf("Output hdf5 file: %s\n", outputFilePath);
    ret = H5Zunregister(H5Z_FILTER_SZ3);
    if (ret < 0) return -1;
    H5close();
    return 0;
}
