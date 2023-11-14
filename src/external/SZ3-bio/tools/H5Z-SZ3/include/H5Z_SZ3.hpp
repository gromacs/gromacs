//
// Created by arham23 on 2/8/22.
//

#ifndef SZ3_H5Z_SZ3_H
#define SZ3_H5Z_SZ3_H

#define H5Z_FILTER_SZ3 32024
#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9

#include "hdf5.h"
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <cstdint>




#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1
#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1

#ifdef __cplusplus
extern "C" {
#endif

extern hid_t H5Z_SZ_ERRCLASS;
#define ERROR(FNAME)                                              \
do {                                                              \
    int _errno = errno;                                           \
    fprintf(stderr, #FNAME " failed at line %d, errno=%d (%s)\n", \
        __LINE__, _errno, _errno?strerror(_errno):"ok");          \
    return 1;                                                     \
} while(0)

#define H5Z_SZ_PUSH_AND_GOTO(MAJ, MIN, RET, MSG)     \
do                                                    \
{                                                     \
	H5Epush(H5E_DEFAULT,__FILE__,_funcname_,__LINE__,H5Z_SZ_ERRCLASS,MAJ,MIN,MSG); \
	return RET;                                     \
} while(0)
extern int sysEndianType;
extern int dataEndianType;

void SZ_refreshDimForCdArray(int dataType, size_t old_cd_nelmts, unsigned int *old_cd_values, size_t* new_cd_nelmts, unsigned int **new_cd_values, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void SZ_errConfigToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int error_bound_mode, double abs_error, double rel_error, double l2normErrorBound, double psnr);

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id);

static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf);

void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1);

void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
                            int* error_bound_mode, double* abs_error, double* rel_error, double* l2norm_error, double* psnr);

void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]);

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

double bytesToDouble(unsigned char* bytes);
void doubleToBytes(unsigned char *b, double num);

void longToBytes_bigEndian(unsigned char *b, uint64_t num) ;

int bytesToInt_bigEndian(unsigned char* bytes);
int64_t bytesToLong_bigEndian(unsigned char* b);

void detectSysEndianType();
void symTransform_8bytes(unsigned char data[8]);

int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension);

#ifdef __cplusplus
}
#endif

#endif //SZ3_H5Z_SZ3_H
