/**
 *  @file H5Z_SZ.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the H5Z_SZ.c.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _H5Z_SZ
#define _H5Z_SZ

#include <stdio.h>
#include <hdf5.h>
#include "sz.h"

#define H5Z_FILTER_SZ 32017
#define MAX_CHUNK_SIZE 4294967295 //2^32-1
static hid_t H5Z_SZ_ERRCLASS = -1;

#ifdef __cplusplus
extern "C" {
#endif


extern int load_conffile_flag;
extern int init_sz_flag;

extern char cfgFile[256];

/* convenience macro to handle errors */
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
	H5Epush(H5E_DEFAULT,__FILE__,_funcname_,__LINE__, \
		H5Z_SZ_ERRCLASS,MAJ,MIN,MSG);                \
	retval = RET;                                     \
	goto done;                                        \
} while(0)

int H5Z_SZ_Init(char* cfgFile);
int H5Z_SZ_Init_Params(sz_params *params);
sz_params* H5Z_SZ_Init_Default();
int H5Z_SZ_Finalize();

void SZ_refreshDimForCdArray(int dataType, size_t old_cd_nelmts, unsigned int *old_cd_values, size_t* new_cd_nelmts, unsigned int **new_cd_values, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1);
void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
int* error_bound_mode, double* abs_error, double* rel_error, double* pw_rel_error, double* psnr);

void SZ_errConfigToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int error_bound_mode, double abs_error, double rel_error, double pw_rel_error, double psnr);

int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]);
static size_t H5Z_filter_sz(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf);
static herr_t H5Z_sz_set_local(hid_t dcpl_id, hid_t type_id, hid_t space_id);


void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _H5Z_SZ_metadata  ----- */
