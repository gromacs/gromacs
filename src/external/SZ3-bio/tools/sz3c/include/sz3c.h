//
// Created by Kai Zhao on 10/27/22.
//

#ifndef SZ3_SZ3C_H
#define SZ3_SZ3C_H

#include <stdio.h>
/** Begin errorbound mode in SZ2 (defines.h) **/
#define ABS 0
#define REL 1
#define VR_REL 1  //alternative name to REL
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4
#define NORM 5

#define PW_REL 10
#define ABS_AND_PW_REL 11
#define ABS_OR_PW_REL 12
#define REL_AND_PW_REL 13
#define REL_OR_PW_REL 14
/** End errorbound mode in SZ2 (defines.h) **/

/** Begin dataType in SZ2 (defines.h) **/
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
/** End dataType in SZ2 (defines.h) **/

#ifdef __cplusplus
extern "C" {
#endif
unsigned char *SZ_compress_args(int dataType, void *data, size_t *outSize,
                                int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio,
                                size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength,
                    size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#ifdef __cplusplus
}
#endif

#endif //SZ3_SZ3C_H
