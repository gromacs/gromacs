#if defined (__cplusplus)
extern "C" {
#endif

#ifndef SZ_API_H
#define SZ_API_H

/* =====   SZLIB_API : control library symbols visibility   ===== */
#ifndef SZLIB_VISIBILITY
#  if defined(__GNUC__) && (__GNUC__ >= 4)
#    define SZLIB_VISIBILITY __attribute__ ((visibility ("default")))
#  else
#    define SZLIB_VISIBILITY
#  endif
#endif
#if defined(SZ_DLL_EXPORT) && (SZ_DLL_EXPORT==1)
#  define SZLIB_API __declspec(dllexport) SZLIB_VISIBILITY
#elif defined(SZ_DLL_IMPORT) && (SZ_DLL_IMPORT==1)
#  define SZLIB_API __declspec(dllimport) SZLIB_VISIBILITY /* It isn't required but allows to generate better code, saving a function pointer load from the IAT and an indirect jump.*/
#else
#  define SZLIB_API SZLIB_VISIBILITY
#endif

#include "defines.h"
#include "ByteToolkit.h"

/* array meta data and compression parameters for SZ_Init_Params() */
typedef struct sz_params
{
    int dataType;
    unsigned int max_quant_intervals; //max number of quantization intervals for quantization
    unsigned int quantization_intervals;
    unsigned int maxRangeRadius;
    int sol_ID;// it's SZ or SZ_Transpose, unless the setting is PASTRI compression mode (./configure --enable-pastri)
    int losslessCompressor;
    int sampleDistance; //2 bytes
    float predThreshold;  // 2 bytes
    int szMode; //* 0 (best speed) or 1 (better compression with Zstd/Gzip) or 3 temporal-dimension based compression
    int gzipMode; //* four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION
    int  errorBoundMode; //4bits (0.5byte), //ABS, REL, ABS_AND_REL, or ABS_OR_REL, PSNR, or PW_REL, PSNR
    double absErrBound; //absolute error bound
    double relBoundRatio; //value range based relative error bound ratio
    double psnr; //PSNR
    double normErr;
    double pw_relBoundRatio; //point-wise relative error bound
    int segment_size; //only used for 2D/3D data compression with pw_relBoundRatio (deprecated)
    int pwr_type; //only used for 2D/3D data compression with pw_relBoundRatio

    int protectValueRange; //0 or 1
    float fmin, fmax;
    double dmin, dmax;

    int snapshotCmprStep; //perform single-snapshot-based compression if time_step == snapshotCmprStep
    int predictionMode;

    int accelerate_pw_rel_compression;
    int plus_bits;

    int randomAccess;
    int withRegression;

} sz_params;

//-------------------key global variables--------------
extern int dataEndianType; //*endian type of the data read from disk
extern int sysEndianType; //*sysEndianType is actually set automatically.

extern sz_params *confparams_cpr;
extern sz_params *confparams_dec;

int SZ_Init(const char *configFilePath);

int SZ_Init_Params(sz_params *params);

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension);

unsigned char *SZ_compress(int dataType, void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound,
double relBoundRatio, double pwrBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
size_t SZ_decompress_args(int dataType, unsigned char *bytes, size_t byteLength, void* decompressed_array, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void SZ_Finalize();

#endif  /* SZ_API_H */

#if defined (__cplusplus)
}
#endif
