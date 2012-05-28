
#ifndef _gmxDetectCpu_h_
#define _gmxDetectCpu_h_

typedef enum 
{
    GMX_DETECTCPU_VENDOR_UNKNOWN = 0,
    GMX_DETECTCPU_VENDOR_INTEL,
    GMX_DETECTCPU_VENDOR_AMD,
    GMX_DETECTCPU_NVENDORS
} 
gmxDetectCpuVendorId_t;
    
extern const char *
gmxDetectCpuVendorIdString[GMX_DETECTCPU_NVENDORS];




/* Always add entries to the end of this list, just before the last NFEATURES line.
 * To keep the length of this list reasonable, we only add flags referring to
 * features that we actually might have to check/use in Gromacs - feel free to add more.
*
 * AMD and Intel are unfortunately gradually diverging, so while we can use the
 * same type of intrinsic instruction functions in the source, the resulting binary
 * is frequently not compatible starting from AVX. 
 */
typedef enum
{
    GMX_DETECTCPU_FEATURE_CANNOTDETECT = 0,  /* Flag set if we could not detect on this CPU  */
    GMX_DETECTCPU_FEATURE_X86_HTT,           /* Hyperthreading technology                    */        
    GMX_DETECTCPU_FEATURE_X86_SSE2,          /* SSE 2                                        */
    GMX_DETECTCPU_FEATURE_X86_SSE4_1,        /* SSE 4.1                                      */
    GMX_DETECTCPU_FEATURE_X86_RDRAND,        /* RDRAND high-quality hardware random numbers  */
    GMX_DETECTCPU_FEATURE_X86_AES,           /* x86 advanced encryption standard accel.      */
    GMX_DETECTCPU_FEATURE_X86_AVX,           /* Advanced vector extensions                   */
    GMX_DETECTCPU_FEATURE_X86_FMA,           /* Fused-multiply add support (mainly for AVX)  */
    GMX_DETECTCPU_FEATURE_X86_FMA4,          /* 4-operand FMA, only on AMD for now           */
    GMX_DETECTCPU_FEATURE_X86_XOP,           /* AMD extended instructions, only AMD for now  */
    GMX_DETECTCPU_FEATURE_X86_AVX2,          /* AVX2 including gather support (not used yet) */
    GMX_DETECTCPU_NFEATURES
}
gmxDetectCpuFeature_t;

extern const char *
gmxDetectCpuFeatureString[GMX_DETECTCPU_NFEATURES];


/* Currently supported acceleration instruction sets, intrinsics or other similar combinations
 * in Gromacs. There is not always a 1-to-1 correspondence with feature flags; on some AMD
 * hardware we prefer to use 128bit AVX instructions (although 256-bit ones could be executed),
 * and we still havent written the AVX2 kernels.
 */
typedef enum
{
    GMX_DETECTCPU_ACCELERATION_NONE = 0,
    GMX_DETECTCPU_ACCELERATION_X86_SSE2,
    GMX_DETECTCPU_ACCELERATION_X86_SSE4_1,
    GMX_DETECTCPU_ACCELERATION_X86_AVX_128_FMA,
    GMX_DETECTCPU_ACCELERATION_X86_AVX_256,
    GMX_DETECTCPU_NACCELERATIONS
}
gmxDetectCpuAcceleration_t;


extern const char *
gmxDetectCpuAccelerationString[GMX_DETECTCPU_NACCELERATIONS];



#define GMX_DETECTCPU_STRLEN  64

typedef struct
{
    gmxDetectCpuVendorId_t     vendorId;
    char                       brand[GMX_DETECTCPU_STRLEN];
    int                        family;
    int                        model;
    int                        stepping;
    
    char                       feature[GMX_DETECTCPU_NFEATURES];
    /*
     Not used right now...
    int                        physicalCoresPerComputeUnit;
    int                        logicalCoresPerComputeUnit;
    int                        physicalCoresPerPackage;
    int                        logicalCoresPerPackage;
     */
}
gmxDetectCpu_t;



/* Fill the entire structure by using CPU detection instructions.
 * Return 0 on success, 1 if something bad happened. 
 */
int
gmxDetectCpu                   (gmxDetectCpu_t *              data);


/* Formats a text string (up to n characters) from the data structure.
 * The output will have max 80 chars between newline characters.
 */
int
gmxDetectCpuFormatString       (gmxDetectCpu_t                data,
                                char *                        s,
                                int                           n);


/* Suggests a suitable gromacs acceleration based on the support in the
 * hardware. 
 */
int
gmxDetectCpuSuggestAcceleration  (gmxDetectCpu_t                data,
                                  gmxDetectCpuAcceleration_t *  acc);

/* Check if this binary was compiled with the same acceleration as we
 * would suggest for the current hardware.
 */
int
gmxDetectCpuCheckAcceleration    (gmxDetectCpu_t                data,
                                  FILE *                        log);

 


#endif /* _gmxDetectCpu_h_ */
