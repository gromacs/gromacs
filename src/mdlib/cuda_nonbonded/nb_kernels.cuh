/*!
 *  Generate kernels for the different electrostatics types:
 *  Cut-off, Reaction-Field, and Ewald/PME.
 *  (No include fence as it can be included multiple times.)
 */

/* Cut-Off */
#define EL_CUTOFF
#define FUNCTION_NAME(x, y) x##_cutoff_##y
#include "nb_kernel_1.cuh"
#include "nb_kernel_2.cuh"
#undef EL_CUTOFF
#undef FUNCTION_NAME

/* Reaction-Field */
#define EL_RF
#define FUNCTION_NAME(x, y) x##_RF_##y
#include "nb_kernel_1.cuh"
#include "nb_kernel_2.cuh"
#undef EL_RF
#undef FUNCTION_NAME

/* Ewald */
#define EL_EWALD
#define FUNCTION_NAME(x, y) x##_ewald_##y
#include "nb_kernel_1.cuh"
#include "nb_kernel_2.cuh"
#undef EL_EWALD
#undef FUNCTION_NAME
