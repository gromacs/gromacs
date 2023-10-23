#ifndef _CUNIT_ARRAY_H
#define _CUNIT_ARRAY_H

#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"

#ifdef __cplusplus
extern "C" {
#endif

void CU_ASSERT_EQUAL_ARRAY_INT(int* actual, int* expected, int count);
void CU_ASSERT_EQUAL_ARRAY_BYTE(unsigned char* actual, unsigned char* expected, int count);
void CU_ASSERT_EQUAL_ARRAY_FLOAT(float* actual, float* expected, int count, double granularity);
void CU_ASSERT_EQUAL_ARRAY_DOUBLE(double* actual, double* expected, int count, double granularity);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _CUNIT_ARRAY_H  ----- */
