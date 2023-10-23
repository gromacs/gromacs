#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"

#include <stdio.h>  // for printf

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_computeRangeSize_double(void)
{
	double valueRangeSize, medianValue;
	double data[7] = {1,2,3,4,5,6,7};
	computeRangeSize_double(data, 7, &valueRangeSize, &medianValue);
	CU_ASSERT_DOUBLE_EQUAL(valueRangeSize, 6, 1E-6);
	CU_ASSERT_DOUBLE_EQUAL(medianValue, 4, 1E-6);
}

void test_computeRangeSize_float(void)
{
	float valueRangeSize, medianValue;
	float data[7] = {1,2,3,4,5,6,7};
	computeRangeSize_float(data, 7, &valueRangeSize, &medianValue);
	CU_ASSERT_DOUBLE_EQUAL(valueRangeSize, 6, 1E-6);
	CU_ASSERT_DOUBLE_EQUAL(medianValue, 4, 1E-6);	
}


void test_computeRangeSize_double_subblock(void)
{
	//TODO
}

void test_computeRangeSize_float_subblock(void)
{
	//TOOD
}

void test_min_d(void)
{
	CU_ASSERT_DOUBLE_EQUAL(min_d(1,2),1,1E-6);
	CU_ASSERT_DOUBLE_EQUAL(min_d(0,1),0,1E-6);
}

void test_max_d(void)
{
	CU_ASSERT_DOUBLE_EQUAL(max_d(1,2),2,1E-6);
	CU_ASSERT_DOUBLE_EQUAL(max_d(0,1),1,1E-6);	
}

void test_min_f(void)
{
	CU_ASSERT_DOUBLE_EQUAL(min_f(1,2),1,1E-6);
	CU_ASSERT_DOUBLE_EQUAL(min_f(0,1),0,1E-6);		
}

void test_max_f(void)
{
	CU_ASSERT_DOUBLE_EQUAL(max_f(1,2),2,1E-6);
	CU_ASSERT_DOUBLE_EQUAL(max_f(0,1),1,1E-6);		
}

void test_getRealPrecision_double(void)
{
	int status;
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, ABS, 0.01, 0.01, &status), 0.01, 1E-6);
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, REL, 0.01, 0.01, &status), 1, 1E-6);	
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, ABS_AND_REL, 0.01, 0.01, &status), 0.01, 1E-6);
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, ABS_OR_REL, 0.01, 0.01, &status), 1, 1E-6);
}

void test_getRealPrecision_float(void)
{
	int status;
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, ABS, 0.01, 0.01, &status), 0.01, 1E-6);
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, REL, 0.01, 0.01, &status), 1, 1E-6);	
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, ABS_AND_REL, 0.01, 0.01, &status), 0.01, 1E-6);
	CU_ASSERT_DOUBLE_EQUAL(getRealPrecision_double(100, ABS_OR_REL, 0.01, 0.01, &status), 1, 1E-6);	
}

void test_symTransform_8bytes(void)
{
	unsigned char bytes[8] = {1,2,3,4,5,6,7,8};
	symTransform_8bytes(bytes);
	unsigned char expected[8] = {8,7,6,5,4,3,2,1};
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytes, expected, 8);
}

void test_symTransform_2bytes(void)
{
	unsigned char bytes[2] = {1,2};
	symTransform_2bytes(bytes);
	unsigned char expected[2] = {2,1};
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytes, expected, 2);	
}

void test_symTransform_4bytes(void)
{
	unsigned char bytes[4] = {1,2,3,4};
	symTransform_4bytes(bytes);
	unsigned char expected[4] = {4,3,2,1};
	CU_ASSERT_EQUAL_ARRAY_BYTE(bytes, expected, 4);	
}

void test_compressSingleFloatValue(void)
{
	//TODO
}

void test_compressSingleDoubleValue(void)
{
	//TODO
}
void test_compIdenticalLeadingBytesCount_double(void)
{
	unsigned char b1[8] = {1,2,3,4,5,6,7,8}, b2[8] = {1,2,3,4,6,7,8,9}, b3[8] = {1,2,4,5,6,7,8,9};
	CU_ASSERT_EQUAL(compIdenticalLeadingBytesCount_double(b1, b2), 3);
	CU_ASSERT_EQUAL(compIdenticalLeadingBytesCount_double(b1, b3), 2);
	
}

void test_compIdenticalLeadingBytesCount_float(void)
{
	unsigned char b1[8] = {1,2,3,4,5,6,7,8}, b2[8] = {1,2,3,4,6,7,8,9}, b3[8] = {1,2,4,5,6,7,8,9};
	CU_ASSERT_EQUAL(compIdenticalLeadingBytesCount_float(b1, b2), 3);
	CU_ASSERT_EQUAL(compIdenticalLeadingBytesCount_float(b1, b3), 2);	
}

void test_addExactData(void)
{
	//TODO
}


/************* Test Runner Code goes here **************/

int main ( void )
{
   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if ( CUE_SUCCESS != CU_initialize_registry() )
      return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite( "test_conf_suite", init_suite, clean_suite );
   if ( NULL == pSuite ) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   if ( (NULL == CU_add_test(pSuite, "test_computeRangeSize_double", test_computeRangeSize_double)) ||
        (NULL == CU_add_test(pSuite, "test_computeRangeSize_float", test_computeRangeSize_float)) ||
        (NULL == CU_add_test(pSuite, "test_min_d", test_min_d)) ||
		(NULL == CU_add_test(pSuite, "test_max_d", test_max_d)) ||
        (NULL == CU_add_test(pSuite, "test_min_f", test_min_f)) ||
        (NULL == CU_add_test(pSuite, "test_max_d", test_max_d)) ||
        (NULL == CU_add_test(pSuite, "test_getRealPrecision_double", test_getRealPrecision_double)) ||
        (NULL == CU_add_test(pSuite, "test_getRealPrecision_float", test_getRealPrecision_float)) ||
        (NULL == CU_add_test(pSuite, "test_symTransform_8bytes", test_symTransform_8bytes)) ||
        (NULL == CU_add_test(pSuite, "test_symTransform_2bytes", test_symTransform_2bytes)) ||
        (NULL == CU_add_test(pSuite, "test_symTransform_4bytes", test_symTransform_4bytes)) ||
        (NULL == CU_add_test(pSuite, "test_compIdenticalLeadingBytesCount_double", test_compIdenticalLeadingBytesCount_double)) ||        
        (NULL == CU_add_test(pSuite, "test_compIdenticalLeadingBytesCount_float", test_compIdenticalLeadingBytesCount_float))
      )
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

   // Run all tests using the basic interface
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   printf("\n");
   CU_basic_show_failures(CU_get_failure_list());
	 unsigned int num_failures = CU_get_number_of_failures();
   printf("\n\n");
/*
   // Run all tests using the automated interface
   CU_automated_run_tests();
   CU_list_tests_to_file();

   // Run all tests using the console interface
   CU_console_run_tests();
*/
   /* Clean up registry and return */
   CU_cleanup_registry();
   return num_failures || CU_get_error();
}
