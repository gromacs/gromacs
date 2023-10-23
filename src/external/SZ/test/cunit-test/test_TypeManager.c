#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"

#include <stdio.h>  // for printf
#include <stdlib.h>

/* Test Suite setup and cleanup functions: */

int
init_suite(void)
{
  exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
  return 0;
}

int
clean_suite(void)
{
  free(exe_params);
  return 0;
}

/************* Test case functions ****************/

void test_convertIntArray2ByteArray_fast_2b(void)
{
	unsigned char b_array[35] = {0,0,0,1,2, 3,2,1,1,0, 0,1,0,0,1, 2,2,2,1,1, 0,3,3,3,1, 2,3,2,1,0, 1,1,1,1,0};
	unsigned char *b_result;
	int b_result_len = convertIntArray2ByteArray_fast_2b(b_array, 35, &b_result);
	unsigned char b_expected[9] = {1,185,65,6,165,63,110,69,80};
	CU_ASSERT_EQUAL(b_result_len, 9);
	CU_ASSERT_EQUAL_ARRAY_BYTE(b_result, b_expected, 9);
	free(b_result);
}

void test_convertByteArray2IntArray_fast_2b(void)
{
	unsigned char bytes[9] = {1,185,65,6,165,63,110,69,80};
	unsigned char *intArray;
	convertByteArray2IntArray_fast_2b(35, bytes, 9, &intArray);
	unsigned char expected[35] = {0,0,0,1,2, 3,2,1,1,0, 0,1,0,0,1, 2,2,2,1,1, 0,3,3,3,1, 2,3,2,1,0, 1,1,1,1,0};
	CU_ASSERT_EQUAL_ARRAY_BYTE(intArray, expected, 35);
	free(intArray);
}

void test_convertIntArray2ByteArray_fast_3b(void)
{
	unsigned char b3_array[35] = {0,0,0,1,2, 6,7,7,7,7, 6,6,6,4,4, 4,2,1,2,1, 2,5,4,0,0, 0,5,6,4,3, 5,5,2,2,1};
	unsigned char *b3_result;
	int b3_result_len = convertIntArray2ByteArray_fast_3b(b3_array, 35, &b3_result);
	unsigned char b3_expected[14] = {0,21,191,255,109, 36,69,21,96,2, 232,237,72,128};
	CU_ASSERT_EQUAL(b3_result_len, 14);
	CU_ASSERT_EQUAL_ARRAY_BYTE(b3_result, b3_expected, 14);	
	free(b3_result);
}

void test_convertByteArray2IntArray_fast_3b()
{
	unsigned char bytes[14] = {0,21,191,255,109, 36,69,21,96,2, 232,237,72,128};
	unsigned char *intArray;
	convertByteArray2IntArray_fast_3b(35, bytes, 14, &intArray);
	unsigned char expected[35] = {0,0,0,1,2, 6,7,7,7,7, 6,6,6,4,4, 4,2,1,2,1, 2,5,4,0,0, 0,5,6,4,3, 5,5,2,2,1};
	CU_ASSERT_EQUAL_ARRAY_BYTE(intArray, expected, 35);
	free(intArray);	
}

void test_getLeftMovingSteps(void)
{
	CU_ASSERT_EQUAL(getLeftMovingSteps(1, 4), 3);
	CU_ASSERT_EQUAL(getLeftMovingSteps(2, 4), 2);
	CU_ASSERT_EQUAL(getLeftMovingSteps(3, 4), 1);
	CU_ASSERT_EQUAL(getLeftMovingSteps(4, 4), 0);	
}

void test_computeBitNumRequired(void)
{
	int i;
	int result[100];
	for(i=0;i<100;i++)
	{
		result[i] = computeBitNumRequired(i);
	}
	int expected[100];
	expected[0] = 0;
	expected[1] = 1;
	expected[2] = expected[3] = 2;
	for(i=4;i<8;i++)
		expected[i] = 3;
	for(i=8;i<16;i++)
		expected[i] = 4;
	for(i=16;i<32;i++)
		expected[i] = 5;
	for(i=32;i<64;i++)
		expected[i] = 6;
	for(i=64;i<100;i++)
		expected[i] = 7;
		
	CU_ASSERT_EQUAL_ARRAY_INT(result, expected, 100);
}

void test_convertIntArray2ByteArray_fast_dynamic(void)
{
	//uneasy to implement for unit testing, because it needs some particular input generated during the compression.
}

void test_decompressBitArraybySimpleLZ77(void)
{
	//decompressBitArraybySimpleLZ77 not being used by the main implementation of SZ (it's used only for the case with 'preserved' data)
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
   if ( (NULL == CU_add_test(pSuite, "test_convertIntArray2ByteArray_fast_2b", test_convertIntArray2ByteArray_fast_2b)) ||
        (NULL == CU_add_test(pSuite, "test_convertByteArray2IntArray_fast_2b", test_convertByteArray2IntArray_fast_2b)) ||
        (NULL == CU_add_test(pSuite, "test_convertIntArray2ByteArray_fast_3b", test_convertIntArray2ByteArray_fast_3b)) ||
		(NULL == CU_add_test(pSuite, "test_convertIntArray2ByteArray_fast_3b", test_convertIntArray2ByteArray_fast_3b)) ||
        (NULL == CU_add_test(pSuite, "test_getLeftMovingSteps", test_getLeftMovingSteps)) ||
        (NULL == CU_add_test(pSuite, "test_computeBitNumRequired", test_computeBitNumRequired))
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
