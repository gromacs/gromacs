#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"
#include "rw.h"

#include <stdio.h>  // for printf
#include <stdlib.h>

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_checkFileSize(void)
{
	int state;
	int size = checkFileSize("../example/testdata/x86/testfloat_8_8_128.dat", &state);
	CU_ASSERT(SZ_SCES==state);
	CU_ASSERT_EQUAL(size, 32768);
}

void test_readByteData(void)
{
	int length, status;
	int x = 1;
	char *y = (char*)&x;
	dataEndianType = LITTLE_ENDIAN_SYSTEM;
    if(*y==1) 
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else
		sysEndianType = BIG_ENDIAN_SYSTEM;
	unsigned char* data = readByteData("../example/testdata/x86/testfloat_8_8_128.dat", &length, &status);
	unsigned char expected[10] = {129, 44, 112, 62, 37, 38, 112, 62, 196, 34};
	CU_ASSERT_EQUAL_ARRAY_BYTE(data, expected, 10);
	free(data);
}

void test_readDoubleData(void)
{
	int length, status;
	int x = 1;
	char *y = (char*)&x;
	dataEndianType = LITTLE_ENDIAN_SYSTEM;
    if(*y==1) 
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else
		sysEndianType = BIG_ENDIAN_SYSTEM;	
	
	double* data = readDoubleData("../example/testdata/x86/testdouble_8_8_128.dat", &length, &status);
	double expected[10] = {	0.225611633422454, 0.225634615576362, 0.225690839068313, 0.225738829973467, 0.225738288637820, 
							0.225691310405107, 0.225623293392676, 0.225563768873762, 0.225611633422454, 0.225634615576362};
	CU_ASSERT_EQUAL_ARRAY_DOUBLE(data, expected, 10, 1E-14);
	free(data);
}

void test_readFloatData(void)
{
	int length, status;
	int x = 1;
	char *y = (char*)&x;
	dataEndianType = LITTLE_ENDIAN_SYSTEM;
    if(*y==1) 
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else
		sysEndianType = BIG_ENDIAN_SYSTEM;	
	
	float* data = readFloatData("../example/testdata/x86/testfloat_8_8_128.dat", &length, &status);
	float expected[10] = {	0.23454477, 0.23452051, 0.23450762, 0.23450902, 0.23451449, 
							0.23453577, 0.23457345, 0.23459189, 0.23454477, 0.23452051};
	CU_ASSERT_EQUAL_ARRAY_FLOAT(data, expected, 10, 1E-8);
	free(data);
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
   if ( (NULL == CU_add_test(pSuite, "test_checkFileSize", test_checkFileSize)) ||
        (NULL == CU_add_test(pSuite, "test_readByteData", test_readByteData)) ||
        (NULL == CU_add_test(pSuite, "test_readFloatData", test_readFloatData)) ||
        (NULL == CU_add_test(pSuite, "test_readDoubleData", test_readDoubleData))
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
   return CU_get_error();
}
