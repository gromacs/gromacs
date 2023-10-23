#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"

#include "sz.h"

#include <stdio.h>  // for printf

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_roundUpToPowerOf2(void)
{
	CU_ASSERT_EQUAL(roundUpToPowerOf2(10), 16);
	CU_ASSERT_EQUAL(roundUpToPowerOf2(16), 16);
	CU_ASSERT_EQUAL(roundUpToPowerOf2(17), 32);
	CU_ASSERT_EQUAL(roundUpToPowerOf2(3), 4);
}

void test_SZ_LoadConf(void)
{
	sz_cfgFile = "sz.config";
	CU_ASSERT_EQUAL(SZ_LoadConf(),SZ_SCES);
}

void test_checkVersion(void)
{
	char verNum[4] = {SZ_VER_MAJOR,SZ_VER_MINOR,SZ_VER_BUILD,SZ_VER_REVISION};
	CU_ASSERT_EQUAL(checkVersion(verNum), 1);
	verNum[0] = SZ_VER_MAJOR-1;
	CU_ASSERT_EQUAL(checkVersion(verNum), 0);
	verNum[0] = SZ_VER_MAJOR;
	verNum[1] = SZ_VER_MINOR+1;
	CU_ASSERT_EQUAL(checkVersion(verNum), 0);
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
   if ( (NULL == CU_add_test(pSuite, "test_roundUpToPowerOf2", test_roundUpToPowerOf2)) ||
        (NULL == CU_add_test(pSuite, "test_SZ_LoadConf", test_SZ_LoadConf)) ||
        (NULL == CU_add_test(pSuite, "test_checkVersion", test_checkVersion))
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
