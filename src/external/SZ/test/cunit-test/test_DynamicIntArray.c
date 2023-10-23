#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"

#include <stdio.h>  // for printf

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_new_DIA(void)
{
	DynamicIntArray *dia = NULL;
	new_DIA(&dia, 1024);
	CU_ASSERT(dia->size==0 && dia->capacity==1024);
	free_DIA(dia);
}

void test_addDIA_Data(void)
{
	DynamicIntArray *dia = NULL;
	new_DIA(&dia, 1024);	
	addDIA_Data(dia, 1);
	addDIA_Data(dia, 2);
	addDIA_Data(dia, 3);
	addDIA_Data(dia, 4);
	addDIA_Data(dia, 5);
	CU_ASSERT_EQUAL(dia->size, 5);
	free_DIA(dia);
}

void test_convertDIAtoInts(void)
{
	DynamicIntArray *dia = NULL;
	new_DIA(&dia, 1024);	
	addDIA_Data(dia, 1);
	addDIA_Data(dia, 2);
	addDIA_Data(dia, 3);
	addDIA_Data(dia, 4);
	addDIA_Data(dia, 5);
	unsigned char *data;	
	convertDIAtoInts(dia, &data);
	
	unsigned char expected[5] = {1,2,3,4,5};
	CU_ASSERT_EQUAL_ARRAY_BYTE(data, expected, 5);
	
	free(data);
	free_DIA(dia);
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
   if ( (NULL == CU_add_test(pSuite, "test_new_DIA", test_new_DIA)) ||
        (NULL == CU_add_test(pSuite, "test_addDIA_Data", test_addDIA_Data)) ||
        (NULL == CU_add_test(pSuite, "test_convertDIAtoInts", test_convertDIAtoInts))
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
