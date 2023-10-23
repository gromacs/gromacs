#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "sz.h"

#include <stdio.h>  // for printf

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_new_DBA(void)
{
	DynamicByteArray *dia = NULL;
	new_DBA(&dia, 1024);
	CU_ASSERT(dia->size==0 && dia->capacity==1024);
	free_DBA(dia);
}

void test_addDBA_Data(void)
{
	DynamicByteArray *dia = NULL;
	new_DBA(&dia, 1024);	
	addDBA_Data(dia, 1);
	addDBA_Data(dia, 2);
	addDBA_Data(dia, 3);
	addDBA_Data(dia, 4);
	addDBA_Data(dia, 5);
	CU_ASSERT_EQUAL(dia->size, 5);
	free_DBA(dia);
}

void test_convertDBAtoBytes(void)
{
	DynamicByteArray *dia = NULL;
	new_DBA(&dia, 1024);	
	addDBA_Data(dia, 1);
	addDBA_Data(dia, 2);
	addDBA_Data(dia, 3);
	addDBA_Data(dia, 4);
	addDBA_Data(dia, 5);
	unsigned char *data;	
	convertDBAtoBytes(dia, &data);
	
	unsigned char expected[5] = {1,2,3,4,5};
	CU_ASSERT_EQUAL_ARRAY_BYTE(data, expected, 5);
	
	free(data);
	free_DBA(dia);
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
   if ( (NULL == CU_add_test(pSuite, "test_new_DBA", test_new_DBA)) ||
        (NULL == CU_add_test(pSuite, "test_addDBA_Data", test_addDBA_Data)) ||
        (NULL == CU_add_test(pSuite, "test_convertDBAtoBytes", test_convertDBAtoBytes))
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
