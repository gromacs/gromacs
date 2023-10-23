#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "CUnit_Array.h"

#include "DynamicDoubleArray.h"

#include <stdio.h>  // for printf
#include <stdlib.h>

/* Test Suite setup and cleanup functions: */

int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void test_new_DDA(void)
{
	DynamicDoubleArray *dia = NULL;
	new_DDA(&dia, 1024);
	CU_ASSERT(dia->size==0 && dia->capacity==1024);
	free_DDA(dia);
}

void test_addDDA_Data(void)
{
	DynamicDoubleArray *dia = NULL;
	new_DDA(&dia, 1024);	
	addDDA_Data(dia, 1);
	addDDA_Data(dia, 2);
	addDDA_Data(dia, 3);
	addDDA_Data(dia, 4);
	addDDA_Data(dia, 5);
	CU_ASSERT_EQUAL(dia->size, 5);
	free_DDA(dia);
}

void test_convertDDAtoDoubles(void)
{
	DynamicDoubleArray *dia = NULL;
	new_DDA(&dia, 1024);	
	addDDA_Data(dia, 1);
	addDDA_Data(dia, 2);
	addDDA_Data(dia, 3);
	addDDA_Data(dia, 4);
	addDDA_Data(dia, 5);
	double *data;	
	convertDDAtoDoubles(dia, &data);
	
	double expected[5] = {1,2,3,4,5};
	CU_ASSERT_EQUAL_ARRAY_DOUBLE(data, expected, 5, 1E-6);
	
	free(data);
	free_DDA(dia);
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
   if ( (NULL == CU_add_test(pSuite, "test_new_DDA", test_new_DDA)) ||
        (NULL == CU_add_test(pSuite, "test_addDDA_Data", test_addDDA_Data)) ||
        (NULL == CU_add_test(pSuite, "test_convertDDAtoDoubles", test_convertDDAtoDoubles))
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
