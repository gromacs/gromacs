#include <stdio.h>
#include <stdlib.h>
#include "adiosReader.h"

int adiosReader_1D (char *filename, size_t r1, int8_t **I8, int16_t **I16, int32_t **I32, int64_t **I64, uint8_t **U8, uint16_t **U16, uint32_t **U32, uint64_t **U64, float **R32, double **R64)
{
	adios_read_init_method(ADIOS_READ_METHOD_BP, 0, "verbose=3");

	// Open the file for reading
	ADIOS_FILE *f = adios_read_open_file(filename, ADIOS_READ_METHOD_BP, 0);

	// Check the variables exist
	ADIOS_VARINFO *var_i8 = adios_inq_var(f, "i8");
	ADIOS_VARINFO *var_i16 = adios_inq_var(f, "i16");
	ADIOS_VARINFO *var_i32 = adios_inq_var(f, "i32");
	ADIOS_VARINFO *var_i64 = adios_inq_var(f, "i64");
	ADIOS_VARINFO *var_u8 = adios_inq_var(f, "u8");
	ADIOS_VARINFO *var_u16 = adios_inq_var(f, "u16");
	ADIOS_VARINFO *var_u32 = adios_inq_var(f, "u32");
	ADIOS_VARINFO *var_u64 = adios_inq_var(f, "u64");
	ADIOS_VARINFO *var_r32 = adios_inq_var(f, "r32");
	ADIOS_VARINFO *var_r64 = adios_inq_var(f, "r64");

	if (var_i8  != NULL) *I8 = (int8_t *)malloc(r1*sizeof(int8_t));
	if (var_i16 != NULL) *I16 = (int16_t *)malloc(r1*sizeof(int16_t));
	if (var_i32 != NULL) *I32 = (int32_t *)malloc(r1*sizeof(int32_t));
	if (var_i64 != NULL) *I64 = (int64_t *)malloc(r1*sizeof(int64_t));
	if (var_u8  != NULL) *U8 = (uint8_t *)malloc(r1*sizeof(uint8_t));
	if (var_u16 != NULL) *U16 = (uint16_t *)malloc(r1*sizeof(uint16_t));
	if (var_u32 != NULL) *U32 = (uint32_t *)malloc(r1*sizeof(uint32_t));
	if (var_u64 != NULL) *U64 = (uint64_t *)malloc(r1*sizeof(uint64_t));
	if (var_r32 != NULL) *R32 = (float *)malloc(r1*sizeof(float));
	if (var_r64 != NULL) *R64 = (double *)malloc(r1*sizeof(double));

	uint64_t start[1] = {0};
	uint64_t count[1] = {r1};
	ADIOS_SELECTION *sel = adios_selection_boundingbox(1, start, count);

	// Read stuff
	//	for (size_t t = 0; t < NSteps; ++t)
	//	{
	size_t t = 1; // Only considr nstep = 1
	
	// Read the current step
	adios_schedule_read_byid(f, sel, var_i8->varid, t, 1, *I8);
	adios_schedule_read_byid(f, sel, var_i16->varid, t, 1, *I16);
	adios_schedule_read_byid(f, sel, var_i32->varid, t, 1, *I32);
	adios_schedule_read_byid(f, sel, var_i64->varid, t, 1, *I64);
	adios_schedule_read_byid(f, sel, var_u8->varid, t, 1, *U8);
	adios_schedule_read_byid(f, sel, var_u16->varid, t, 1, *U16);
	adios_schedule_read_byid(f, sel, var_u32->varid, t, 1, *U32);
	adios_schedule_read_byid(f, sel, var_u64->varid, t, 1, *U64);
	adios_schedule_read_byid(f, sel, var_r32->varid, t, 1, *R32);
	adios_schedule_read_byid(f, sel, var_r64->varid, t, 1, *R64);
	adios_perform_reads(f, 1);
	//	}

	adios_selection_delete(sel);

	// Cleanup variable structures
	adios_free_varinfo(var_i8);
	adios_free_varinfo(var_i16);
	adios_free_varinfo(var_i32);
	adios_free_varinfo(var_i64);
	adios_free_varinfo(var_u8);
	adios_free_varinfo(var_u16);
	adios_free_varinfo(var_u32);
	adios_free_varinfo(var_u64);
	adios_free_varinfo(var_r32);
	adios_free_varinfo(var_r64);

	// Cleanup file
	adios_read_close(f);

	adios_read_finalize_method(ADIOS_READ_METHOD_BP);

	return 0;
}

int adiosReader_2D (char *filename, size_t r1, size_t r2, int8_t **I8, int16_t **I16, int32_t **I32, int64_t **I64, uint8_t **U8, uint16_t **U16, uint32_t **U32, uint64_t **U64, float **R32, double **R64)
{
	adios_read_init_method(ADIOS_READ_METHOD_BP, 0, "verbose=3");
	// Open the file for reading
	ADIOS_FILE *f = adios_read_open_file(filename, ADIOS_READ_METHOD_BP, 0);
	// Check the variables exist
	ADIOS_VARINFO *var_i8 = adios_inq_var(f, "i8");
	ADIOS_VARINFO *var_i16 = adios_inq_var(f, "i16");
	ADIOS_VARINFO *var_i32 = adios_inq_var(f, "i32");
	ADIOS_VARINFO *var_i64 = adios_inq_var(f, "i64");
	ADIOS_VARINFO *var_u8 = adios_inq_var(f, "u8");
	ADIOS_VARINFO *var_u16 = adios_inq_var(f, "u16");
	ADIOS_VARINFO *var_u32 = adios_inq_var(f, "u32");
	ADIOS_VARINFO *var_u64 = adios_inq_var(f, "u64");
	ADIOS_VARINFO *var_r32 = adios_inq_var(f, "r32");
	ADIOS_VARINFO *var_r64 = adios_inq_var(f, "r64");

	// If the size of the array is smaller than the data
	// the result is weird... double and uint64_t would get completely
	// garbage data
	
	if (var_i8  != NULL) *I8 = (int8_t *)malloc(r1*r2*sizeof(int8_t));
	if (var_i16 != NULL) *I16 = (int16_t *)malloc(r1*r2*sizeof(int16_t));
	if (var_i32 != NULL) *I32 = (int32_t *)malloc(r1*r2*sizeof(int32_t));
	if (var_i64 != NULL) *I64 = (int64_t *)malloc(r1*r2*sizeof(int64_t));
	if (var_u8  != NULL) *U8 = (uint8_t *)malloc(r1*r2*sizeof(uint8_t));
	if (var_u16 != NULL) *U16 = (uint16_t *)malloc(r1*r2*sizeof(uint16_t));
	if (var_u32 != NULL) *U32 = (uint32_t *)malloc(r1*r2*sizeof(uint32_t));
	if (var_u64 != NULL) *U64 = (uint64_t *)malloc(r1*r2*sizeof(uint64_t));
	if (var_r32 != NULL) *R32 = (float *)malloc(r1*r2*sizeof(float));
	if (var_r64 != NULL) *R64 = (double *)malloc(r1*r2*sizeof(double));

	uint64_t start[2] = {0, 0};
	uint64_t count[2] = {r2, r1};
	ADIOS_SELECTION *sel = adios_selection_boundingbox(2, start, count);

	// Read stuff
	//	for (size_t t = 0; t < NSteps; ++t)
	//	{
	size_t t = 1; // Only considr nstep = 1
	
	// Read the current step
	adios_schedule_read_byid(f, sel, var_i8->varid, t, 1, *I8);
	adios_schedule_read_byid(f, sel, var_i16->varid, t, 1, *I16);
	adios_schedule_read_byid(f, sel, var_i32->varid, t, 1, *I32);
	adios_schedule_read_byid(f, sel, var_i64->varid, t, 1, *I64);
	adios_schedule_read_byid(f, sel, var_u8->varid, t, 1, *U8);
	adios_schedule_read_byid(f, sel, var_u16->varid, t, 1, *U16);
	adios_schedule_read_byid(f, sel, var_u32->varid, t, 1, *U32);
	adios_schedule_read_byid(f, sel, var_u64->varid, t, 1, *U64);
	adios_schedule_read_byid(f, sel, var_r32->varid, t, 1, *R32);
	adios_schedule_read_byid(f, sel, var_r64->varid, t, 1, *R64);
	adios_perform_reads(f, 1);
	//	}

	adios_selection_delete(sel);

	// Cleanup variable structures
	adios_free_varinfo(var_i8);
	adios_free_varinfo(var_i16);
	adios_free_varinfo(var_i32);
	adios_free_varinfo(var_i64);
	adios_free_varinfo(var_u8);
	adios_free_varinfo(var_u16);
	adios_free_varinfo(var_u32);
	adios_free_varinfo(var_u64);
	adios_free_varinfo(var_r32);
	adios_free_varinfo(var_r64);

	// Cleanup file
	adios_read_close(f);

	adios_read_finalize_method(ADIOS_READ_METHOD_BP);

	return 0;
}

int adiosReader_3D (char *filename, size_t r1, size_t r2, size_t r3, int8_t **I8, int16_t **I16, int32_t **I32, int64_t **I64, uint8_t **U8, uint16_t **U16, uint32_t **U32, uint64_t **U64, float **R32, double **R64)
{
	adios_read_init_method(ADIOS_READ_METHOD_BP, 0, "verbose=3");
	// Open the file for reading
	ADIOS_FILE *f = adios_read_open_file(filename, ADIOS_READ_METHOD_BP, 0);
	// Check the variables exist
	ADIOS_VARINFO *var_i8 = adios_inq_var(f, "i8");
	ADIOS_VARINFO *var_i16 = adios_inq_var(f, "i16");
	ADIOS_VARINFO *var_i32 = adios_inq_var(f, "i32");
	ADIOS_VARINFO *var_i64 = adios_inq_var(f, "i64");
	ADIOS_VARINFO *var_u8 = adios_inq_var(f, "u8");
	ADIOS_VARINFO *var_u16 = adios_inq_var(f, "u16");
	ADIOS_VARINFO *var_u32 = adios_inq_var(f, "u32");
	ADIOS_VARINFO *var_u64 = adios_inq_var(f, "u64");
	ADIOS_VARINFO *var_r32 = adios_inq_var(f, "r32");
	ADIOS_VARINFO *var_r64 = adios_inq_var(f, "r64");

	// If the size of the array is smaller than the data
	// the result is weird... double and uint64_t would get completely
	// garbage data
	
	if (var_i8  != NULL) *I8 = (int8_t *)malloc(r1*r2*r3*sizeof(int8_t));
	if (var_i16 != NULL) *I16 = (int16_t *)malloc(r1*r2*r3*sizeof(int16_t));
	if (var_i32 != NULL) *I32 = (int32_t *)malloc(r1*r2*r3*sizeof(int32_t));
	if (var_i64 != NULL) *I64 = (int64_t *)malloc(r1*r2*r3*sizeof(int64_t));
	if (var_u8  != NULL) *U8 = (uint8_t *)malloc(r1*r2*r3*sizeof(uint8_t));
	if (var_u16 != NULL) *U16 = (uint16_t *)malloc(r1*r2*r3*sizeof(uint16_t));
	if (var_u32 != NULL) *U32 = (uint32_t *)malloc(r1*r2*r3*sizeof(uint32_t));
	if (var_u64 != NULL) *U64 = (uint64_t *)malloc(r1*r2*r3*sizeof(uint64_t));
	if (var_r32 != NULL) *R32 = (float *)malloc(r1*r2*r3*sizeof(float));
	if (var_r64 != NULL) *R64 = (double *)malloc(r1*r2*r3*sizeof(double));

	uint64_t start[3] = {0, 0, 0};
	uint64_t count[3] = {r3, r2, r1};
	ADIOS_SELECTION *sel = adios_selection_boundingbox(3, start, count);

	// Read stuff
	//	for (size_t t = 0; t < NSteps; ++t)
	//	{
	size_t t = 1; // Only considr nstep = 1
	
	// Read the current step
	adios_schedule_read_byid(f, sel, var_i8->varid, t, 1, *I8);
	adios_schedule_read_byid(f, sel, var_i16->varid, t, 1, *I16);
	adios_schedule_read_byid(f, sel, var_i32->varid, t, 1, *I32);
	adios_schedule_read_byid(f, sel, var_i64->varid, t, 1, *I64);
	adios_schedule_read_byid(f, sel, var_u8->varid, t, 1, *U8);
	adios_schedule_read_byid(f, sel, var_u16->varid, t, 1, *U16);
	adios_schedule_read_byid(f, sel, var_u32->varid, t, 1, *U32);
	adios_schedule_read_byid(f, sel, var_u64->varid, t, 1, *U64);
	adios_schedule_read_byid(f, sel, var_r32->varid, t, 1, *R32);
	adios_schedule_read_byid(f, sel, var_r64->varid, t, 1, *R64);
	adios_perform_reads(f, 1);
	//	}

	adios_selection_delete(sel);

	// Cleanup variable structures
	adios_free_varinfo(var_i8);
	adios_free_varinfo(var_i16);
	adios_free_varinfo(var_i32);
	adios_free_varinfo(var_i64);
	adios_free_varinfo(var_u8);
	adios_free_varinfo(var_u16);
	adios_free_varinfo(var_u32);
	adios_free_varinfo(var_u64);
	adios_free_varinfo(var_r32);
	adios_free_varinfo(var_r64);

	// Cleanup file
	adios_read_close(f);

	adios_read_finalize_method(ADIOS_READ_METHOD_BP);

	return 0;
}


int adiosReader_4D (char *filename, size_t r1, size_t r2, size_t r3, size_t r4, int8_t **I8, int16_t **I16, int32_t **I32, int64_t **I64, uint8_t **U8, uint16_t **U16, uint32_t **U32, uint64_t **U64, float **R32, double **R64)
{
	adios_read_init_method(ADIOS_READ_METHOD_BP, 0, "verbose=3");
	// Open the file for reading
	ADIOS_FILE *f = adios_read_open_file(filename, ADIOS_READ_METHOD_BP, 0);
	// Check the variables exist
	ADIOS_VARINFO *var_i8 = adios_inq_var(f, "i8");
	ADIOS_VARINFO *var_i16 = adios_inq_var(f, "i16");
	ADIOS_VARINFO *var_i32 = adios_inq_var(f, "i32");
	ADIOS_VARINFO *var_i64 = adios_inq_var(f, "i64");
	ADIOS_VARINFO *var_u8 = adios_inq_var(f, "u8");
	ADIOS_VARINFO *var_u16 = adios_inq_var(f, "u16");
	ADIOS_VARINFO *var_u32 = adios_inq_var(f, "u32");
	ADIOS_VARINFO *var_u64 = adios_inq_var(f, "u64");
	ADIOS_VARINFO *var_r32 = adios_inq_var(f, "r32");
	ADIOS_VARINFO *var_r64 = adios_inq_var(f, "r64");

	// If the size of the array is smaller than the data
	// the result is weird... double and uint64_t would get completely
	// garbage data
	
	if (var_i8  != NULL) *I8 = (int8_t *)malloc(r1*r2*r3*r4*sizeof(int8_t));
	if (var_i16 != NULL) *I16 = (int16_t *)malloc(r1*r2*r3*r4*sizeof(int16_t));
	if (var_i32 != NULL) *I32 = (int32_t *)malloc(r1*r2*r3*r4*sizeof(int32_t));
	if (var_i64 != NULL) *I64 = (int64_t *)malloc(r1*r2*r3*r4*sizeof(int64_t));
	if (var_u8  != NULL) *U8 = (uint8_t *)malloc(r1*r2*r3*r4*sizeof(uint8_t));
	if (var_u16 != NULL) *U16 = (uint16_t *)malloc(r1*r2*r3*r4*sizeof(uint16_t));
	if (var_u32 != NULL) *U32 = (uint32_t *)malloc(r1*r2*r3*r4*sizeof(uint32_t));
	if (var_u64 != NULL) *U64 = (uint64_t *)malloc(r1*r2*r3*r4*sizeof(uint64_t));
	if (var_r32 != NULL) *R32 = (float *)malloc(r1*r2*r3*r4*sizeof(float));
	if (var_r64 != NULL) *R64 = (double *)malloc(r1*r2*r3*r4*sizeof(double));

	uint64_t start[4] = {0, 0, 0, 0};
	uint64_t count[4] = {r4, r3, r2, r1};
	ADIOS_SELECTION *sel = adios_selection_boundingbox(4, start, count);

	// Read stuff
	//	for (size_t t = 0; t < NSteps; ++t)
	//	{
	size_t t = 1; // Only considr nstep = 1
	
	// Read the current step
	adios_schedule_read_byid(f, sel, var_i8->varid, t, 1, *I8);
	adios_schedule_read_byid(f, sel, var_i16->varid, t, 1, *I16);
	adios_schedule_read_byid(f, sel, var_i32->varid, t, 1, *I32);
	adios_schedule_read_byid(f, sel, var_i64->varid, t, 1, *I64);
	adios_schedule_read_byid(f, sel, var_u8->varid, t, 1, *U8);
	adios_schedule_read_byid(f, sel, var_u16->varid, t, 1, *U16);
	adios_schedule_read_byid(f, sel, var_u32->varid, t, 1, *U32);
	adios_schedule_read_byid(f, sel, var_u64->varid, t, 1, *U64);
	adios_schedule_read_byid(f, sel, var_r32->varid, t, 1, *R32);
	adios_schedule_read_byid(f, sel, var_r64->varid, t, 1, *R64);
	adios_perform_reads(f, 1);
	//	}

	adios_selection_delete(sel);

	// Cleanup variable structures
	adios_free_varinfo(var_i8);
	adios_free_varinfo(var_i16);
	adios_free_varinfo(var_i32);
	adios_free_varinfo(var_i64);
	adios_free_varinfo(var_u8);
	adios_free_varinfo(var_u16);
	adios_free_varinfo(var_u32);
	adios_free_varinfo(var_u64);
	adios_free_varinfo(var_r32);
	adios_free_varinfo(var_r64);

	// Cleanup file
	adios_read_close(f);

	adios_read_finalize_method(ADIOS_READ_METHOD_BP);

	return 0;
}


int adiosReader_5D (char *filename, size_t r1, size_t r2, size_t r3, size_t r4, size_t r5, int8_t **I8, int16_t **I16, int32_t **I32, int64_t **I64, uint8_t **U8, uint16_t **U16, uint32_t **U32, uint64_t **U64, float **R32, double **R64)
{
	adios_read_init_method(ADIOS_READ_METHOD_BP, 0, "verbose=3");
	// Open the file for reading
	ADIOS_FILE *f = adios_read_open_file(filename, ADIOS_READ_METHOD_BP, 0);
	// Check the variables exist
	ADIOS_VARINFO *var_i8 = adios_inq_var(f, "i8");
	ADIOS_VARINFO *var_i16 = adios_inq_var(f, "i16");
	ADIOS_VARINFO *var_i32 = adios_inq_var(f, "i32");
	ADIOS_VARINFO *var_i64 = adios_inq_var(f, "i64");
	ADIOS_VARINFO *var_u8 = adios_inq_var(f, "u8");
	ADIOS_VARINFO *var_u16 = adios_inq_var(f, "u16");
	ADIOS_VARINFO *var_u32 = adios_inq_var(f, "u32");
	ADIOS_VARINFO *var_u64 = adios_inq_var(f, "u64");
	ADIOS_VARINFO *var_r32 = adios_inq_var(f, "r32");
	ADIOS_VARINFO *var_r64 = adios_inq_var(f, "r64");

	// If the size of the array is smaller than the data
	// the result is weird... double and uint64_t would get completely
	// garbage data
	
	if (var_i8  != NULL) *I8 = (int8_t *)malloc(r1*r2*r3*r4*r5*sizeof(int8_t));
	if (var_i16 != NULL) *I16 = (int16_t *)malloc(r1*r2*r3*r4*r5*sizeof(int16_t));
	if (var_i32 != NULL) *I32 = (int32_t *)malloc(r1*r2*r3*r4*r5*sizeof(int32_t));
	if (var_i64 != NULL) *I64 = (int64_t *)malloc(r1*r2*r3*r4*r5*sizeof(int64_t));
	if (var_u8  != NULL) *U8 = (uint8_t *)malloc(r1*r2*r3*r4*r5*sizeof(uint8_t));
	if (var_u16 != NULL) *U16 = (uint16_t *)malloc(r1*r2*r3*r4*r5*sizeof(uint16_t));
	if (var_u32 != NULL) *U32 = (uint32_t *)malloc(r1*r2*r3*r4*r5*sizeof(uint32_t));
	if (var_u64 != NULL) *U64 = (uint64_t *)malloc(r1*r2*r3*r4*r5*sizeof(uint64_t));
	if (var_r32 != NULL) *R32 = (float *)malloc(r1*r2*r3*r4*r5*sizeof(float));
	if (var_r64 != NULL) *R64 = (double *)malloc(r1*r2*r3*r4*r5*sizeof(double));

	uint64_t start[5] = {0, 0, 0, 0, 0};
	uint64_t count[5] = {r5, r4, r3, r2, r1};
	ADIOS_SELECTION *sel = adios_selection_boundingbox(5, start, count);

	// Read stuff
	//	for (size_t t = 0; t < NSteps; ++t)
	//	{
	size_t t = 1; // Only considr nstep = 1
	
	// Read the current step
	adios_schedule_read_byid(f, sel, var_i8->varid, t, 1, *I8);
	adios_schedule_read_byid(f, sel, var_i16->varid, t, 1, *I16);
	adios_schedule_read_byid(f, sel, var_i32->varid, t, 1, *I32);
	adios_schedule_read_byid(f, sel, var_i64->varid, t, 1, *I64);
	adios_schedule_read_byid(f, sel, var_u8->varid, t, 1, *U8);
	adios_schedule_read_byid(f, sel, var_u16->varid, t, 1, *U16);
	adios_schedule_read_byid(f, sel, var_u32->varid, t, 1, *U32);
	adios_schedule_read_byid(f, sel, var_u64->varid, t, 1, *U64);
	adios_schedule_read_byid(f, sel, var_r32->varid, t, 1, *R32);
	adios_schedule_read_byid(f, sel, var_r64->varid, t, 1, *R64);
	adios_perform_reads(f, 1);
	//	}

	adios_selection_delete(sel);

	// Cleanup variable structures
	adios_free_varinfo(var_i8);
	adios_free_varinfo(var_i16);
	adios_free_varinfo(var_i32);
	adios_free_varinfo(var_i64);
	adios_free_varinfo(var_u8);
	adios_free_varinfo(var_u16);
	adios_free_varinfo(var_u32);
	adios_free_varinfo(var_u64);
	adios_free_varinfo(var_r32);
	adios_free_varinfo(var_r64);

	// Cleanup file
	adios_read_close(f);

	adios_read_finalize_method(ADIOS_READ_METHOD_BP);

	return 0;
}
