#include <stdio.h>
#include <stdlib.h>
#include "NetCDFReader.h"

int netcdfReader(void *data, char *filename, char *dataset, int dataType)
{
	int ncid, varid, retval;

	/* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.*/
	if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
	{
		printf("Error: %s file cannot be open!\n", filename);
		exit(0);
	}

	/* Get the varid of the data variable, based on its name. */
	if ((retval = nc_inq_varid(ncid, dataset, &varid)))
	{
		printf("Error: %s dataset cannot be open!\n", dataset);
		exit(0);
	}

	/* Read the data. */
	if (dataType == 0)
	{
		if ((retval = nc_get_var_float(ncid, varid, (float*)data)))
		{
			printf("Error: %s dataset cannot be read!\n", dataset);
			exit(0);
		}
	}
	else
	{
		if ((retval = nc_get_var_double(ncid, varid, (double*)data)))
		{
			printf("Error: %s dataset cannot be read!\n", dataset);
			exit(0);
		}
	}

	/* Close the file, freeing all resources. */
	retval = nc_close(ncid);

	return 0;
}
