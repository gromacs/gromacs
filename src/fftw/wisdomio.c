/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <fftw-int.h>

/**************** import/export using file ***************/

static void file_emitter(char c, void *data)
{
     putc(c, (FILE *) data);
}

void fftw_export_wisdom_to_file(FILE *output_file)
{
     if (output_file)
	  fftw_export_wisdom(file_emitter, (void *) output_file);
}

static int file_get_input(void *data)
{
     return getc((FILE *) data);
}

fftw_status fftw_import_wisdom_from_file(FILE *input_file)
{
     if (!input_file)
	  return FFTW_FAILURE;
     return fftw_import_wisdom(file_get_input, (void *) input_file);
}

/*************** import/export using string **************/

static void emission_counter(char c, void *data)
{
     int *counter = (int *) data;

     ++*counter;
}

static void string_emitter(char c, void *data)
{
     char **output_string = (char **) data;

     *((*output_string)++) = c;
     **output_string = 0;
}

char *fftw_export_wisdom_to_string(void)
{
     int string_length = 0;
     char *s, *s2;

     fftw_export_wisdom(emission_counter, (void *) &string_length);

     s = (char *) fftw_malloc(sizeof(char) * (string_length + 1));
     if (!s)
	  return 0;
     s2 = s;

     fftw_export_wisdom(string_emitter, (void *) &s2);

     if (s + string_length != s2)
	  fftw_die("Unexpected output string length!");

     return s;
}

static int string_get_input(void *data)
{
     char **input_string = (char **) data;

     if (**input_string)
	  return *((*input_string)++);
     else
	  return 0;
}

fftw_status fftw_import_wisdom_from_string(const char *input_string)
{
     const char *s = input_string;

     if (!input_string)
	  return FFTW_FAILURE;
     return fftw_import_wisdom(string_get_input, (void *) &s);
}
