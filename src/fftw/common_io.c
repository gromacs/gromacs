/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 */

#include <stdio.h>
#include <stdlib.h>

#include <fftw.h>

/**************** import/export using file ***************/

static void file_emitter(char c, void *data)
{
     putc(c,(FILE *) data);
}

void fftw_export_wisdom_to_file(FILE *output_file)
{
     if (output_file)
	  fftw_export_wisdom(file_emitter,(void *) output_file);
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

     s = fftw_malloc(sizeof(char) * (string_length + 1));
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
