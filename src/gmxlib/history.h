/*
 * history.h
 *
 *  Created on: Nov 27, 2009
 *      Author: rschulz
 */

#ifndef HISTORY_H_
#define HISTORY_H_
void init_history(int argc, char** argv);
void histopenfile(FILE* file, const char* fn, const char* mode);
int histclosefile(FILE** file);
void histaddinput(char* str);
void print_history();

#endif /* HISTORY_H_ */
