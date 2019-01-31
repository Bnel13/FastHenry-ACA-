#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "induct.h"
#include <stdio.h>

extern FILE* current_output_file;
extern char* current_output_fname;

void start_write_currents(void);
void dump_current_output(char* outp);
void end_write_currents(void);
void place_current_output(EXTERNAL* in_port, EXTERNAL* out_port, CX current);
void declare_current_output_file(char* outfile);

#endif /* OUTPUT_H_ */
