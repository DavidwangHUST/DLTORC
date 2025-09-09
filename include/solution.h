#ifndef PRINT_DEF
#define PRINT_DEF
#include <string.h>	//for strings
#include <stdio.h>	//for printf,scanf
#include <stdlib.h>	//for atoi, atof
#endif

#ifndef SUNDIALS_DEF
#define SUNDIALS_DEF
#include <sundials/sundials_types.h>
#include <nvector/nvector_openmp.h>
#endif

int allocateSolution(size_t neq, int nThreads, N_Vector *y, N_Vector *ydot, N_Vector *atolv, N_Vector *constraints);
void freeSolution(N_Vector *y, N_Vector *ydot, N_Vector *atolv, N_Vector *constraints);
