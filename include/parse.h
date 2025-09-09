#ifndef PRINT_DEF
#define PRINT_DEF
#include <string.h>	//for strings
#include <stdio.h>	//for printf,scanf
#include <stdlib.h>	//for atoi, atof
#endif

#ifndef PARSE_DEF
#define PARSE_DEF
#define MAXBUFLEN 200

void getFromString (const char* buf, int* n);
void getFromString (const char* buf, size_t* n);
void getFromString (const char* buf, double* n);
void getFromString (const char* buf, char* n);

int parseString(FILE* input, const char* keyword, const size_t bufLen, char* n);

template<typename T>
int parseNumber(FILE* input, const char* keyword, const size_t bufLen, T* n);

template<typename T>
int parseArray(FILE* input, const char* keyword, const size_t bufLen,
		const size_t arrLen, T y[]);


#include "parse.hpp"

#endif


