#include "parse.h"
void getFromString (const char* buf, int* n){
	*n=atoi(buf);
	printf("%d\n",*n);
}

void getFromString (const char* buf, size_t* n){
	*n=(size_t)(atoi(buf));
	printf("%lu\n",*n);
}

void getFromString (const char* buf, double* n){
	*n=(double)(atof(buf));
	printf("%15.6e\n",*n);
}

void getFromString (const char* buf, char* n){
	sscanf(buf,"%s",n);
	printf("%s\n",n);
}
