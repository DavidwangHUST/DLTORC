#include "solution.h"
int allocateSolution(size_t neq, 
		     int nThreads,
		     N_Vector *y, 
		     N_Vector *ydot, 
		     N_Vector *atolv,
		     N_Vector *constraints){


  	*y = N_VNew_OpenMP(neq,nThreads);
  	if(*y==NULL){
		printf("Allocation Failed!\n");
		return(1);
	}
  	N_VConst((realtype)(0.0e0), *y);

  	*ydot = N_VNew_OpenMP(neq,nThreads);
  	if(*ydot==NULL){
		printf("Allocation Failed!\n");
		return(1);
	}
  	N_VConst((realtype)(0.0e0), *ydot);


  	*atolv = N_VNew_OpenMP(neq,nThreads);
  	if(*atolv==NULL){
		printf("Allocation Failed!\n");
		return(1);
	}
  	N_VConst((realtype)(0.0e0), *atolv);

  	*constraints = N_VNew_OpenMP(neq,nThreads);
  	if(*constraints==NULL){
		printf("Allocation Failed!\n");
		return(1);
	}
  	N_VConst((realtype)(0.0e0), *constraints);

	return(0);
}

void freeSolution(N_Vector *y, 
		  N_Vector *ydot, 
		  N_Vector *atolv,
		  N_Vector *constraints){

  	if(*y!=NULL){
		N_VDestroy_OpenMP(*y);
		printf("y Destroyed!\n");
	}
	if(*ydot!=NULL){
  		N_VDestroy_OpenMP(*ydot);
		printf("ydot Destroyed!\n");
	}
	if(*atolv!=NULL){
  		N_VDestroy_OpenMP(*atolv);
		printf("atolv Destroyed!\n");
	}

	if(*constraints!=NULL){
  		N_VDestroy_OpenMP(*constraints);
		printf("constraints Destroyed!\n");
	}
	printf("\n\n");
}
