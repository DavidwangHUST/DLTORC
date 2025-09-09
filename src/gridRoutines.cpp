#include "gridRoutines.h"
#include <stdio.h>
inline double l(const double* x, 
	 	const double* a,
	 	const double* w,
		const double* fac,
		const int* refineLeft){
	if(*refineLeft==0){
		return(tanh(-*a*(*x+*w*100.0e0)));
	}else{
		return(tanh(-*a*(*x-*w*(*fac))));
	}
}

inline double r(const double* x, 
		const double* a, 
		const double* w,
		const double* fac,
		const int* refineRight){
	if(*refineRight==0){
		return(tanh(*a*(*x-(1.0e0+*w*100.0e0))));
	}else{
		return(tanh(*a*(*x-(1.0e0-*w*(*fac)))));
	}
}

inline double f(const double* x,
		const double* a,
		const double* c,
		const double* w){
	return(tanh(-*a*(*x-(*c+*w)))
		+tanh(-*a*((*x-1.0e0)-(*c+*w)))
		+tanh(-*a*((*x+1.0e0)-(*c+*w))));
}

inline double g(const double* x, 
		const double* a, 
		const double* c, 
		const double* w){
	return(tanh(*a*(*x-(*c-*w)))
	       +tanh(*a*((*x-1.0e0)-(*c-*w)))
	       +tanh(*a*((*x+1.0e0)-(*c-*w))));
}

inline double rho(const double* x,
	   	  const double* a,
	   	  const double* c,
	   	  const double* w,
	   	  const double* mag,
	   	  const double* leftFac,
	   	  const double* rightFac,
		  const int* refineLeft,
		  const int* refineRight){

	return(((2.0e0+f(x,a,c,w)
		 +g(x,a,c,w)
		 +l(x,a,w,leftFac,refineLeft)
		 +r(x,a,w,rightFac,refineRight))*0.5e0)
		  *(*mag-1.0e0)+1.0e0);
}

size_t maxPoints(const size_t basePts,
	   	 const double* a,
	   	 const double* w,
	   	 const double* mag,
	   	 const double* leftFac,
	   	 const double* rightFac,
	   	 const int* refineLeft,
	   	 const int* refineRight){
	double dx=1.0e0/((double)(basePts)-1.0e0);
	double y=0.0e0;
	size_t i=0;
	double r=0.0e0;
	double t=0.5e0;
	while(y<=1.0e0){
		r=rho(&y,a,&t,w,mag,leftFac,rightFac,refineLeft,refineRight);
		y=y+(dx/r);
		i++;
	}
	return(i);
}

void fillGrid(const size_t* basePts,
	      const size_t* nPts,
	      const double* a,
	      const double* c,
	      const double* w,
	      const double* mag,
	      const double* leftFac,
	      const double* rightFac,
	      const int* refineLeft,
	      const int* refineRight,
	      double x[]){

	FILE* out;out=fopen("tmp.dat","w");

	double y=0.0e0;
	double r=0.0e0;
	double dx=1.0e0/((double)(*basePts)-1.0e0);
	for(size_t j=0;j<*nPts;j++){
		r=rho(&y,a,c,w,mag,leftFac,rightFac,refineLeft,refineRight);
		fprintf(out, "%15.15e\n",dx/r);
		y=y+(dx/r);
	}
	fclose(out);
	
	double dxp[*nPts-1];
	for (size_t j = 0; j < *nPts; j++) {
		dxp[j]=0.0e0;
	}

	FILE* tmp;tmp=fopen("tmp.dat","r");
	char buf[MAXBUFLEN];
	size_t i=0;
	while (fgets(buf,MAXBUFLEN, tmp)!=NULL){
		sscanf(buf, "%lf", &y);
		dxp[i]=y;
		i++;
	}
	fclose(tmp);

	double sum=0.0e0;
	double err=0.0e0;
	double fix=0.0e0;
	for(size_t j=0;j<*nPts-1;j++){
		sum+=dxp[j];
	}
	err=1.0e0-sum;
	printf("sum before correction: %15.6e\n",sum);
	printf("err before correction: %15.6e\n",err);
	fix=err/((double)(*nPts));
	sum=0.0e0;
	for(size_t j=0;j<*nPts-1;j++){
		dxp[j]+=fix;
		sum+=dxp[j];
	}
	err=1.0e0-sum;
	printf("sum after correction:%15.6e\n",sum);
	printf("err after correction: %15.6e\n",err);

	x[0]=0.0e0;
	for(size_t j=0;j<*nPts-1;j++){
		x[j+1]=x[j]+dxp[j];
	}
	x[*nPts-1]=1.0e0;

}

double safePosition(double c, double w){
	if(c<w){
		return(w);
	}
	else if(c>1.0e0-w){
		return(1.0e0-w);
	}
	else{
		return(c);
	}
}

int reGrid(UserGrid grid, double position){

	printf("before regrid: %ld\n", grid->nPts);
	double xx[grid->nPts];

	fillGrid(&grid->basePts,
		 &grid->nPts,
		 &grid->a,
		 &position,
		 &grid->w,
		 &grid->mag,
		 &grid->leftFac,
		 &grid->rightFac,
		 &grid->refineLeft,
		 &grid->refineRight,
		 xx);

	for (size_t i = 0; i < grid->nPts; i++) {
		grid->x[i]=xx[i];
	}
	return(0);
}

void storeGrid(const double* x, double *y, const size_t nPts){
	for(size_t i=0;i<nPts;i++){
		y[i]=x[i];
	}
}

int initializeGrid(UserGrid grid){
	grid->nPts=maxPoints(grid->basePts,
			     &grid->a,
			     &grid->w,
			     &grid->mag,
			     &grid->leftFac,
			     &grid->rightFac,
			     &grid->refineLeft,
			     &grid->refineRight);
	printf("nPts: %ld\n",grid->nPts);
	grid->leastMove=grid->w;

	grid->x    = new double [grid->nPts];
	grid->xOld = new double [grid->nPts];
	for (size_t i = 0; i < grid->nPts; i++) {
		grid->x[i]=0.0e0;
		grid->xOld[i]=0.0e0;
	}
	return(0);
}

int getGridSettings(FILE *input, UserGrid grid){

	int ier=0;

	ier=parseNumber<size_t>(input, "basePts"	, MAXBUFLEN, &grid->basePts);
	if(ier==-1)return(-1);

	ier=parseNumber<double>(input, "gridDensitySlope", MAXBUFLEN, &grid->a);
	if(ier==-1)return(-1);

	ier=parseNumber<double>(input, "fineGridHalfWidth", MAXBUFLEN, &grid->w);
	if(ier==-1)return(-1);

	ier=parseNumber<double>(input, "gridRefinement", MAXBUFLEN, &grid->mag);
	if(ier==-1)return(-1);

	ier=parseNumber<double>(input, "leftRefineFactor", MAXBUFLEN, &grid->leftFac);
	if(ier==-1)return(-1);

	ier=parseNumber<double>(input, "rightRefineFactor", MAXBUFLEN, &grid->rightFac);
	if(ier==-1)return(-1);

	ier=parseNumber<int>(input, "refineLeft"	, MAXBUFLEN, &grid->refineLeft);
	if(ier==-1)return(-1);

	ier=parseNumber<int>(input, "refineRight"	, MAXBUFLEN, &grid->refineRight);
	if(ier==-1)return(-1);

	ier=parseNumber<double>(input, "position"	, MAXBUFLEN, &grid->position);
	if(ier==-1)return(-1);

	return(0);
}

