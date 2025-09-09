#include "parse.h"

#ifndef GSL_DEF
#define GSL_DEF
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#endif

#ifndef GRID_DEF
#define GRID_DEF

typedef struct gridTag{

	size_t basePts;
	size_t nPts;
	double position;
	double leastMove;
	double a;
	double w;
	double mag;
	double leftFac;
	double rightFac;
	int refineLeft;
	int refineRight;
	double* x;
	double* xOld;

} *UserGrid;

inline double l(const double* x, 
	 	const double* a,
	 	const double* w,
		const double* fac,
		const int* refineLeft);

inline double r(const double* x, 
		const double* a, 
		const double* w,
		const double* fac,
		const int* refineRight);

inline double f(const double* x,
		const double* a,
		const double* c,
		const double* w);

inline double g(const double* x, 
		const double* a, 
		const double* c, 
		const double* w);

inline double rho(const double* x,
	   	  const double* a,
	   	  const double* c,
	   	  const double* w,
	   	  const double* mag,
	   	  const double* leftFac,
	   	  const double* rightFac,
		  const int* refineLeft,
		  const int* refineRight);

size_t maxPoints(const size_t basePts,
	   	 const double* a,
	   	 const double* w,
	   	 const double* mag,
	   	 const double* leftFac,
	   	 const double* rightFac,
	   	 const int* refineLeft,
	   	 const int* refineRight);


double safePosition(double c, double w);

int reGrid(UserGrid grid, double position);

void storeGrid(const double* x, double *y, const size_t nPts);

int initializeGrid(UserGrid grid);

int getGridSettings(FILE *input, UserGrid grid);

#endif
