#pragma  once  //avoid re-definition error

#ifndef SUNDIALS_DEF
#define SUNDIALS_DEF
#include <sundials/sundials_types.h>
#include <nvector/nvector_openmp.h>
#endif

#ifndef PRINT_DEF
#define PRINT_DEF
#include <string.h>	//for strings
#include <stdio.h>	//for printf,scanf
#include <stdlib.h>	//for atoi, atof
#endif

#ifndef CANTERA_DEF
#define CANTERA_DEF
#include <cantera/IdealGasMix.h>
#include <cantera/transport.h>
#include <cantera/zerodim.h>
#endif

#include "UserData.h"
// #include "interface_solver.h"
#include <ctime>
double maxTemperature(double *ydata, UserData data);

void setGasToUnity(UserData data, double temperature,
	double pressure,double* YArrayPtr);

void setGasToUnityMole(UserData data, double temperature,
	double pressure,double* XArrayPtr);

double maxGradPosition(const double* y, const size_t nt, 
	               const size_t nvar, const double* x, size_t nPts);

int maxGradIndex(const double* y, const size_t nt, 
	         const size_t nvar, const double* x, size_t nPts);

double maxCurvPosition(const double* y, const size_t nt, 
	               const size_t nvar, const double* x, size_t nPts);

int maxCurvIndex(const double* y, const size_t nt, 
	               const size_t nvar, const double* x, size_t nPts);

double isothermPosition(const double* y, const double T, const size_t nt, 
	                const size_t nvar, const double* x, const size_t nPts);

//int setAlgebraicVariables(N_Vector *id,UserData data);

inline double calc_area(double x,int* i);

void updateSolution(double* y, double* ydot, const size_t nvar,
		    const double xOld[],const double xNew[],const size_t npts);

//void readInitialCondition(FILE* input, double* ydata, const size_t nvar, const size_t nr, const size_t nPts);
void readInitialCondition(FILE* input, double* ydata, const size_t nvar, const size_t nPts, UserData data);
void readInitialCondition1(FILE* input,FILE* input1, double* ydata, const size_t nvar, const size_t nlPts, const size_t nPts, UserData data);

void readHeatFile(FILE* input, double* ydata, const size_t nvar, const size_t nr, const size_t nPts);

double systemMass(double* ydata, UserData data);

int initializePsiGrid(double* ydata, double* psidata, UserData data);

int setInitialCondition(N_Vector* y, 
	    		N_Vector* ydot, 
            		UserData data);


//inline void setGas(UserData data, double *ydata, size_t gridPoint, size_t phase);

void getTransport(UserData data, 
		  double *ydata, 
		  size_t gridPoint,
		  double *rho,
		  double *lambda,
		  double *YV);

void getInterfaceTransportWithState(UserData data,
					const double TLeft,
					const double TRight,
					const double P,
					const double* YLeft,
					const double* YRight,
					const double deltaR,
					double *rho,
					double *lambda,
					double YV[]);

int fun(double t, 
	    N_Vector y, 
	    N_Vector ydot, 
            void *user_data);

// void trackFlameOH(N_Vector y,UserData data);
// void trackFlame(N_Vector y,UserData data);
// size_t BathGasIndex(UserData data);
// size_t oxidizerIndex(UserData data);

//inline double Qdot(double* t, 
//	           double* x, 
//	           double* ignTime, 
//	           double* kernelSize, 
//	           double* maxQdot);

//double Qrad(UserData data,
//            double *ydata,
//            size_t gridPoint);

void getR(N_Vector y, UserData data);
void getR(double* ydata, UserData data);
void printSpaceTimeHeader(UserData data, FILE* output );
void printSpaceTimeOutput(double t, N_Vector* y, FILE* output, UserData data);
void printSpaceTimeOutput(double t, double* ydata, FILE* output, UserData data);
void printSpaceTimeRates(double t, N_Vector ydot, UserData data);
//void printGlobalHeader(UserData data);
//void printGlobalVariables(double t, N_Vector* y, N_Vector* ydot, UserData data);
void printSpaceTimeOutputInterpolated(double t, N_Vector y, UserData data);

void writeRestart(double t, N_Vector* y, FILE* output, UserData data);
void readRestart(N_Vector* y, N_Vector* ydot, FILE* input, UserData data);
void isobaricAdvance(double* ydata, UserData data);
void isobaricAdvance2(double* ydata, UserData data);

//
//struct isobaricReactor {
//    void *obj;
//};
//
//struct reactorNet {
//    void *obj;
//};
//
////typedef struct idealGas idealGas_t;
////typedef struct transport transport_t;
////typedef struct isochoricReactor isochoricReactor_t;
//typedef struct isobaricReactor isobaricReactor_t;
//typedef struct reactorNet reactorNet_t;
//
////idealGas_t *idealGasCreate(const char* mech);
////transport_t *transportCreate(idealGas_t *gas);
////isochoricReactor_t *isochoricReactorCreate();
//isobaricReactor_t *isobaricReactorCreate();
//reactorNet_t *reactorNetCreate();
//
////void idealGasDestroy(idealGas_t *gas);
////void isochoricReactorDestroy(isochoricReactor_t *r);
//void isobaricReactorDestroy(isobaricReactor_t *r);
//void reactorNetDestroy(reactorNet_t *sim);
//
////void addGasToIsochoricReactor(idealGas_t *gas, isochoricReactor_t *r);
////void addGasToIsobaricReactor(idealGas_t *gas, isobaricReactor_t *r);
////void setIsochoricReactorState(isochoricReactor_t *r, double* y);
//void setIsobaricReactorState(isobaricReactor_t *r, double* y);
////void getIsochoricReactorState(isochoricReactor_t *r, double* y);
//void getIsobaricReactorState(isobaricReactor_t *r, double* y);
////void addIsochoricReactorToNet(isochoricReactor_t *r, reactorNet_t *sim);
//void addIsobaricReactorToNet(isobaricReactor_t *r, reactorNet_t *sim);
//
//void reactorNetAdvance(reactorNet_t *sim, double t);
////void reactorNetReinitialize(reactorNet_t *sim);
//void reactorNetSetTolerances(reactorNet_t *sim, double rtol, double atol);
//void setNetTime(reactorNet_t *sim, double t);
//
///* modified version for LTORC */
//void addGasToIsobaricReactor(Cantera::IdealGasMix *gas, isobaricReactor_t *r);

void addGasToIsobaricReactor(UserData  data);
void addIsobaricReactorToNet(UserData data);
void reactorNetSetTolerances(UserData data);
void setIsobaricReactorState(UserData data, double* y);
void setNetTime(UserData data, double t);
void reactorNetAdvance(UserData data,double t);
void getIsobaricReactorState(UserData data,double *y);

void debugFuncReadInitialFile(double* ydata, UserData data);
int initializePsiEtaGrid(double* ydata, double* psidata, UserData data);
void debugFuncInitialPsiEtaGrid(double* psidata, UserData data);
int initializeRGridNew(double* rNew, double* ydata, double* psidata, UserData data);
void debugYdata(double* ydata, UserData data);
void debugRdata(UserData data);
int funNew(double t, N_Vector y, N_Vector ydot, void *user_data);
void getRNew(double* ydata, UserData data);
void getInterfaceTransport(UserData data,
                           double *ydata,
                           double *rho,
                           double *lambda,
                           double YV[]);

void updateInterfaceParameter(double* ydata, UserData data);
double heptaneVaporPressure(double temperature);
void updateInterfaceMassFracArray(double temperature, double PAmbience, double* massFracArrayOld,
                                  double* massFracArrayNew, double* MWArray, const size_t nsp, size_t fuelIndex);
void updateInterfaceState(double* ydata, UserData data, double delta_t);
void printSpaceTimeOutputNew(double t,double* ydata,FILE* output,UserData data);
void updateDropletMass(double* ydata, UserData data, double delta_t);
void printDropletGlobalHeader(UserData data,FILE* output);
void printDropletGlobalOutput(UserData data,FILE* output,double t);
double heptaneLatentHeat(double temperature);
void getTransportNew(UserData data,
                     double *ydata,
                     size_t gridPoint,
                     double *rho,
                     double *lambda,
                     double YV[]);
void getInterfaceMassFlux(double TLeft,double TRight,double P,double YArrayLeft[],double YArrayRight[],double deltaR,double YV[]);
double computeDerivative(const std::vector<double>& x, const std::vector<double>& y, double targetX);
double octaneVaporPressure(double temperature);
double octaneLatentHeat(double T);

void updateInterfaceState(double* ydata, UserData data, double delta_t);
void updateInterfaceCell(double* ydata, UserData data, double delta_t);