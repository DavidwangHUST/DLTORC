#ifndef MACRO_DEF
#define MACRO_DEF
//#define SUNDIALS_DOUBLE_PRECISION 1
//#define SUNDIALS_SINGLE_PRECISION 1

#define ZERO    RCONST(0.0)
#define HALF    RCONST(0.5)
#define ONE     RCONST(1.0)
#define TWO     RCONST(2.0)
#define THREE   RCONST(3.0)
#define FOUR    RCONST(4.0)
#define TEN     RCONST(10.0)
#define KB_eV   RCONST(8.6173e-5) // eV/K
#define KB_J    RCONST(1.380649e-23) // J/K
#define EC      RCONST(1.6022e-19) // C
#define EM      RCONST(9.1094e-31) // kg
#define PI      RCONST(3.14159)

#define GAS 0
#define TSP 1

/* In order to keep begin the index numbers from 1 instead of 0, we define
 * macros here. Also, we define macros to ease referencing various variables in
 * the sundials nvector.
 */
//#define psi(i)   psidata[i-1]
#define psi(i)  psidata[data->nlpts+i-1]
#define eta(i)  psidata[i-1]

//#define leftGhostR data->Rdata[0]
#define liquidR(i) data->Rdata[i-1]
#define interfaceR data->Rdata[data->nlpts]
#define gasR(i) data->Rdata[data->nlpts+i-1]
//#define rightGhostR data->Rdata[data->nlpts+data->npts+2]

#define T(i)   ydata[((i-1)*data->nvar)+data->nt]
#define Y(i,k) ydata[((i-1)*data->nvar)+data->ny+k-1]
#define P(i)   ydata[((i-1)*data->nvar)+data->np]
#define R(i)   data->Rdata[i-1]

#define Tdot(i)   ydotdata[((i-1)*data->nvar)+data->nt]
#define Ydot(i,k) ydotdata[((i-1)*data->nvar)+data->ny+k-1]
#define Pdot(i)   ydotdata[((i-1)*data->nvar)+data->np]

#define Yav(i) Yav[i-1]
#define YAvg(i) YAvg[i-1]
#define YVmhalf(i) YVmhalf[i-1]
#define YVphalf(i) YVphalf[i-1]
#define X(i) X[i-1]
#define Xp(i) Xp[i-1]
#define Xgradhalf(i) Xgradhalf[i-1]
#define XLeft(i) XLeft[i-1]
#define XRight(i) XRight[i-1]
#define gradX(i) gradX[i-1]
#define wdot(i) wdot[i-1]
#define enthalpy(i) enthalpy[i-1]
#define energy(i) energy[i-1]
#define Cp(i) Cp[i-1]

#define atolT(i)   atolvdata[((i-1)*data->nvar)+data->nt]
#define atolY(i,k)   atolvdata[((i-1)*data->nvar)+data->ny+k-1]
// #define atolR(i)   atolvdata[((i-1)*data->nvar)+data->nr]
#define atolP(i)   atolvdata[((i-1)*data->nvar)+data->np]

#define constraintsY(i,k)   constraintsdata[((i-1)*data->nvar)+data->ny+k-1]

#endif
