#ifndef GSL_DEF
#define GSL_DEF
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#endif
#include "residue.h"
#include "macros.h"
#include <cmath>
#include <stdio.h>
#include "timing.hpp"
#include "math.h"

#include <vector>
#include <omp.h>
#include <memory>
#include <cstdlib>

double maxTemperature(double *ydata, UserData data) {
    double Tmax = ZERO;
    for (size_t II = 1; II <= data->nlpts+data->npts; II++) {
        if (T(II) > Tmax) {
            Tmax = T(II);
        }
    }
    return Tmax;
}

double maxGradPosition(const double* y, const size_t nt, 
	               const size_t nvar, const double* x, size_t nPts){
	double maxGradT=0.0e0;
	double gradT=0.0e0;
	double pos=0.0e0;
	size_t j,jm;
	for (size_t i = 1; i <nPts; i++) {
		j=i*nvar+nt;
		jm=(i-1)*nvar+nt;
		gradT=fabs((y[j]-y[jm])/(x[i]-x[i-1]));
		if (gradT>=maxGradT) {
			maxGradT=gradT;

			pos=x[i];
		}
	}
	return(pos);
}

int maxGradIndex(const double* y, const size_t nt, 
	         const size_t nvar, const double* x, size_t nPts){
	double maxGradT=0.0e0;
	double gradT=0.0e0;
	int pos=0.0e0;
	size_t j,jm;
	for (size_t i = 1; i <nPts; i++) {
		j=i*nvar+nt;
		jm=(i-1)*nvar+nt;
		gradT=fabs((y[j]-y[jm])/(x[i]-x[i-1]));
		if (gradT>=maxGradT) {
			maxGradT=gradT;
			pos=i;
		}
	}
	return(pos);
}

double maxCurvPosition(const double* y, const size_t nt, 
	               const size_t nvar, const double* x, size_t nPts){
	double maxCurvT=0.0e0;
	double gradTp=0.0e0;
	double gradTm=0.0e0;
	double curvT=0.0e0;
	double dx=0.0e0;
	double pos=0.0e0;
	size_t j,jm,jp;
	for (size_t i = 1; i <nPts-1; i++) {
		j=i*nvar+nt;
		jm=(i-1)*nvar+nt;
		jp=(i+1)*nvar+nt;
		gradTp=fabs((y[jp]-y[j])/(x[i+1]-x[i]));
		gradTm=fabs((y[j]-y[jm])/(x[i]-x[i-1]));
		dx=0.5e0*((x[i]+x[i+1])-(x[i-1]+x[i]));
		curvT=(gradTp-gradTm)/dx;
		if (curvT>=maxCurvT) {
			maxCurvT=curvT;
			pos=x[i];
		}
	}
	return(pos);
}

int maxCurvIndex(const double* y, const size_t nt, 
	               const size_t nvar, const double* x, size_t nPts){
	double maxCurvT=0.0e0;
	double gradTp=0.0e0;
	double gradTm=0.0e0;
	double curvT=0.0e0;
	double dx=0.0e0;
	int pos=0;
	size_t j,jm,jp;
	for (size_t i = 1; i <nPts-1; i++) {
		j=i*nvar+nt;
		jm=(i-1)*nvar+nt;
		jp=(i+1)*nvar+nt;
		gradTp=fabs((y[jp]-y[j])/(x[i+1]-x[i]));
		gradTm=fabs((y[j]-y[jm])/(x[i]-x[i-1]));
		dx=0.5e0*((x[i]+x[i+1])-(x[i-1]+x[i]));
		curvT=(gradTp-gradTm)/dx;
		if (curvT>=maxCurvT) {
			maxCurvT=curvT;
			pos=i;
		}
	}
	return(pos);
}

double isothermPosition(const double* y, const double T, const size_t nt, 
	                const size_t nvar, const double* x, const size_t nPts){
	double pos=x[nPts-1];
	size_t j;
	for (size_t i = 1; i <nPts; i++) {
		j=i*nvar+nt;
		if (y[j]<=T) {
			pos=x[i];
			break;
		}
	}
	return(pos);
}

void updateSolution(double* y, double* ydot, const size_t nvar,
		    const double xOld[],const double xNew[],const size_t nPts){

	double ytemp[nPts],ydottemp[nPts];

	gsl_interp_accel* acc;
	gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, nPts);

	gsl_interp_accel* accdot;
	gsl_spline* splinedot;
	accdot    = gsl_interp_accel_alloc();
	splinedot = gsl_spline_alloc(gsl_interp_steffen, nPts);

	for (size_t j = 0; j < nvar; j++) {

		for (size_t i = 0; i < nPts; i++) {
			ytemp[i]=y[j+i*nvar];
			ydottemp[i]=ydot[j+i*nvar];
		}

		gsl_spline_init(spline,xOld,ytemp,nPts);
		gsl_spline_init(splinedot,xOld,ydottemp,nPts);

		for (size_t i = 0; i < nPts; i++) {

		      	y[j+i*nvar]=gsl_spline_eval(spline,xNew[i],acc);
		      	ydot[j+i*nvar]=gsl_spline_eval(splinedot,xNew[i],accdot);
		}
	}

	//Exploring "fixing" boundary conditions:
	//for (size_t j = 1; j < nvar; j++) {
	//	//printf("%15.6e\t%15.6e\n", y[j],y[j+nvar]);
	//	y[j]=y[j+nvar];
	//	//y[j+(nPts-1)*nvar]=y[j+(nPts-2)*nvar];
	//	//ydot[j+nvar]=ydot[j];
	//}
	//y[0]=0.0e0;
	
	gsl_interp_accel_free(acc);
	gsl_spline_free(spline);

	gsl_interp_accel_free(accdot);
	gsl_spline_free(splinedot);
}

inline double calc_area(double x,int* i){
	switch (*i) {
		case 0:
			return(ONE);
		case 1:
			return(x);
		case 2:
			return(x*x);
		default:
			return(ONE);
	}
}

void readInitialCondition(FILE* input, double* ydata, const size_t nvar, const size_t nPts, UserData data){

   //TODO Needs testing
	FILE* output;output=fopen("test.dat","w");

	size_t bufLen=10000;
	size_t nRows=0;
	size_t nColumns=nvar+1; //"1" represent spatial coordinate "R"

	char buf[bufLen];
	char buf1[bufLen];
	char comment[1];
	char *ret;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)!=0){
			nRows++;
		}
	}
	rewind(input);

	printf("nRows: %ld\n", nRows);
	double y[nRows*nColumns];

	size_t i=0;
	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
		}
		else{
			ret=strtok(buf," \t");
			size_t j=0;
			y[i*nColumns+j]=(double)(atof(ret));
			j++;
			while(ret!=NULL){
				ret=strtok(NULL," \t");
				if(j<nColumns){
					y[i*nColumns+j]=(double)(atof(ret));
				}
				j++;
			}
			i++;
		}
	}

	for (i = 0; i < nRows; i++) {
		for (size_t j = 0; j < nColumns; j++) {
			fprintf(output, "%15.6e\t",y[i*nColumns+j]);
		}
		fprintf(output, "\n");
	}
	fclose(output);

	double xOld[nRows],xNew[nPts],ytemp[nPts];

	for (size_t j = 0; j < nRows; j++) {
        double temp = ZERO;
        temp = y[j*nColumns] ;
		xOld[j]=temp;
	}

	double dx=(xOld[nRows-1] - xOld[0]) /((double)(nPts)-1.0e0);
	for (size_t j = 0; j < nPts; j++) {
		xNew[j]=data->Rd+(double)j*dx;
	}

	gsl_interp_accel* acc;
	gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, nRows);

	for (size_t j = 0; j < nColumns; j++) {

		for (size_t k = 0; k < nRows; k++) {
			ytemp[k]=y[j+k*nColumns];
		}

		gsl_spline_init(spline,xOld,ytemp,nRows);

      if(j == 0){
		   for (size_t k = 0; k < nPts; k++) {
               data->Rdata[k] = xNew[k];
//		      data->Rdata[k]=gsl_spline_eval(spline,xNew[k],acc);
		   }
      }
      else{
		   for (size_t k = 0; k < nPts; k++) {
		      ydata[j-1+k*(nColumns-1)]=gsl_spline_eval(spline,xNew[k],acc);
		   }
      }
	}

	gsl_interp_accel_free(acc);
	gsl_spline_free(spline);

}

void readInitialCondition1(FILE* input,FILE* input1, double* ydata, const size_t nvar, const size_t nlPts, const size_t nPts, UserData data){

    //TODO Needs testing
    FILE* output;output=fopen("test.dat","w");
    FILE* output1;output1=fopen("test1.dat","w");

    size_t bufLen=10000;
    size_t nRows=0;
    size_t nRows1=0;
    size_t nColumns=nvar+1; //"1" represent spatial coordinate "R"

    char buf[bufLen];
    char buf1[bufLen];
    char comment[1];
    char *ret;

    while (fgets(buf,bufLen, input)!=NULL){
        comment[0]=buf[0];
        if(strncmp(comment,"#",1)!=0){
            nRows++;
        }
    }

    rewind(input);

    printf("Initialcondition data nRows: %ld\n", nRows);

    double y[nRows*nColumns];
//    double yGlobal[(nRows+nRows1)*nColumns];
    size_t i=0;

    while (fgets(buf,bufLen, input)!=NULL){
        comment[0]=buf[0];
        if(strncmp(comment,"#",1)==0){
        }
        else{
            ret=strtok(buf," \t");
            size_t j=0;
            y[i*nColumns+j]=(double)(atof(ret));
            j++;
            while(ret!=NULL){
                ret=strtok(NULL," \t");
                if(j<nColumns){
                    y[i*nColumns+j]=(double)(atof(ret));
                }
                j++;
            }
            i++;
        }
    }


    //print y[] to output file(i.e., test.out) to validate the function
    for (i = 0; i < nRows; i++) {
        for (size_t j = 0; j < nColumns; j++) {
            fprintf(output, "%15.6e\t",y[i*nColumns+j]);
        }
        fprintf(output, "\n");
    }
    fclose(output);


//    double xOld[nRows],xNew[nPts],ytemp[nPts];
    double xOld[nRows],xNewTemp[nPts+1],xNew[nPts],ytemp[nRows];

    //assemble the xOld[] array : index range:[0,nRow-1]
    for (size_t j = 0; j < nRows; j++) {
        double temp = ZERO;
        temp = y[j*nColumns] ;
        xOld[j]=temp;
    }

    // assemble the xNew[] array
//    double dx=(xOld[nRows-1] - xOld[0]) /((double)(nPts)-1.0e0);
    double dx=(xOld[nRows-1] - xOld[0]) /(double)nPts; //
    for (size_t j = 0; j < nPts; j++) {
        xNewTemp[j]=xOld[0]+(double)j*dx;
//        xNew[j] = xOld[0] + (double)(j+1)*dx; // only include the gas phase internal grid points
    }
    xNewTemp[nPts] = xOld[nRows-1];
    for (size_t j=0; j<nPts;j++){
        xNew[j] = xNewTemp[j+1] ;// assign to xNew[] array for solver calculation,
    }

    interfaceR = xOld[0]; // = Rd
    // assign gas phase interface array
    for (size_t ii = 0; ii < nColumns; ii++) {
        data->interfaceGasCellArr[ii] = y[ii];// sequence: R,T,Ys,P
//        data->rightGhostCellArr[ii] = y[(nRows - 1) * nColumns + ii];
    }

    //allocate gsl interpolation method
    gsl_interp_accel* acc;
    gsl_spline* spline;
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_steffen, nRows);


    for (size_t j = 0; j < nColumns; j++) {
        for (size_t k = 0; k < nRows; k++) {
            ytemp[k]=y[j+k*nColumns]; //ytemp: temporary y: radius, temperature , mass fraction and pressure
        }

        gsl_spline_init(spline,xOld,ytemp,nRows);

        if(j == 0){
            for (size_t k = 0; k <= nPts; k++) {
//                data->Rdata[k] = xNew[k];
//                gasR(k+1) = xNew[k]; //book-keeping
                gasR(k+1) = xNewTemp[k];
//		      data->Rdata[k]=gsl_spline_eval(spline,xNew[k],acc);
            }
        }
        else{
            for (size_t k = 0; k < nPts; k++) {
//                ydata[j-1+k*(nColumns-1)]=gsl_spline_eval(spline,xNew[k],acc);
                ydata[j - 1 + (k + nlPts) * (nColumns - 1)] = gsl_spline_eval(spline, xNew[k], acc);
            }
        }
    }

    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);

    // now we interpolate the liquid phase points
    while (fgets(buf,bufLen, input1)!=NULL){
        comment[0]=buf[0];
        if(strncmp(comment,"#",1)!=0){
            nRows1++;
        }
    }
    rewind(input1);
    printf("Liquid initialCondition data nRows: %ld\n", nRows1);

    double y1[nRows1*nColumns];

    i = 0 ;
    //load y1[] with liquid initialcondition file data
    while (fgets(buf,bufLen, input1)!=NULL){
        comment[0]=buf[0];
        if(strncmp(comment,"#",1)==0){
        }
        else{
            ret=strtok(buf," \t");
            size_t j=0;
            y1[i*nColumns+j]=(double)(atof(ret));
            j++;
            while(ret!=NULL){
                ret=strtok(NULL," \t");
                if(j<nColumns){
                    y1[i*nColumns+j]=(double)(atof(ret));
                }
                j++;
            }
            i++;
        }
    }

    //print y[] to output1 file(i.e., test1.out) to validate the function
    for (i = 0; i < nRows1; i++) {
        for (size_t j = 0; j < nColumns; j++) {
            fprintf(output1, "%15.6e\t",y1[i*nColumns+j]);
        }
        fprintf(output1, "\n");
    }
    fclose(output1);

    double xOld1[nRows1],xNewTemp1[nlPts+1],xNew1[nlPts],ytemp1[nRows1];

    //assemble the xOld1[] array
    for (size_t j = 0; j < nRows1; j++) {
        double temp = ZERO;
        temp = y1[j*nColumns] ;
        xOld1[j]=temp;
    }

    //assemble the xNew[] array
//    double dx=(xOld[nRows-1] - xOld[0]) /((double)(nPts)-1.0e0);
    double dx1=(xOld1[nRows1-1] - xOld1[0]) / (double)nlPts;
    for (size_t j = 0; j < nlPts; j++) {
//        xNew[j]=data->Rd+(double)j*dx;
        xNewTemp1[j] = xOld1[0] + (double)j*dx1; // only include the liquid phase internal grid points
    }
    xNewTemp1[nlPts] = xOld1[nRows1-1] ;
    for(size_t j=0;j<nlPts;j++) {
        xNew1[j] = xNewTemp1[j]; // extract the first nlPts points for CVODE
    }


    // assign liquid phase interface array
    for (size_t ii = 0; ii < nColumns; ii++) {
//        data->leftGhostCellArr[ii] = y1[ii];
        data->interfaceLiquidCellArr[ii] = y1[(nRows1 - 1) * nColumns + ii];
    }

    //allocate gsl interpolation method
    gsl_interp_accel* acc1;
    gsl_spline* spline1;
    acc1 = gsl_interp_accel_alloc();
    spline1 = gsl_spline_alloc(gsl_interp_steffen, nRows1);

    for (size_t j = 0; j < nColumns; j++) {

        for (size_t k = 0; k < nRows1; k++) {
            ytemp1[k]=y1[j+k*nColumns]; //ytemp: temporary y: radius, temperature and mass fraction
        }

        gsl_spline_init(spline1,xOld1,ytemp1,nRows1);

        if(j == 0){
            for (size_t k = 0; k < nlPts; k++) {
//                data->Rdata[k] = xNew[k];
                liquidR(k+1) = xNew1[k]; //book-keeping
//		      data->Rdata[k]=gsl_spline_eval(spline,xNew[k],acc);
            }
        }
        else{
            for (size_t k = 0; k < nlPts; k++) {
//                ydata[j-1+k*(nColumns-1)]=gsl_spline_eval(spline,xNew[k],acc);
                ydata[j - 1 + k * (nColumns - 1)] = gsl_spline_eval(spline1, xNew1[k], acc1);
            }
        }
    }

    // DEBUG:need to print out the ydata and Rdata
//    debugFuncReadInitialFile(ydata,data);

    gsl_interp_accel_free(acc1);
    gsl_spline_free(spline1);
}


double systemMass(double* ydata,UserData data){
	double mass=0.0e0;
	double rho;
	for (size_t i = 2; i <=data->npts; i++) {
  		gas->setState_TPY(T(i), P(i), &Y(i,1));
		rho=gas->density();
		//psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(HALF*(R(i)+R(i-1)),&m);
		mass+=rho*(R(i)-R(i-1))*calc_area(R(i),&data->metric);
	}
	return(mass);
}

/*update the spatial coordinate R as well as ydata array(T,Ys,P) */
int initializeRGrid(double* rNew, double* ydata, double* psidata, UserData data){
    //Convert spatial cordinate R to psi
    size_t nvar = data->nvar;
    size_t npts = data->npts;
    double rOld[npts];
    double ytemp[npts];

	gsl_interp_accel* acc;
	gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, npts);

	for (size_t j = 1; j <= npts; j++) {
		rOld[j-1]=R(j);
	}

	for (size_t j = 0; j < nvar; j++) {

		for (size_t k = 0; k < npts; k++) {
			ytemp[k]=ydata[j+k*nvar];
		}

		gsl_spline_init(spline,rOld,ytemp,npts);

		for (size_t k = 0; k < npts; k++) {
		  	ydata[j+k*nvar]=gsl_spline_eval(spline,rNew[k],acc);
		}
	}

	for (size_t j = 1; j <= npts; j++) {
		R(j)=rNew[j-1];
	}
	
	gsl_interp_accel_free(acc);
	gsl_spline_free(spline);

//    //DEBUG
//    //last check for ydata
//    // Open a file in write mode
//    FILE* outYdataFile = fopen("checkYdata.txt", "w");
//    if (outYdataFile == nullptr) {
//        printf(("allocate outYdataFile failed ! \n"));
//        return 1;
//    }
//    // Write array data to file using fprintf
//    for (int i = 0; i < npts; ++i) {
//        for (int j = 0; j < nvar; ++j) {
//            fprintf(outYdataFile, "%15.6e\t", ydata[i*nvar+j]);  // Print each element with 2 decimal places
//        }
//        fprintf(outYdataFile, "\n");                     // Newline after each row
//    }
//
//    // Close the file
//    fclose(outYdataFile);



	double rho;
	/*Create a psi grid that corresponds CONSISTENTLY to the spatial grid
	 * "R" created above. Note that the Lagrangian variable psi has  units
	 * of kg. */
	psi(1)=ZERO;
	for (size_t i = 2; i <=data->npts; i++) {
  		gas->setState_TPY(T(i), P(i), &Y(i,1));
		rho=gas->density();
		//psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(HALF*(R(i)+R(i-1)),&data->metric);
		psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(R(i),&data->metric);
	}

	/*The mass of the entire system is the value of psi at the last grid
	 * point. Normalize psi by this mass so that it varies from zero to
	 * one. This makes psi dimensionless. So the mass needs to be
	 * multiplied back in the approporiate places in the governing
	 * equations so that units match.*/
	data->mass=psi(data->npts);
	for (size_t i = 1; i <=data->npts; i++) {
		psi(i)=psi(i)/data->mass;
	}
	return(0);
}

/*update the spatial coordinate R as well as ydata array(T,Ys,P) */
int initializeRGridNew(double* rNew, double* ydata, double* psidata, UserData data){
    //Convert spatial cordinate R to psi
    size_t nvar = data->nvar;
    size_t npts = data->npts;
    size_t nlpts = data->nlpts;
//    double rOld[npts];
    double rOld[npts+1];
    double ytemp[npts];
    double ytemp_1[npts+1];

    gsl_interp_accel* acc;
    gsl_spline* spline;
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_steffen, npts+1);

    for (size_t j = 1; j <= npts+1; j++) {
        rOld[j-1] = gasR(j);
//        rOld[j-1]=R(j);
    }

    for (size_t j = 0; j < nvar; j++) {

        for (size_t k = 0; k < npts; k++) {
//            ytemp[k]=ydata[j+k*nvar]; // not include the interface cell
            ytemp[k] = ydata[j+(k+nlpts)*nvar] ;
        }

        // re-assign using ytemp_1
        ytemp_1[0] = data->interfaceGasCellArr[j + 1];
        for (size_t k = 0; k < npts; k++) {
            ytemp_1[k + 1] = ytemp[k];
        }

        gsl_spline_init(spline,rOld,ytemp_1,npts+1);

        for (size_t k = 0; k < npts; k++) {
//            ydata[j+k*nvar]=gsl_spline_eval(spline,rNew[k],acc);
                ydata[j+(k+nlpts)*nvar] = gsl_spline_eval(spline,rNew[k+1],acc) ;
        }
    }

    for (size_t j = 1; j <= npts+1; j++) {
//        R(j)=rNew[j-1];
        gasR(j) = rNew[j-1] ; //size of rNew : nlpts+1
    }

    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);

//    //DEBUG
//    //last check for ydata
//    // Open a file in write mode
//    FILE* outYdataFile = fopen("checkYdata.txt", "w");
//    if (outYdataFile == nullptr) {
//        printf(("allocate outYdataFile failed ! \n"));
//        return 1;
//    }
//    // Write array data to file using fprintf
//    for (int i = 0; i < npts; ++i) {
//        for (int j = 0; j < nvar; ++j) {
//            fprintf(outYdataFile, "%15.6e\t", ydata[i*nvar+j]);  // Print each element with 2 decimal places
//        }
//        fprintf(outYdataFile, "\n");                     // Newline after each row
//    }
//
//    // Close the file
//    fclose(outYdataFile);

    double rho;
    /*Create a psi grid that corresponds CONSISTENTLY to the spatial grid
     * "R" created above. Note that the Lagrangian variable psi has  units
     * of kg. */
    psi(1) = ZERO;
    for (size_t i =2; i<= (data->npts +1);i++){
        size_t gasII = i+data->nlpts-1;
        gas->setState_TPY(T(gasII), P(gasII),&Y(gasII,1));
        rho = gas->density();
        psi(i) = psi(i-1)+rho*(gasR(i)-gasR(i-1))*calc_area(HALF*(gasR(i)+ gasR(i-1)),&data->metric);
    }

    /*The mass of the entire system is the value of psi at the last grid
     * point. Normalize psi by this mass so that it varies from zero to
     * one. This makes psi dimensionless. So the mass needs to be
     * multiplied back in the approporiate places in the governing
     * equations so that units match.*/
    data->mass = psi(data->npts+1) ;
    for (size_t i = 1; i <=(data->npts+1); i++) {
        psi(i)=psi(i)/data->mass;
    }

//    //DEBUG
//    debugFuncInitialPsiEtaGrid(psidata, data);
//    debugYdata(ydata,data);
//    debugRdata(data);

    return(0);
}

int initializePsiGrid(double* ydata, double* psidata, UserData data){

	double rho;
	/*Create a psi grid that corresponds CONSISTENTLY to the spatial grid
	 * "R" created above. Note that the Lagrangian variable psi has  units
	 * of kg. */
	psi(1)=ZERO;
	for (size_t i = 2; i <=data->npts; i++) {
  		gas->setState_TPY(T(i), P(i), &Y(i,1));
		rho=gas->density();
		//psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(HALF*(R(i)+R(i-1)),&data->metric);
		psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(R(i),&data->metric);
	}

	/*The mass of the entire system is the value of psi at the last grid
	 * point. Normalize psi by this mass so that it varies from zero to
	 * one. This makes psi dimensionless. So the mass needs to be
	 * multiplied back in the approporiate places in the governing
	 * equations so that units match.*/
	data->mass=psi(data->npts);
	for (size_t i = 1; i <=data->npts; i++) {
		psi(i)=psi(i)/data->mass;
	}
	return(0);
}

// Initialize eta(liquid) and psi(gas) in below function
int initializePsiEtaGrid(double* ydata, double* psidata, UserData data){

    // step1 : we assign the gas phase psi
    double rho;
    /*Create a psi grid that corresponds CONSISTENTLY to the spatial grid
     * "R" created above. Note that the Lagrangian variable psi has  units
     * of kg. */
//    psi(1)=ZERO;
//    for (size_t i = 2; i <=data->npts; i++) {
//        data->gas->setState_TPY(T(i), P(i), &Y(i,1));
//        rho=data->gas->density();
//        //psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(HALF*(R(i)+R(i-1)),&data->metric);
//        psi(i)=psi(i-1)+rho*(R(i)-R(i-1))*calc_area(R(i),&data->metric);
//    }
    psi(1) = ZERO;
    for (size_t i =2; i<= (data->npts +1);i++){
        size_t gasII = i+data->nlpts-1;
        gas->setState_TPY(T(gasII), P(gasII),&Y(gasII,1));
        rho = gas->density();
        psi(i) = psi(i-1)+rho*(gasR(i)-gasR(i-1))*calc_area(HALF*(gasR(i)+ gasR(i-1)),&data->metric);
    }

    /*The mass of the entire system is the value of psi at the last grid
     * point. Normalize psi by this mass so that it varies from zero to
     * one. This makes psi dimensionless. So the mass needs to be
     * multiplied back in the approporiate places in the governing
     * equations so that units match.*/
//    data->mass=psi(data->npts);
    data->mass = psi(data->npts+1) ;
    for (size_t i = 1; i <=(data->npts+1); i++) {
        psi(i)=psi(i)/data->mass;
    }

    // step2 : we assign the liquid phase eta
    eta(data->nlpts+1) = ZERO;
    for (size_t i = data->nlpts;i>=1;i--) {
        singleLiquidFuel->setLiquidState(T(i), P(i));
        rho = singleLiquidFuel->calculateDensity();
        eta(i)= eta(i+1) + rho*(liquidR(i+1)-liquidR(i))* calc_area(HALF*(liquidR(i)+ liquidR(i+1)),&data->metric);
    }
    data->dropletMass = eta(1) ;
    for (size_t i = data->nlpts;i>=1;i--){
        eta(i) = eta(i)/data->dropletMass;
    }

    //DEBUG
//    debugFuncInitialPsiEtaGrid(psidata,data);

    return(0);
}


int setInitialCondition(N_Vector* y, 
	    		N_Vector* ydot, 
            		UserData data){

	double* ydata;
	double* ydotdata;
	double* psidata;
	double* innerMassFractionsData;
	double f=ZERO;
	double g=ZERO;

	double perturb,rho;
	double epsilon=ZERO;
	int m,ier;
  	ydata    = N_VGetArrayPointer_OpenMP(*y);
  	ydotdata = N_VGetArrayPointer_OpenMP(*ydot);
	innerMassFractionsData =  data->innerMassFractions;

	if(data->adaptiveGrid){
  		psidata  = data->grid->xOld;
	}
	else{
  		psidata  = data->uniformGrid;
	}

	m=data->metric;

	data->innerTemperature=data->initialTemperature;
	for (size_t k = 1; k <=data->nsp; k++) {
  		innerMassFractionsData[k-1]=gas->massFraction(k-1);
	}

	/*following code snippet to be ready for data->problemType = 0 and 1*/
	//Define Grid:
	//double Rmin=1e-03*data->domainLength;
	double Rmin=0.0e0;
	double dR=(data->domainLength-Rmin)/((double)(data->npts)-1.0e0);
	double dv=(pow(data->domainLength,1+data->metric)-pow(data->firstRadius*data->domainLength,1+data->metric))/((double)(data->npts)-1.0e0);
	for (size_t i = 1; i <=data->npts; i++) {
		if(data->metric==0){
			R(i)=Rmin+(double)((i-1)*dR);
		}else{
			if(i==1){
				R(i)=ZERO;
			}else if(i==2){
				R(i)=data->firstRadius*data->domainLength;
			}else{
				R(i)=pow(pow(R(i-1),1+data->metric)+dv,1.0/((double)(1+data->metric)));
			}
		}
		T(i)=data->initialTemperature;
		for (size_t k = 1; k <=data->nsp; k++) {
  	        	Y(i,k)=gas->massFraction(k-1);	//Indexing different in Cantera
		}
  	        P(i)=data->initialPressure*Cantera::OneAtm;
	}	
	R(data->npts)=data->domainLength;


	double Tmax;
	double Tmin=data->initialTemperature;
	double w=data->mixingWidth;
	double YN2=ZERO;
	double YO2=ZERO;
	double YFuel,YOxidizer,sum;

	if(data->problemType==0){
		gas->equilibrate("HP");
		data->maxTemperature=gas->temperature();
	}
	else if(data->problemType==1){
		/*Premixed Combustion: Equilibrium products comprise ignition
		 * kernel at t=0. The width of the kernel is "mixingWidth"
		 * shifted by "shift" from the center.*/
		gas->equilibrate("HP");
		Tmax=gas->temperature();
		for (size_t i = 1; i <=data->npts; i++) {
			g=HALF*(tanh((R(i)-data->shift)/w)+ONE);	//increasing function of x
			f=ONE-g;					//decreasing function of x
			T(i)=(Tmax-Tmin)*f+Tmin;
			for (size_t k = 1; k <=data->nsp; k++) {
				Y(i,k)=(gas->massFraction(k-1)-Y(i,k))*f+Y(i,k);
			}
		}
		if(data->dirichletOuter){
			T(data->npts)=data->wallTemperature;
		}
	}
	else if(data->problemType==2){
		FILE* input = fopen("initialCondition.dat","r");
        FILE* input1 = fopen("initialLiquidCondition.dat","r");
		if(input != NULL && input1 != NULL){
//			readInitialCondition(input, ydata, data->nvar, data->npts, data);
            readInitialCondition1(input,input1,ydata,data->nvar,data->nlpts,data->npts,data);
			fclose(input);
		}
		else{
			printf("file initialCondition.dat(gas) or initialLiquidCondition.dat(liquid) not found!\n");
			return(-1);
		}
	}

//	initializePsiGrid(ydata,psidata,data);
    initializePsiEtaGrid(ydata,psidata,data);

	if(data->adaptiveGrid){
	//	if(data->problemType!=0){
	//		data->grid->position=maxGradPosition(ydata, data->nt, data->nvar,
	//						 data->grid->xOld, data->npts);
	//		//data->grid->position=maxCurvPosition(ydata, data->nt, data->nvar,
	//		//				 data->grid->xOld, data->npts);
	//	}
	//	else{
	//	}
		if(data->problemType!=3){
			data->grid->position=0.0e0;
			double x=data->grid->position+data->gridOffset*data->grid->leastMove;
			printf("New grid center:%15.6e\n",x);
			ier=reGrid(data->grid, x);
			if(ier==-1)return(-1);
			updateSolution(ydata, ydotdata, data->nvar,
			               data->grid->xOld,data->grid->x,data->npts);
			storeGrid(data->grid->x,data->grid->xOld,data->npts);
		}
	}
	else{
//        double Rg = data->Rg, dR0, rNew[data->npts];
        double Rg = data-> Rg, dR0, rNew[data->npts+1];
        size_t NN = data->npts;
//        size_t NN = data->npts-1;

        if( Rg == 1.0){
            dR0 = 1.0/((double)NN);
        }else{
            dR0 = (data->domainLength) * (pow(Rg,1.0/((double)NN-1.0)) - 1.0) / (pow(Rg,(double)(NN)/((double)NN-1)) - 1.0);
        }
        rNew[0] = data->Rd;
        for(size_t i = 1; i < data->npts; i++){
            rNew[i] = rNew[i-1] + dR0 * pow(Rg,(double)(i-1)/((double)NN-1));
        }
        rNew[data->npts] = data->Rd + data->domainLength;

        //DEBUG

        initializeRGridNew(rNew, ydata, psidata, data); // re-assign the Rdata and corresponding Ydata in gas phase
	}

	if(data->reflectProblem){
		double temp;
		int j=1;
		while (data->npts+1-2*j>=0) {
			temp=T(j);
			T(j)=T(data->npts+1-j);
			T(data->npts+1-j)=temp;
			for (size_t k = 1; k <=data->nsp; k++) {
				temp=Y(j,k);
				Y(j,k)=Y(data->npts+1-j,k);
				Y(data->npts+1-j,k)=temp;
			}
			j=j+1;
		}
	}

//	/*Floor small values to zero*/
//	for (size_t i = 1; i <=data->npts; i++) {
//		for (size_t k = 1; k <=data->nsp; k++) {
//			if(fabs(Y(i,k))<=data->massFractionTolerance){
//				Y(i,k)=0.0e0;
//			}
//		}
//	}

    /*Floor small values to zero*/
//    for (size_t i = 1; i <=data->npts; i++) {
//        for (size_t k = 1; k <=data->nsp; k++) {
//            if(fabs(Y(i,k))<=data->massFractionTolerance){
//                Y(i,k)=0.0e0;
//            }
//        }
//    }

    /*Floor small values to zero*/
    for (size_t i = 1; i <= (data->nlpts + data->npts); i++) {
        for (size_t k = 1; k <=data->nsp; k++) {
            if(fabs(Y(i,k))<=data->massFractionTolerance){
                Y(i,k)=0.0e0;
            }
        }
    }

	//Set grid to location of maximum curvature:
	if(data->adaptiveGrid){
		data->grid->position=maxCurvPosition(ydata, data->nt, data->nvar,
				  		data->grid->x, data->npts);
		ier=reGrid(data->grid, data->grid->position);
		updateSolution(ydata, ydotdata, data->nvar,
		               data->grid->xOld,data->grid->x,data->npts);
		storeGrid(data->grid->x,data->grid->xOld,data->npts);
	}

    //	/*Ensure consistent boundary conditions*/
//	T(1)=T(2);
//	for (size_t k = 1; k <=data->nsp; k++) {
//		Y(1,k)=Y(2,k);
//		Y(data->npts,k)=Y(data->npts-1,k);
//	}

    //DEBUG
//    getRNew(ydata,data);
//    debugRdata(data) ;

	return(0);
}

//inline double Qdot(double* t,
//	    	       double* x,
//	    	       UserData data,
//                   size_t gridPoint,
//		           double *ydata){
//   double qdot; //W/m3
//   double* ignTime = &data->ignTime;
//   double* maxQDot = &data->maxQDot;
//   double kernelSize = data->kernelSize;
//   double electrodeRadius = data->electrodeRadius;
//   if((*t)<=(*ignTime)){
//      if(data->heatType == 0){
//	      if(*x<=kernelSize){
//	   	   qdot = (*maxQDot);
//         }else{
//            qdot = 0.0;
//         }
//      }else if(data->heatType == 1){
//		   qdot = gsl_spline_eval(data->heatSpline,*t,data->heatAcc);
//         //printf("time: %15.6e, interp: %15.6e, ",*t,qdot);
//         qdot = qdot/(4.0/3.0 * PI * pow(kernelSize,3.0)) * exp( - PI / 4.0 * pow((*x)/kernelSize,6.0));
//         //printf("qdot: %15.6e\n",*t,qdot);
//      }
//   }else{
//      qdot = 0.0;
//   }
//
//
//
//	//if(*x<=*kernelSize){
//	//	if((*t)<=(*ignTime)){
//   //         if(data->heatType == 0){
//	//		    qdot = (*maxQDot);
//   //         }else if(data->heatType == 1){
//	//	      	qdot = gsl_spline_eval(data->heatSpline,*t,data->heatAcc);
//   //         }
//	//	}
//	//	else{
//	//		qdot=0.0e0;
//	//	}
//	//}else{
//	//	qdot=0.0e0;
//	//}
//
//    if(data->electrodeLoss && *x<=electrodeRadius){
//        //Estimate Heat losses to electrode
//        //Assume the electrode temperature is equal to the initial temperature
//        //Take the average temperature between grid point and electrode to find film temperature Tavg
//        //Get thermal conductivity, lambda, at film temperature
//        double Tavg;
//        double lambda,heatFlux,Qe;
//        double delta = data->eLossDelta;
//        double Yt[data->nsp-1];
//        size_t k;
//	    for (size_t tii = 0; tii < (data->nsp-1); tii++) {
//            k = data->k_transport_[tii];
//	    	Yt[tii]=Y(gridPoint,k+1);
//	    }
//
//        Tavg = (T(gridPoint) + data->initialTemperature)/TWO;
//  	    data->transport->setTemperature(Tavg);
//  	    data->transport->setMassFractions_NoNorm(Yt);
//  	    //data->transport->setMassFractions(Yt);
//  	    data->transport->setPressure(P(gridPoint));
//
//	    lambda=data->trmix->thermalConductivity();
//        //Compute heat flux assuming some small length scale, delta
//        heatFlux=lambda * (T(gridPoint) - data->initialTemperature)/delta;
//        //Qe represents heat loss in W/m3
//        //The total heat loss occurs at the top and bottom areas of the cell, in a direction we are not resolving
//        //Multiply heatFlux by these areas for the total heat loss
//        //Then, divide this by the volume of the cell, for cylindrical coordinates results in the following
//        //Approximate the height of the cell as the kernelSize
//        //TODO Replace 2.0e-4 with electrodeGap option
//        Qe = TWO  * heatFlux / 2.0e-4;
//
//        ////setGas(data, ydata, gridPoint, TSP);
//
//	    return(qdot - Qe);
//    }else{
//	    return(qdot);
//    }
//}

//double Qrad(UserData data,
//			double* ydata,
//			size_t gridPoint){
//	double Qrad;
//	Qrad = gsl_spline_eval(data->radSpline,T(gridPoint),data->radAcc);
//    return(Qrad);
//}

//double therCondAirPlasma(UserData data,
//			double* ydata,
//			size_t gridPoint){
//	double therCond;
//	therCond = gsl_spline_eval(data->therCondSpline,T(gridPoint),data->therCondAcc);
//    return(therCond);
//}

inline void setGas(UserData data, 
	    double *ydata, 
	    size_t gridPoint){

  	gas->setTemperature(T(gridPoint));
  	gas->setMassFractions_NoNorm(&Y(gridPoint,1));
  	gas->setPressure(P(gridPoint));
}

inline void setGasToUnity(UserData data, double temperature,
                          double pressure,double* YArrayPtr){
    gas->setTemperature(temperature);
    gas->setPressure(pressure);
    gas->setMassFractions_NoNorm(YArrayPtr);
}

// inline void setGas(UserData data, 
// 	    double *ydata, 
// 	    size_t gridPoint,
//        size_t phase){
//     try{
//         if(phase == GAS){
//   	        gas->setTemperature(T(gridPoint));
//   	        data->gas->setMassFractions_NoNorm(&Y(gridPoint,1)); 	
//   	        data->gas->setPressure(P(gridPoint));
//         }else if(phase == TSP){
//            double Ytransport[data->nsp-1];
//            size_t k;
//            for(size_t tii=0; tii < (data->nsp-1); tii++){
//               k = data->k_transport_[tii];
//               Ytransport[tii] = Y(gridPoint, k+1);
//            }
//   	        data->transport->setTemperature(T(gridPoint)); 	        
//   	        data->transport->setMassFractions_NoNorm(Ytransport); 	
//   	        //data->transport->setMassFractions(Ytransport); 	
//   	        data->transport->setPressure(P(gridPoint));
//         } else {
//           printf("setGas: a phase other than the gas or transport phase was selected!\n");
//           return;
//         }
//     }catch(Cantera::CanteraError& err){

//         printf("Printing err.dat\n");
//         FILE* errOutput;
//         errOutput = fopen("err.dat","w");
//         printSpaceTimeHeader(data, errOutput);
//         printSpaceTimeOutput(0.0, ydata, errOutput, data);
//         fclose(errOutput);
//         printf("done!\n");
//         return;
//     }
// }

void getInterfaceMassFlux(double TLeft,double TRight,double P,double YArrayLeft[],double YArrayRight[],double deltaR,double YV[]){
    size_t nsp = gas->nSpecies();
    double gradT;
    double YAvg[nsp],XArrayLeft[nsp],XArrayRight[nsp],gradX[nsp];
    double TAvg = HALF*(TLeft+TRight);

    gas->setTemperature(TLeft);
    gas->setPressure(P);
    gas->setMassFractions_NoNorm(YArrayLeft);
    gas->getMoleFractions(XArrayLeft);

    gas->setTemperature(TRight);
    gas->setPressure(P);
    gas->setMassFractions_NoNorm(YArrayRight);
    gas->getMoleFractions(XArrayRight);

    for (size_t k=1;k<=nsp;k++){
        YAvg[k-1] = YArrayLeft[k-1] + YArrayRight[k-1];
        gradX[k-1] = XArrayRight[k-1] - XArrayLeft[k-1];
    }
    gradT = (TRight-TLeft)/deltaR;
    trmix->getSpeciesFluxes(1,&gradT,nsp,gradX,nsp,YV) ;
}

void getTransport(UserData data, 
		  double *ydata, 
		  size_t gridPoint,
		  double *rho,
		  double *lambda,
		  double YV[]){

	double YAvg[data->nsp],
		 XLeft[data->nsp],
		 XRight[data->nsp],
		 gradX[data->nsp];

	setGas(data,ydata,gridPoint);
	gas->getMoleFractions(XLeft);
	setGas(data,ydata,gridPoint+1);
	gas->getMoleFractions(XRight);
                                                                             
	for (size_t k = 1; k <=data->nsp; k++) {                                          
		YAvg(k)=HALF*(Y(gridPoint,k)+
			      Y(gridPoint+1,k));                                 
		gradX(k)=(XRight(k)-XLeft(k))/
			 (R(gridPoint+1)-R(gridPoint));                       
	}                                                                    
	double TAvg = HALF*(T(gridPoint)+T(gridPoint+1));
	double gradT=(T(gridPoint+1)-T(gridPoint))/
		       (R(gridPoint+1)-R(gridPoint));                                  


  	gas->setTemperature(TAvg);
  	gas->setMassFractions_NoNorm(YAvg);
  	gas->setPressure(P(gridPoint));

	*rho=gas->density();
	*lambda=trmix->thermalConductivity();
	trmix->getSpeciesFluxes(1,&gradT,data->nsp,
				      gradX,data->nsp,YV);                
	//setGas(data,ydata,gridPoint);
}

void getTransportNew(UserData data,
                  double *ydata,
                  size_t gridPoint,
                  double *rho,
                  double *lambda,
                  double YV[]){

    double YAvg[data->nsp],
            XLeft[data->nsp],
            XRight[data->nsp],
            gradX[data->nsp];
    double dR = R(gridPoint+2)-R(gridPoint+1);

    setGas(data,ydata,gridPoint);
    gas->getMoleFractions(XLeft);
    setGas(data,ydata,gridPoint+1);
    gas->getMoleFractions(XRight);

    for (size_t k = 1; k <=data->nsp; k++) {
        YAvg(k)=HALF*(Y(gridPoint,k)+
                      Y(gridPoint+1,k));
        gradX(k)=(XRight(k)-XLeft(k))/ dR;
    }
    double TAvg = HALF*(T(gridPoint)+T(gridPoint+1));
    double gradT=(T(gridPoint+1)-T(gridPoint))/ dR;


    gas->setTemperature(TAvg);
    gas->setMassFractions_NoNorm(YAvg);
    gas->setPressure(P(gridPoint));

    *rho=gas->density();
    *lambda=trmix->thermalConductivity();
    trmix->getSpeciesFluxes(1,&gradT,data->nsp,
                            gradX,data->nsp,YV);
    //setGas(data,ydata,gridPoint);
}

void getInterfaceTransport(UserData data,
                  double *ydata,
                  double *rho,
                  double *lambda,
                  double YV[]){
    double YAvg[data->nsp],
            YLeft[data->nsp],
            XLeft[data->nsp],
            XRight[data->nsp],
            gradX[data->nsp];
    double interfaceArray[data->nvar+1] ;
    for(size_t i =0 ;i<= data->nvar;i++) {
        interfaceArray[i] = data->interfaceGasCellArr[i] ;
    }
    for(size_t k=0;k<data->nsp;k++){
        YLeft[k] = interfaceArray[2+k] ;
    }

//    setGas(data,ydata,gridPoint);
    gas->setState_TPY(interfaceArray[1],interfaceArray[data->nvar],&interfaceArray[2]);
    gas->getMoleFractions(XLeft);
    setGas(data,ydata,data->nlpts+1);
    gas->getMoleFractions(XRight);

    for (size_t k = 1; k <=data->nsp; k++) {
        YAvg(k)=HALF*(YLeft[k-1]+
                      Y(data->nlpts+1,k));
        gradX(k)=(XRight(k)-XLeft(k))/
                 (gasR(2)-gasR(1));
    }
    double TAvg = HALF*(interfaceArray[1]+T(data->nlpts+1));
//    double gradT=(T(gridPoint+1)-T(gridPoint))/
//                 (R(gridPoint+1)-R(gridPoint));
    double gradT = (T(data->nlpts+1)-interfaceArray[1]) / (gasR(2) - gasR(1));

    gas->setTemperature(TAvg);
    gas->setMassFractions_NoNorm(YAvg);
    gas->setPressure(P(data->nlpts+1));

    *rho=gas->density();
    *lambda=trmix->thermalConductivity();
    trmix->getSpeciesFluxes(1,&gradT,data->nsp,
                                  gradX,data->nsp,YV);// units: kg/(m^2*s)
    //setGas(data,ydata,gridPoint);
}

// Compute transport properties and species diffusive fluxes at the gas-side of the interface
// using an explicit interface state (TLeft, YLeft) and the first gas-cell state (TRight, YRight).
// All fluxes are returned per unit area, units of kg/(m^2*s).
static inline void getInterfaceTransportWithState(UserData data,
				 const double TLeft,
				 const double TRight,
				 const double P,
				 const double* YLeft,
				 const double* YRight,
				 const double deltaR,
				 double *rho,
				 double *lambda,
				 double YV[]){

	const size_t nsp = data->nsp;
    double YAvg[nsp], XLeft[nsp], XRight[nsp], gradX[nsp];

	// Left (interface) state
	gas->setTemperature(TLeft);
	gas->setPressure(P);
	gas->setMassFractions_NoNorm(YLeft);
	gas->getMoleFractions(XLeft);

	// Right (first gas cell) state
	gas->setTemperature(TRight);
	gas->setPressure(P);
	gas->setMassFractions_NoNorm(YRight);
	gas->getMoleFractions(XRight);

	for (size_t k = 1; k <= nsp; k++) {
		YAvg[k-1] = HALF * (YLeft[k-1] + YRight[k-1]);
		gradX[k-1] = (XRight[k-1] - XLeft[k-1]) / deltaR;
	}
	const double TAvg = HALF * (TLeft + TRight);
	const double gradT = (TRight - TLeft) / deltaR;

	gas->setTemperature(TAvg);
	gas->setMassFractions_NoNorm(YAvg);
	gas->setPressure(P);

	*rho = gas->density();
	*lambda = trmix->thermalConductivity();
	trmix->getSpeciesFluxes(1, &gradT, nsp, gradX, nsp, YV);
}

// Context for multi-variable interface solve enforcing non-penetration constraints (fuel only crosses)
struct InterfaceNPContext {
	UserData data;
	size_t nsp;
	size_t fuelIndex; // 1-based index
	double P;
	// geometry/neighboring states
	double interface_R;
	double right_R;
	double left_R;
	double left_T;
	double right_T;
	std::vector<size_t> nonFuelIdx; // 0-based species indices in x ordering
	std::vector<double> YRight;     // Y at first gas cell
	std::vector<double> MW;         // molecular weights
};

// Residual for the nonlinear system with unknowns: x[0]=T_interface; x[1..]=Y_nonfuel at interface
// Equations:
// - For each non-fuel species i: J_i + Y_i * m" = 0 (non-penetration)
// - Vapor-pressure equilibrium: X_fuel - Psat(T)/P = 0; with Y_fuel = 1 - sum(nonFuel)
static int interface_np_f(const gsl_vector* x, void* params, gsl_vector* f){
	InterfaceNPContext* ctx = static_cast<InterfaceNPContext*>(params);
	UserData data = ctx->data;
	const size_t nsp = ctx->nsp;
	const size_t fuelII = ctx->fuelIndex; // 1-based

	// Unpack unknowns
	const double T_int = gsl_vector_get(x, 0);
	const double mdot_area = gsl_vector_get(x, 1); // kg/(m^2*s)
	std::vector<double> Yint(nsp, 0.0);
	double sumNonFuel = 0.0;
	for(size_t j=0;j<ctx->nonFuelIdx.size();j++){
		double val = gsl_vector_get(x, 2 + j);
		if(val < 0.0) val = 0.0;
		Yint[ ctx->nonFuelIdx[j] ] = val;
		sumNonFuel += val;
	}
	// Fuel mass fraction by closure
	const size_t fuel0 = fuelII - 1;
	double Yfuel = ONE - sumNonFuel;
	if(Yfuel < 0.0) Yfuel = 0.0;
	if(Yfuel > 1.0) Yfuel = 1.0;
	Yint[fuel0] = Yfuel;

	// Transport and diffusive fluxes on the gas side of the interface
	double rho_g, lambda_g;
	std::vector<double> YV(nsp,0.0);
	const double deltaR_g = std::max(ctx->right_R - ctx->interface_R, 1e-12);
	getInterfaceTransportWithState(data,
		T_int, ctx->right_T, ctx->P,
		&Yint[0], &ctx->YRight[0], deltaR_g,
		&rho_g, &lambda_g, &YV[0]);

	// Non-penetration residuals for all non-fuel species
	size_t eq_idx = 0;
	for(size_t j=0;j<ctx->nonFuelIdx.size();j++){
		const size_t k0 = ctx->nonFuelIdx[j];
		const double res_np = YV[k0] + Yint[k0]*mdot_area;
		gsl_vector_set(f, eq_idx++, res_np);
	}

	// Vapor-pressure equilibrium residual in mole fraction
	double denom_x = 0.0;
	for(size_t k=0;k<nsp;k++) denom_x += Yint[k]/ctx->MW[k];
	denom_x = std::max(denom_x, 1e-300);
	double Xfuel = (Yfuel/ctx->MW[fuel0]) / denom_x;
	const double Xfuel_eq = heptaneVaporPressure(T_int)/ctx->P;
	const double res_vp = Xfuel - Xfuel_eq;
	gsl_vector_set(f, eq_idx++, res_vp);

	// Energy conservation across interface: q_l - q_g - m" * L = 0
	const double deltaR_l = std::max(ctx->interface_R - ctx->left_R, 1e-12);
	singleLiquidFuel->setLiquidState(T_int, ctx->P);
	const double lambda_l = singleLiquidFuel->calculateThermalConductivity();
	const double q_l = lambda_l * (T_int - ctx->left_T) / deltaR_l; // W/m^2
	const double q_g = lambda_g * (ctx->right_T - T_int) / deltaR_g; // W/m^2
	const double L = heptaneLatentHeat(T_int); // J/kg
	const double res_energy = q_l - q_g - mdot_area * L;
	gsl_vector_set(f, eq_idx++, res_energy);

	// Optional debug prints: set env DEBUG_IFACE=1 to enable
	static int dbg = -1;
	if (dbg < 0) {
		const char* v = std::getenv("DEBUG_IFACE");
		dbg = (v && v[0] != '\0') ? 1 : 0;
	}
	if (dbg) {
		printf("[interface_np_f] T_int=%.6e, mdot_area=%.6e, deltaR_g=%.6e, deltaR_l=%.6e\n",
			T_int, mdot_area, deltaR_g, deltaR_l);
		for (size_t j = 0; j < ctx->nonFuelIdx.size(); j++) {
			const size_t k0 = ctx->nonFuelIdx[j];
			const double res_np = YV[k0] + Yint[k0]*mdot_area;
			printf("  res_np[%zu](%s) = %.6e  (YV=%.6e, Y=%.6e)\n",
				k0+1, gas->speciesName(k0).c_str(), res_np, YV[k0], Yint[k0]);
		}
		printf("  res_vp = %.6e  (Xfuel=%.6e, Xeq=%.6e)\n", res_vp, Xfuel, Xfuel_eq);
		printf("  res_energy = %.6e  (q_l=%.6e, q_g=%.6e, L=%.6e)\n",
			res_energy, q_l, q_g, L);
	}

	return GSL_SUCCESS;
}

// Solve interface (non-penetration for non-fuel, VLE for fuel) and update interface state
static int solveInterfaceNonPenetration(UserData data, double* ydata,double delta_t){
	const size_t nsp = data->nsp;
	const size_t nlpts = data->nlpts;
	const size_t fuelII = data->dropII; // 1-based
	const double P = data->initialPressure * Cantera::OneAtm;

	InterfaceNPContext ctx;
	ctx.data = data;
	ctx.nsp = nsp;
	ctx.fuelIndex = fuelII;
	ctx.P = P;
	ctx.interface_R = data->interfaceGasCellArr[0];
	ctx.right_R = R(nlpts+1);
	if(nlpts > 1){
		ctx.left_R = R(nlpts);
	}else{
		ctx.left_R = std::max(data->interfaceGasCellArr[0] - (data->interfaceGasCellArr[0]-data->interfaceGasCellArr[1]), 0.0);
	}
	ctx.left_T = T(nlpts);
	ctx.right_T = T(nlpts+1);
	ctx.YRight.assign(nsp, 0.0);
	for(size_t k=0;k<nsp;k++) ctx.YRight[k] = Y(nlpts+1, k+1);
	ctx.MW.assign(nsp, 0.0);
	for(size_t k=0;k<nsp;k++) ctx.MW[k] = gas->molecularWeight(k);
	ctx.nonFuelIdx.clear();
	for(size_t k0=0;k0<nsp;k0++) if((k0+1) != fuelII) ctx.nonFuelIdx.push_back(k0); // 0-based

	// Unknowns: [T_interface, mdot_area, Yi_nonfuel...]
	const size_t nUnknown = 2 + ctx.nonFuelIdx.size();
	gsl_vector* x = gsl_vector_alloc(nUnknown);
	gsl_vector_set(x, 0, data->interfaceGasCellArr[1]);
	// initial guess for Yi from current interface; fallback to right cell
	double sumNF = 0.0;
	for(size_t j=0;j<ctx.nonFuelIdx.size();j++){
		const size_t k0 = ctx.nonFuelIdx[j];
		double v = data->interfaceGasCellArr[2 + k0];
		if(!(v==v) || v<0.0){ v = Y(nlpts+1, k0+1); }
		gsl_vector_set(x, 2+j, v);
		sumNF += v;
	}
	if(sumNF <= 0.0){
		for(size_t j=0;j<ctx.nonFuelIdx.size();j++){
			const size_t k0 = ctx.nonFuelIdx[j];
			gsl_vector_set(x, 2+j, Y(nlpts+1, k0+1));
		}
	}
	// Initial guess for mdot_area from current interface state
	{
		std::vector<double> Yint0(nsp,0.0);
		double sumNF0 = 0.0;
		for(size_t j=0;j<ctx.nonFuelIdx.size();j++){
			const size_t k0 = ctx.nonFuelIdx[j];
			Yint0[k0] = gsl_vector_get(x, 2+j);
			sumNF0 += Yint0[k0];
		}
		const size_t fuel0 = fuelII-1;
		Yint0[fuel0] = std::max(0.0, 1.0 - sumNF0);
		double rho_g0, lambda_g0; std::vector<double> YV0(nsp,0.0);
		const double dRg = std::max(ctx.right_R - ctx.interface_R, 1e-12);
		getInterfaceTransportWithState(data,
			data->interfaceGasCellArr[1], ctx.right_T, P,
			&Yint0[0], &ctx.YRight[0], dRg,
			&rho_g0, &lambda_g0, &YV0[0]);
		const double denom = std::max(1.0 - Yint0[fuel0], 1e-12);
		gsl_vector_set(x, 1, YV0[fuel0]/denom);
	}

	const gsl_multiroot_fsolver_type* Tsolver = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(Tsolver, nUnknown);
	gsl_multiroot_function F;
	F.f = &interface_np_f;
	F.n = nUnknown; // (nsp-1 non-fuel) + vapor pressure + energy
	F.params = &ctx;
	gsl_multiroot_fsolver_set(solver, &F, x);

	int status = GSL_CONTINUE;
	size_t iter = 0; const size_t max_iter = 100;
	while (status == GSL_CONTINUE && iter < max_iter){
		iter++;
		status = gsl_multiroot_fsolver_iterate(solver);
		if(status) break;
		status = gsl_multiroot_test_residual(solver->f, 1e-9);
	}

	if(status != GSL_SUCCESS){
		gsl_multiroot_fsolver_free(solver);
		gsl_vector_free(x);
		return -1;
	}

	// Extract solution
	gsl_vector* xsol = solver->x;
	const double T_int = gsl_vector_get(xsol, 0);
	std::vector<double> Yint(nsp, 0.0);
	double sumNonFuelSol = 0.0;
	for(size_t j=0;j<ctx.nonFuelIdx.size();j++){
		double v = gsl_vector_get(xsol, 2+j);
		if(v < 0.0) v = 0.0;
		Yint[ ctx.nonFuelIdx[j] ] = v;
		sumNonFuelSol += v;
	}
	Yint[fuelII-1] = std::max(0.0, 1.0 - sumNonFuelSol);

	// Update interface arrays: sequence R, T, Ys, P
	data->interfaceGasCellArr[1] = T_int;
	for(size_t k0=0;k0<nsp;k0++) data->interfaceGasCellArr[2+k0] = Yint[k0];

    double mdot_area = gsl_vector_get(xsol, 1) * calc_area(ctx.interface_R, &data->metric); // kg/s
    data->dropletMass -= (mdot_area) * delta_t;

	gsl_multiroot_fsolver_free(solver);
	gsl_vector_free(x);
	return 0;
}

// void getTransport(UserData data, 
// 		  double *ydata, 
// 		  size_t gridPoint,
// 		  double *rho,
// 		  double *lambda,
// 		  double YV[]){

// 	double YAvg[data->nsp-1],
// 		   XLeft[data->nsp-1],
// 		   XRight[data->nsp-1],
// 		   gradX[data->nsp-1],
//            YVtmp[data->nsp-1],
//            sigma, epsilon, mu, gamma, De, gradYe;
       
//    size_t k;

// 	setGas(data,ydata,gridPoint,TSP);
// 	data->transport->getMoleFractions(XLeft);                                     
// 	setGas(data,ydata,gridPoint+1,TSP);
// 	data->transport->getMoleFractions(XRight);                                     

//     //DEBUG
//     //if(gridPoint == 2){
//     //    std::cout << data->transport->report();
//     //    printf("PRINT REPORT\N");
//     //}
                                                                             
// 	for (size_t tii = 0; tii < (data->nsp-1); tii++) {                                          
//         k = data->k_transport_[tii]+1;
// 		YAvg[tii]=HALF*(Y(gridPoint,k)+
// 			      Y(gridPoint+1,k));                                 
// 		gradX[tii]=(XRight[tii]-XLeft[tii])/
// 			 (R(gridPoint+1)-R(gridPoint));                       
// 	}                                                                    
// 	double TAvg = HALF*(T(gridPoint)+T(gridPoint+1));
// 	double gradT=(T(gridPoint+1)-T(gridPoint))/
// 		       (R(gridPoint+1)-R(gridPoint));                                  

//   	data->transport->setTemperature(TAvg); 	        
//   	data->transport->setMassFractions_NoNorm(YAvg); 	
//   	//data->transport->setMassFractions(YAvg); 	
//   	data->transport->setPressure(P(gridPoint));

// 	*rho=data->transport->density();                                       
// 	//*lambda=data->trmix->thermalConductivity();  
//    *lambda = therCondAirPlasma(data,ydata,gridPoint);
   
//    //DEBUG
//    //if(gridPoint == 5){
//    //   printf("T=%0.5e, lambda=%0.5e, ckLabmda=%0.5e\n",T(gridPoint),*lambda,data->trmix->thermalConductivity());
//    //}

// 	data->trmix->getSpeciesFluxes(1,&gradT,data->nsp-1,                     
// 				      gradX,data->nsp-1,YVtmp);                
//    for(size_t tii = 0; tii < (data->nsp-1); tii++){
//       k = data->k_transport_[tii];
//       YV[k] = YVtmp[tii];
//       //DEBUG
//       //printf("k=%d  YV=%0.5e\n",k,YV[k]); 
//    }


//    //Evaluate electron transport
//    gradYe = (Y(gridPoint+1,data->k_e) - Y(gridPoint,data->k_e))/(R(gridPoint+1) - R(gridPoint));
//    gamma = sqrt(TWO * EC / EM);
//    epsilon = 3.0 * KB_eV * TAvg / 2.0; // average electron energy in eV
//    sigma = gsl_spline_eval(data->crossSpline,epsilon,data->crossAcc); // Cross-section in m2
//    mu = TWO * gamma * EC / THREE * sqrt(KB_J / PI) / sigma * sqrt(TAvg) / P(gridPoint); // electron mobility
//    //mu = TWO * gamma * EC / THREE * sqrt(KB_eV / PI) / sigma * sqrt(TAvg) / P(gridPoint); // electron mobility
//    De = mu * KB_J * TAvg / EC; 
//    YV[data->k_e-1] = De * gradYe * *rho;

//    //DEBUG
//    //printf("k=%d  YV=%0.5e\n",data->k_e-1,YV[data->k_e-1]); 
//    //printf("\n\n");

//    //YV[data->k_e-1] = 0.0;

//    //printf("T=%e, epsilon=%e,  sigma=%e,   mu=%e\n",TAvg,epsilon,sigma,mu);
// }

int fun(double t,
	    N_Vector y,
	    N_Vector ydot,
       void *user_data){

	/*Declare and fetch nvectors and user data:*/

  	double *ydata, *ydotdata, *psidata, *innerMassFractionsData;

	UserData data;
	data = (UserData)user_data;
	size_t npts=data->npts;
	size_t nsp=data->nsp;
	size_t k_bath = data->k_bath;

//	getR(y, data);

	ydata	= N_VGetArrayPointer_OpenMP(y);
	ydotdata = N_VGetArrayPointer_OpenMP(ydot);
	if(data->adaptiveGrid==1){
		psidata = data->grid->x;
	}else{
		psidata = data->uniformGrid;
	}

    getR(ydata,data);

	innerMassFractionsData =  data->innerMassFractions;

	/* Grid stencil:*/

	/*-------|---------*---------|---------*---------|-------*/
	/*-------|---------*---------|---------*---------|-------*/
	/*-------|---------*---------|---------*---------|-------*/
	/*-------m-------mhalf-------j-------phalf-------p-------*/
	/*-------|---------*---------|---------*---------|-------*/
	/*-------|---------*---------|---------*---------|-------*/
	/*-------|<=======dxm=======>|<=======dxp=======>|-------*/
	/*-------|---------*<======dxav=======>*---------|-------*/
	/*-------|<================dxpm=================>|-------*/

	/* Various variables defined for book-keeping and storing previously
	 * calculated values:
	 * rho		: densities at points  m, mhalf, j, p, and phalf.
	 * area		: the matric at points m, mhalf, j, p, and phalf.
	 * m 		: exponent that determines geometry;
	 * lambda	: thermal conductivities at mhalf and phalf.
	 * mdot		: mass flow rate at m, j, and p.
	 * X		: mole fractions at j and p.
	 * YV		: diffusion fluxes at mhalf and phalf.
	 * Tgrad	: temperature gradient at mhalf and phalf.
	 * Tav		: average temperature between two points.
	 * Pav		: average pressure between two points.
	 * Yav		: average mass fractions between two points.
	 * Xgradhalf	: mole fraction gradient at j.
	 * Cpb		: mass based bulk specific heat.
	 * tranTerm	: transient terms.
	 * advTerm	: advection terms.
	 * diffTerm	: diffusion terms.
	 * srcTerm	: source terms.
	 */

	double rhomhalf,  rhom, lambdamhalf, YVmhalf[nsp],
		 rho,
		 rhophalf, lambdaphalf, YVphalf[nsp],
		 Cpb, Cvb,  Cp[nsp],   wdot[nsp],   enthalpy[nsp], energy[nsp],
		 tranTerm, diffTerm,    srcTerm, advTerm,
		 area,areamhalf,areaphalf,aream,areamhalfsq,areaphalfsq;

	/*Aliases for difference coefficients:*/
	double cendfm, cendfc, cendfp;
	cendfm=cendfc=cendfp=ZERO;
	/*Aliases for various grid spacings:*/
	double dpsip, dpsiav, dpsipm, dpsim, dpsimm;
	dpsip=dpsiav=dpsipm=dpsim=dpsimm=ONE;
	double mass, mdotIn;
	double sum, sum1, sum2, sum3;

	size_t j,k;
//	int m;
    int m=data->metric;					//Unitless
	mass=data->mass;				//Units: kg
	mdotIn=data->mdot*calc_area(R(npts),&m);	//Units: kg/s

//	/*evaluate properties at j=1*************************/
    //printf("Set gas at j=1\n");
	// setGas(data,ydata,1,GAS);
	setGas(data,ydata,1);
	rhom=gas->density();
        Cpb=gas->cp_mass();       	//J/kg/K
        Cvb=gas->cv_mass();       	//J/kg/K
	aream= calc_area(R(1),&m);

	/*******************************************************************/
	/*Calculate values at j=2's m and mhalf*****************************/

    //printf("Get transport at j=1\n");
	getTransport(data, ydata, 1, &rhomhalf,&lambdamhalf,YVmhalf);
	areamhalf= calc_area(HALF*(R(1)+R(2)),&m);
	areamhalfsq= areamhalf*areamhalf;

	/*Calculate the droplet vaporization rate: kg/s */
	data->Mdot = YVmhalf(data->dropII) * areamhalf / (1 - Y(1, data->dropII)) ;

	/*******************************************************************/

	/*Fill up res with left side (center) boundary conditions:**********/
	/*We impose zero fluxes at the center:*/

	/*Mass:*/
	//Rdot(1) = 0.0;

	/*Energy:*/

	if (data->dirichletInner){
		Tdot(1) = 0.0;
		// TEST:try to fix the boundary mass fraction
		for (k = 1; k <= nsp; k++){
			Ydot(1, k) = 0.0;
		}
	}
	//else{
	//	Tres(1)=T(2)-T(1);
	//	//Tres(1)=Tdot(1) - (Pdot(1)/(rhom*Cpb))
	//	//	+(double)(data->metric+1)*(rhomhalf*lambdamhalf*areamhalfsq*(T(2)-T(1))/psi(2)-psi(1));
	//}

	/*Species:*/
	sum=ZERO;
	for (k = 1; k <=nsp; k++) {
		if(k!=k_bath){
			if(fabs(mdotIn)>1e-14){
               //TODO Account for mdotin
		      	//	Ydot(1,k)=innerMassFractionsData[k-1]-
					//  Y(1,k)-
					//  (YVmhalf(k)*areamhalf)/mdotIn;
			}
		}
	}
	//Yres(1,k_bath)=ONE-sum-Y(1,k_bath);


	/*Pressure:*/
   //TODO implement constant volume case
	Pdot(1)=data->dPdt;

	/*Fill up res with governing equations at inner points:*************/
	for (j = 2; j < npts; j++) {

		/*evaluate various mesh differences*///
        	dpsip =        (psi(j+1) - psi(j)  )*mass;
        	dpsim =        (psi(j)   - psi(j-1))*mass;
        	dpsiav =  HALF*(psi(j+1) - psi(j-1))*mass;
        	dpsipm =       (psi(j+1) - psi(j-1))*mass;
		/***********************************///

		/*evaluate various central difference coefficients*/
        	cendfm = - dpsip / (dpsim*dpsipm);
        	cendfc =   (dpsip-dpsim) / (dpsip*dpsim);
        	cendfp =   dpsim / (dpsip*dpsipm);
		/**************************************************/


		/*evaluate properties at j*************************/
        //printf("Set gas at j=%d\n",j);
		// setGas(data,ydata,j,GAS);
		setGas(data,ydata,j);
		rho=gas->density();		//kg/m^3
        	Cpb=gas->cp_mass();       	//J/kg/K
        	Cvb=gas->cv_mass();       	//J/kg/K
		// gas->getNetProductionRates(wdot); //kmol/m^3

		if (data->rxn == 0){
			for (size_t k = 1; k <= nsp; k++)
			{
				wdot(k) = 0.0;
			}
		}else{
			gas->getNetProductionRates(wdot); //kmol/m^3
		}

        //DEBUG
        //if(T(j) > 0.0){
		//    data->gas->getNetProductionRates(wdot); //kmol/m^3
        //}else{
        //    for(size_t specII=0; specII < data->nsp; specII++){
        //        wdot[specII] = 0.0;
        //    }
        //}
		gas->getEnthalpy_RT(enthalpy);	//unitless
		gas->getCp_R(Cp);			//unitless
		area = calc_area(R(j),&m);	        //m^2

		/*evaluate properties at p*************************/
        //printf("get transport at j=%d\n",j);
		getTransport(data, ydata, j, &rhophalf,&lambdaphalf,YVphalf);
		areaphalf= calc_area(HALF*(R(j)+R(j+1)),&m);
		areaphalfsq= areaphalf*areaphalf;
		/**************************************************///

		/*Energy:*/
		/* ∂T/∂t = - ṁ(∂T/∂ψ)
		 * 	   + (1/cₚ)(∂/∂ψ)(λρA²∂T/∂ψ)
		 * 	   - (A/cₚ) ∑ YᵢVᵢcₚᵢ(∂T/∂ψ)
		 * 	   - (1/ρcₚ)∑ ώᵢhᵢ
		 * 	   + (1/ρcₚ)(∂P/∂t) */
		/*Notes:
		 * λ has units J/m/s/K.
		 * YᵢVᵢ has units kg/m^2/s.
		 * hᵢ has units J/kmol, so we must multiply the enthalpy
		 * defined above (getEnthalpy_RT) by T (K) and the gas constant
		 * (J/kmol/K) to get the right units.
		 * cₚᵢ has units J/kg/K, so we must multiply the specific heat
		 * defined above (getCp_R) by the gas constant (J/kmol/K) and
		 * divide by the molecular weight (kg/kmol) to get the right
		 * units.
		 * */

		//enthalpy formulation:
      if(j == 2){
		   sum=ZERO;
		   sum1=ZERO;
		   for (k = 1; k <=nsp; k++) {
		   	sum=sum+wdot(k)*enthalpy(k);
//		   	 sum1=sum1+(Cp(k)/data->gas->molecularWeight(k-1))
//		   	      *HALF*(YVphalf(k));
			sum1=sum1+(Cp(k)/gas->molecularWeight(k-1))
			     *HALF*(YVmhalf(k)+YVphalf(k));
		   }
		   sum=sum*Cantera::GasConstant*T(j);
		   sum1=sum1*Cantera::GasConstant;
//		    diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip))
//                                        /(dpsiav*Cpb) )
//		    	        -(sum1*area*(cendfp*T(j+1)
//		    	                    +cendfc*T(j)
//		    	                    +cendfm*T(j))/Cpb);
			diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip)
							-(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
							/(dpsiav*Cpb) )
						-(sum1*area*(cendfp*T(j+1)
									+cendfc*T(j)
									+cendfm*T(j-1))/Cpb);
		//    srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata)+Qrad(data,ydata,j))/(rho*Cpb);
			srcTerm = sum/(rho*Cpb);
		   //srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata))/(rho*Cpb);
		   //advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
		   advTerm   = data->Mdot*(T(j)-T(j-1))/dpsim;
		   Tdot(j)= data->dPdt/(rho*Cpb)
		   	      -advTerm
		   	      +diffTerm
		   	      -srcTerm;
           if(!data->dirichletInner) { //neumann inner
               Tdot(j - 1) = Tdot(j);
           }
      }else if(j == npts-1){
		   sum=ZERO;
		   sum1=ZERO;
		   for (k = 1; k <=nsp; k++) {
		   	sum=sum+wdot(k)*enthalpy(k);
               sum1=sum1+(Cp(k)/gas->molecularWeight(k-1))
		   	      *HALF*(YVmhalf(k)+YVphalf(k));
//			sum1=sum1+(Cp(k)/data->gas->molecularWeight(k-1))
//			     *HALF*(YVmhalf(k));
		   }
		   sum=sum*Cantera::GasConstant*T(j);
		   sum1=sum1*Cantera::GasConstant;
		   diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip)
			      -(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
			     /(dpsiav*Cpb) )
			  -(sum1*area*(cendfp*T(j+1)
			              +cendfc*T(j)
			              +cendfm*T(j-1))/Cpb);
//		    diffTerm  =((-(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
//		    	                         /(dpsiav*Cpb) )
//		    	  -(sum1*area*(cendfp*T(j)
//		    	              +cendfc*T(j)
//		    	              +cendfm*T(j-1))/Cpb);
		//    srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata)+Qrad(data,ydata,j))/(rho*Cpb);
			srcTerm = sum/(rho*Cpb);
		   //srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata))/(rho*Cpb);
		//    advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
			advTerm   = data->Mdot*(T(j)-T(j-1))/dpsim;
		   Tdot(j)= data->dPdt/(rho*Cpb)
		   	      -advTerm
		   	      +diffTerm
		   	      -srcTerm;
           if(!data->dirichletOuter){//nermann temperature outer BC
               Tdot(j+1) = Tdot(j);
           }
      }else{
		   sum=ZERO;
		   sum1=ZERO;
		   for (k = 1; k <=nsp; k++) {
		   	sum=sum+wdot(k)*enthalpy(k);
		   	sum1=sum1+(Cp(k)/gas->molecularWeight(k-1))
		   	     *HALF*(YVmhalf(k)+YVphalf(k));
		   }
		   sum=sum*Cantera::GasConstant*T(j);
		   sum1=sum1*Cantera::GasConstant;
		   diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip)
		   	      -(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
		   	     /(dpsiav*Cpb) )
		   	  -(sum1*area*(cendfp*T(j+1)
		   	              +cendfc*T(j)
		   	              +cendfm*T(j-1))/Cpb);
		//    srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata)+Qrad(data,ydata,j))/(rho*Cpb);
		   //srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata))/(rho*Cpb);
			srcTerm = sum/(rho*Cpb);
		//    advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
			advTerm   = data->Mdot*(T(j)-T(j-1))/dpsim;
		   Tdot(j)= data->dPdt/(rho*Cpb)
		   	      -advTerm
		   	      +diffTerm
		   	      -srcTerm;
      }

	//	//energy formulation:
	//	tranTerm = Tdot(j);
	//	sum=ZERO;
	//	sum1=ZERO;
	//	sum2=ZERO;
	//	sum3=ZERO;
	//	for (k = 1; k <=nsp; k++) {
	//		energy(k)=enthalpy(k)-ONE;
	//		sum=sum+wdot(k)*energy(k);
	//		sum1=sum1+(Cp(k)/data->gas->molecularWeight(k-1))*rho
	//		     *HALF*(YVmhalf(k)+YVphalf(k));
	//		sum2=sum2+(YVmhalf(k)/data->gas->molecularWeight(k-1));
	//		sum3=sum3+(YVphalf(k)/data->gas->molecularWeight(k-1));
	//	}
	//	sum=sum*Cantera::GasConstant*T(j);
	//	sum1=sum1*Cantera::GasConstant;
	//	diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip)
	//		      -(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
	//		     /(dpsiav*Cvb) )
	//		  -(sum1*area*(cendfp*T(j+1)
	//		              +cendfc*T(j)
	//		              +cendfm*T(j-1))/Cvb);
	//	srcTerm   = (sum-Qdot(&t,&R(j),&data->ignTime,&data->kernelSize,&data->maxQDot))/(rho*Cvb);
	//	advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
	//	advTerm = advTerm + (Cantera::GasConstant*T(j)*area/Cvb)*((sum3-sum2)/dpsiav);
	//	Tres(j)= tranTerm
	//		+advTerm
	//		-diffTerm
	//		+srcTerm;

		/*Species:*/
		/* ∂Yᵢ/∂t = - ṁ(∂Yᵢ/∂ψ)
		 * 	    - (∂/∂ψ)(AYᵢVᵢ)
		 * 	    + (ώᵢWᵢ/ρ)  */
      if(j == 2){
		   sum=ZERO;
		   for (k = 1; k <=nsp; k++) {
		   	if(k!=k_bath){
		   		// diffTerm  = (YVphalf(k)*areaphalf)/dpsiav;
				diffTerm  = (YVphalf(k)*areaphalf
					    -YVmhalf(k)*areamhalf)/dpsiav;
		   		srcTerm   = wdot(k)
		   	     		    *(gas->molecularWeight(k-1))/rho;
		   		// advTerm   = 0.0;
				advTerm   = (data->Mdot*(Y(j,k)-Y(j-1,k))/dpsim);
		   		Ydot(j,k)= -advTerm
		   			        -diffTerm
		   			        +srcTerm;

            //    Ydot(j-1,k) = Ydot(j,k);

               //DEBUG
               //if(k == data->k_e){
               //   printf("tranTerm= %0.5e\n diffTerm= %0.5e\n srcTerm= %0.5e\n advTerm= %0.5e\n\n",
               //           tranTerm,         diffTerm,         srcTerm,         advTerm);
               //}

		        	sum=sum+Ydot(j,k);
		   	}
		   }
		   Ydot(j,k_bath)=-sum;

          if (!data->dirichletInner){ //neumann inner BC
              for (k=1;k<=nsp;k++){
                  if (k!=k_bath){
                      Ydot(j-1,k) = Ydot(j,k);
                  }else{
                      Ydot(j-1,k_bath)=Ydot(j,k_bath);
                  }
              }
          }

		//    Ydot(j-1,k_bath)=Ydot(j,k_bath);

      }else if(j == npts-1){
		   sum=ZERO;
		   for (k = 1; k <=nsp; k++) {
		   	if(k!=k_bath){
//		   		 diffTerm  = (-YVmhalf(k)*areamhalf)/dpsiav;
				diffTerm  = (YVphalf(k)*areaphalf
					    -YVmhalf(k)*areamhalf)/dpsiav;
		   		srcTerm   = wdot(k)
		   	     		    *(gas->molecularWeight(k-1))/rho;
		   		// advTerm   = (mdotIn*(Y(j,k)-Y(j-1,k))/dpsim);
				advTerm   = (data->Mdot*(Y(j,k)-Y(j-1,k))/dpsim);
		   		Ydot(j,k)= -advTerm
		   			        -diffTerm
		   			        +srcTerm;

               //DEBUG
               //if(k == data->k_e){
               //   printf("tranTerm= %0.5e\n diffTerm= %0.5e\n srcTerm= %0.5e\n advTerm= %0.5e\n\n",
               //           tranTerm,         diffTerm,         srcTerm,         advTerm);
               //}

		        	sum=sum+Ydot(j,k);
		   	}
		   }
		   Ydot(j,k_bath)=-sum;

           if(!data->dirichletOuter){ //neumann temperature outer BC
               for (k=1;k<=nsp;k++){
                   if (k!=k_bath){
                       Ydot(j+1,k)=Ydot(j,k);
                   }else{
                       Ydot(j+1,k_bath)=Ydot(j,k_bath);
                   }
               }
           }
      }else{
		   sum=ZERO;
		   for (k = 1; k <=nsp; k++) {
		   	if(k!=k_bath){
		   		diffTerm  = (YVphalf(k)*areaphalf
		   			    -YVmhalf(k)*areamhalf)/dpsiav;
		   		srcTerm   = wdot(k)
		   	     		    *(gas->molecularWeight(k-1))/rho;
		   		// advTerm   = (mdotIn*(Y(j,k)-Y(j-1,k))/dpsim);
				advTerm   = (data->Mdot*(Y(j,k)-Y(j-1,k))/dpsim);
		   		Ydot(j,k)= -advTerm
		   			        -diffTerm
		   			        +srcTerm;
               //DEBUG
               //if(k == data->k_e){
               //   printf("tranTerm= %0.5e\n diffTerm= %0.5e\n srcTerm= %0.5e\n advTerm= %0.5e\n\n",
               //           tranTerm,         diffTerm,         srcTerm,         advTerm);
               //}

		        	sum=sum+Ydot(j,k);
		   	}
		   }
		   Ydot(j,k_bath)=-sum;

      }

      //TODO Account for constant volume case
		/*Pressure:*/
		Pdot(j) = data->dPdt;

		/*Assign values evaluated at p and phalf to m
		 * and mhalf to save some cpu cost:****************/
		areamhalf=areaphalf;
		areamhalfsq=areaphalfsq;
		aream=area;
		rhom=rho;
		rhomhalf=rhophalf;
		lambdamhalf=lambdaphalf;
		for (k = 1; k <=nsp; k++) {
			YVmhalf(k)=YVphalf(k);
		}
		/**************************************************/

	}
	/*******************************************************************///


	/*Fill up res with right side (wall) boundary conditions:***********/
	/*We impose zero fluxes at the wall:*/
    //printf("set gas at j=N\n");
	// setGas(data,ydata,npts,GAS);
	setGas(data,ydata,npts);
	rho=gas->density();
	area = calc_area(R(npts),&m);

	/*Mass:*/
	dpsim=(psi(npts)-psi(npts-1))*mass;

	/*Energy:*/
	if(data->dirichletOuter){
		Tdot(npts)=0.0;
        for (k = 1; k <= nsp; k++){
            Ydot(1, k) = 0.0;
        }
	}

	/*Pressure:*/
	if(data->constantPressure){
		Pdot(npts) = data->dPdt;
	}
	else{
      //TODO Account for this and pressure across all points for isochoric case
		//Pdot(npts)=R(npts)-data->domainLength;
		//Pres(npts)=Rdot(npts);
	}


	/*******************************************************************/

	//for (j = 1; j <=npts; j++) {
	//	//for (k = 1; k <=nsp; k++) {
	//	//	Yres(j,k)=Ydot(j,k);
	//	//}
	//	//Tres(j)=Tdot(j);
	//}

	return(0);
}



int funNew(double t,
        N_Vector y,
        N_Vector ydot,
        void *user_data){

    /*Declare and fetch nvectors and user data:*/

    double *ydata, *ydotdata, *psidata, *innerMassFractionsData;

    UserData data;
    data = (UserData)user_data;
    const size_t npts=data->npts;
    const size_t nlpts = data->nlpts;
    const size_t nsp=data->nsp;
    const size_t k_bath = data->k_bath;

//	getR(y, data);

    ydata	= N_VGetArrayPointer_OpenMP(y);
    ydotdata = N_VGetArrayPointer_OpenMP(ydot);
    if(data->adaptiveGrid==1){
        psidata = data->grid->x;
    }else{
        psidata = data->uniformGrid;
    }

//    getR(ydata,data);
    getRNew(ydata,data);// update R array

    innerMassFractionsData =  data->innerMassFractions;

    /* Grid stencil:*/

    /*-------|---------*---------|---------*---------|-------*/
    /*-------|---------*---------|---------*---------|-------*/
    /*-------|---------*---------|---------*---------|-------*/
    /*-------m-------mhalf-------j-------phalf-------p-------*/
    /*-------|---------*---------|---------*---------|-------*/
    /*-------|---------*---------|---------*---------|-------*/
    /*-------|<=======dxm=======>|<=======dxp=======>|-------*/
    /*-------|---------*<======dxav=======>*---------|-------*/
    /*-------|<================dxpm=================>|-------*/

    /* Various variables defined for book-keeping and storing previously
     * calculated values:
     * rho		: densities at points  m, mhalf, j, p, and phalf.
     * area		: the matric at points m, mhalf, j, p, and phalf.
     * m 		: exponent that determines geometry;
     * lambda	: thermal conductivities at mhalf and phalf.
     * mdot		: mass flow rate at m, j, and p.
     * X		: mole fractions at j and p.
     * YV		: diffusion fluxes at mhalf and phalf.
     * Tgrad	: temperature gradient at mhalf and phalf.
     * Tav		: average temperature between two points.
     * Pav		: average pressure between two points.
     * Yav		: average mass fractions between two points.
     * Xgradhalf	: mole fraction gradient at j.
     * Cpb		: mass based bulk specific heat.
     * tranTerm	: transient terms.
     * advTerm	: advection terms.
     * diffTerm	: diffusion terms.
     * srcTerm	: source terms.
     */



//    dpsip=dpsiav=dpsipm=dpsim=dpsimm=ONE;
    double mass, mdotIn;
    double sum, sum1;
    double Mdot,dropletMass;
    double interfaceYArray[nsp] ;

    for (size_t k = 0; k < nsp; k++) {
        interfaceYArray[k] = data->interfaceGasCellArr[2 + k];
    }

//    size_t j,k;
    int m;
    m=data->metric;					//Unitless
    mass=data->mass;				//Units: kg
    mdotIn=data->mdot*calc_area(R(npts),&m);	//Units: kg/s
//    Mdot = data->Mdot; //units: kg/s
    dropletMass = data->dropletMass; //unit:kg

    /* GAS PHASE PART GOVERNING EQUATIONS CODE SECTION */

//	/*evaluate properties at j=1*************************/
    //printf("Set gas at j=1\n");
    // setGas(data,ydata,1,GAS);
//    setGas(data,ydata,1);
//    gas->setState_TPY(data->interfaceGasCellArr[1],data->interfaceGasCellArr[data->nvar],&data->interfaceGasCellArr[2]);
//    rhom=gas->density();
//    Cpb=gas->cp_mass();       	//J/kg/K
//    Cvb=gas->cv_mass();       	//J/kg/K
////    aream= calc_area(gasR(1),&m);

    /*******************************************************************/
    /*Calculate values at j=2's m and mhalf*****************************/

    //printf("Get transport at j=1\n");
//    getTransport(data, ydata, 1, &rhomhalf,&lambdamhalf,YVmhalf);

    // flux variables of size : npts
    double rhohalf[npts], lambdahalf[npts],YVhalf[npts][nsp],areahalf[npts],areahalfsq[npts];

    // handle the first element of each flux variables array
    getInterfaceTransport(data,ydata,&rhohalf[0],&lambdahalf[0],YVhalf[0]) ;
    areahalf[0] = calc_area(HALF*(gasR(1)+gasR(2)),&m);
    areahalfsq[0] = areahalf[0]*areahalf[0];

#pragma omp parallel for default(none) shared(rhohalf,lambdahalf,YVhalf,areahalf,nlpts, npts,areahalfsq,data,ydata,m)
    for (size_t j = nlpts+1;j<nlpts+npts;j++){
        size_t gasII = j- nlpts +1 ; // gasII begin with 2
        areahalf[gasII-1] = calc_area(HALF*(gasR(gasII)+ gasR(gasII+1)),&m);
        areahalfsq[gasII-1] = areahalf[gasII-1] * areahalf[gasII-1];
//        getTransport(data,ydata,j+1,&rhohalf[gasII-1],&lambdahalf[gasII-1],YVhalf[gasII-1]) ;
        getTransportNew(data,ydata,j,&rhohalf[gasII-1],&lambdahalf[gasII-1],YVhalf[gasII-1]) ;
    }
//    getInterfaceTransport(data,ydata,&rhomhalf,&lambdamhalf,YVmhalf);
//    areamhalf= calc_area(HALF*(gasR(1)+gasR(2)),&m);
//    areamhalfsq= areamhalf*areamhalf;

    /*Calculate the droplet vaporization rate: kg/s */
//    data->Mdot = YVmhalf(data->dropII) * areamhalf / (1 - Y(1, data->dropII)) ;
    Mdot = YVhalf[0][data->dropII-1] * areahalf[0]/ (1-interfaceYArray[data->dropII-1]); //units: kg/s
//    Mdot = YVmhalf(data->dropII) * areamhalf / (1- interfaceYArray[data->dropII-1]); //units: kg/s
    /*******************************************************************/

    /*Fill up res with left side (center) boundary conditions:**********/
    /*We impose zero fluxes at the center:*/

    /*Mass:*/
    //Rdot(1) = 0.0;

    /*Energy:*/
//
//    if (data->dirichletInner){
//        Tdot(1) = 0.0;
//        // TEST:try to fix the boundary mass fraction
//        for (k = 1; k <= nsp; k++){
//            Ydot(1, k) = 0.0;
//        }
//    }
    //else{
    //	Tres(1)=T(2)-T(1);
    //	//Tres(1)=Tdot(1) - (Pdot(1)/(rhom*Cpb))
    //	//	+(double)(data->metric+1)*(rhomhalf*lambdamhalf*areamhalfsq*(T(2)-T(1))/psi(2)-psi(1));
    //}

    /*Species:*/
//    sum=ZERO;
//    for (k = 1; k <=nsp; k++) {
//        if(k!=k_bath){
//            if(fabs(mdotIn)>1e-14){
//                //TODO Account for mdotin
//                //	Ydot(1,k)=innerMassFractionsData[k-1]-
//                //  Y(1,k)-
//                //  (YVmhalf(k)*areamhalf)/mdotIn;
//            }
//        }
//    }
    //Yres(1,k_bath)=ONE-sum-Y(1,k_bath);


    /*Pressure:*/
    //TODO implement constant volume case
    Pdot(1+nlpts)=data->dPdt;

#pragma  omp parallel default(none) shared(nsp,nlpts,npts,k_bath,data,ydata,ydotdata,psidata,\
    m,mass,dropletMass,Mdot,mdotIn,rhohalf,lambdahalf,YVhalf,areahalf,areahalfsq ,\
    interfaceYArray)
    {
        double rho,Cpb, Cvb,  Cp[nsp],   wdot[nsp],   enthalpy[nsp],
                diffTerm,    srcTerm, advTerm,
                area;
        double sum, sum1;
        /*Aliases for difference coefficients:*/
        double cendfm, cendfc, cendfp;
        cendfm=cendfc=cendfp=ZERO;
        /*Aliases for various grid spacings:*/
        double dpsip, dpsiav, dpsipm, dpsim, dpsimm;
        dpsip=dpsiav=dpsipm=dpsim=dpsimm=ONE;

#pragma omp for
        /*Fill up res with governing equations at inner points:*************/
        for (size_t j = 1+nlpts; j < nlpts+npts; j++) {

            size_t gasII = j-nlpts+1; // begins with gasII = 2, for gasR() ans psi() index ONLY
            size_t fluxArrayIndex = j - nlpts; // begins with fluxArrayIndex = 1, for rhohalf,lambdahalf,YVhalf,areahalf,areahalfsq index


            /*evaluate various mesh differences*///
//        dpsip =        (psi(j+1) - psi(j)  )*mass;
//        dpsim =        (psi(j)   - psi(j-1))*mass;
//        dpsiav =  HALF*(psi(j+1) - psi(j-1))*mass;
//        dpsipm =       (psi(j+1) - psi(j-1))*mass;
            dpsip =        (psi(gasII+1) - psi(gasII)  )*mass;
            dpsim =        (psi(gasII)   - psi(gasII-1))*mass;
            dpsiav =  HALF*(psi(gasII+1) - psi(gasII-1))*mass;
            dpsipm =       (psi(gasII+1) - psi(gasII-1))*mass;


            /***********************************///

            /*evaluate various central difference coefficients*/
            cendfm = - dpsip / (dpsim*dpsipm);
            cendfc =   (dpsip-dpsim) / (dpsip*dpsim);
            cendfp =   dpsim / (dpsip*dpsipm);
            /**************************************************/


            /*evaluate properties at j*************************/
            //printf("Set gas at j=%d\n",j);
            // setGas(data,ydata,j,GAS);
            setGas(data,ydata,j);
            rho=gas->density();		//kg/m^3
            Cpb=gas->cp_mass();       	//J/kg/K
            Cvb=gas->cv_mass();       	//J/kg/K
            // data->gas->getNetProductionRates(wdot); //kmol/m^3

            if (data->rxn == 0){
                for (size_t k = 1; k <= nsp; k++)
                {
                    wdot(k) = 0.0;
                }
            }else{
                gas->getNetProductionRates(wdot); //kmol/m^3
            }

            //DEBUG
            //if(T(j) > 0.0){
            //    data->gas->getNetProductionRates(wdot); //kmol/m^3
            //}else{
            //    for(size_t specII=0; specII < data->nsp; specII++){
            //        wdot[specII] = 0.0;
            //    }
            //}
            gas->getEnthalpy_RT(enthalpy);	//unitless
            gas->getCp_R(Cp);			//unitless
//        area = calc_area(R(j),&m);	        //m^2
            area = calc_area(gasR(gasII),&m);
            /*evaluate properties at p*************************/
            //printf("get transport at j=%d\n",j);
//            getTransport(data, ydata, j, &rhophalf,&lambdaphalf,YVphalf);
//        areaphalf= calc_area(HALF*(R(j)+R(j+1)),&m);
//            areaphalf = calc_area(HALF*(gasR(gasII)+ gasR(gasII+1)),&m) ;
//            areaphalfsq= areaphalf*areaphalf;
            /**************************************************///

            /*Energy:*/
            /* ∂T/∂t = - ṁ(∂T/∂ψ)
             * 	   + (1/cₚ)(∂/∂ψ)(λρA²∂T/∂ψ)
             * 	   - (A/cₚ) ∑ YᵢVᵢcₚᵢ(∂T/∂ψ)
             * 	   - (1/ρcₚ)∑ ώᵢhᵢ
             * 	   + (1/ρcₚ)(∂P/∂t) */
            /*Notes:
             * λ has units J/m/s/K.
             * YᵢVᵢ has units kg/m^2/s.
             * hᵢ has units J/kmol, so we must multiply the enthalpy
             * defined above (getEnthalpy_RT) by T (K) and the gas constant
             * (J/kmol/K) to get the right units.
             * cₚᵢ has units J/kg/K, so we must multiply the specific heat
             * defined above (getCp_R) by the gas constant (J/kmol/K) and
             * divide by the molecular weight (kg/kmol) to get the right
             * units.
             * */

            //enthalpy formulation:
            if (j == nlpts + 1) {
                sum = ZERO;
                sum1 = ZERO;
                for (size_t k = 1; k <= nsp; k++) {
                    sum = sum + wdot(k) * enthalpy(k);
//		   	 sum1=sum1+(Cp(k)/data->gas->molecularWeight(k-1))
//		   	      *HALF*(YVphalf(k));
//                    sum1=sum1+(Cp(k)/gas->molecularWeight(k-1))
//                              *HALF*(YVmhalf(k)+YVphalf(k));
                    sum1 = sum1 + (Cp(k) / gas->molecularWeight(k - 1)) *
                                  HALF * (YVhalf[fluxArrayIndex - 1][k - 1] + YVhalf[fluxArrayIndex][k - 1]);
                }
                sum = sum * Cantera::GasConstant * T(j);
                sum1 = sum1 * Cantera::GasConstant;
//		    diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip))
//                                        /(dpsiav*Cpb) )
//		    	        -(sum1*area*(cendfp*T(j+1)
//		    	                    +cendfc*T(j)
//		    	                    +cendfm*T(j))/Cpb);
//                diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip)
//                              -(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-data->interfaceGasCellArr[1])/dpsim) )
//                            /(dpsiav*Cpb) )
//                           -(sum1*area*(cendfp*T(j+1)
//                                        +cendfc*T(j)
//                                        +cendfm*data->interfaceGasCellArr[1])/Cpb);
                diffTerm = (((rhohalf[fluxArrayIndex] * areahalfsq[fluxArrayIndex] * lambdahalf[fluxArrayIndex] *
                              (T(j + 1) - T(j)) / dpsip)
                             - (rhohalf[fluxArrayIndex - 1] * areahalfsq[fluxArrayIndex - 1] *
                                lambdahalf[fluxArrayIndex - 1] * (T(j) - data->interfaceGasCellArr[1]) / dpsim))
                            / (dpsiav * Cpb))
                           - (sum1 * area * (cendfp * T(j + 1)
                                             + cendfc * T(j)
                                             + cendfm * data->interfaceGasCellArr[1]) / Cpb);
                //    srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata)+Qrad(data,ydata,j))/(rho*Cpb);
                srcTerm = sum / (rho * Cpb);
                //srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata))/(rho*Cpb);
                //advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
                advTerm = Mdot * (T(j) - data->interfaceGasCellArr[1]) / dpsim;
                Tdot(j) = data->dPdt / (rho * Cpb)
                          - advTerm
                          + diffTerm
                          - srcTerm;
//            if(!data->dirichletInner) { //neumann inner
//                Tdot(j - 1) = Tdot(j);
//            }
            } else if (j == (nlpts + npts - 1)) {
                sum = ZERO;
                sum1 = ZERO;
                for (size_t k = 1; k <= nsp; k++) {
                    sum = sum + wdot(k) * enthalpy(k);
                    sum1 = sum1 + (Cp(k) / gas->molecularWeight(k - 1)) *
                                  HALF * (YVhalf[fluxArrayIndex - 1][k - 1] + YVhalf[fluxArrayIndex][k - 1]);
//			sum1=sum1+(Cp(k)/data->gas->molecularWeight(k-1))
//			     *HALF*(YVmhalf(k));
                }
                sum = sum * Cantera::GasConstant * T(j);
                sum1 = sum1 * Cantera::GasConstant;
                diffTerm = (((rhohalf[fluxArrayIndex] * areahalfsq[fluxArrayIndex] * lambdahalf[fluxArrayIndex] *
                              (T(j + 1) - T(j)) / dpsip)
                             - (rhohalf[fluxArrayIndex - 1] * areahalfsq[fluxArrayIndex - 1] *
                                lambdahalf[fluxArrayIndex - 1] * (T(j) - T(j - 1)) / dpsim))
                            / (dpsiav * Cpb))
                           - (sum1 * area * (cendfp * T(j + 1)
                                             + cendfc * T(j)
                                             + cendfm * T(j - 1)) / Cpb);
//		    diffTerm  =((-(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
//		    	                         /(dpsiav*Cpb) )
//		    	  -(sum1*area*(cendfp*T(j)
//		    	              +cendfc*T(j)
//		    	              +cendfm*T(j-1))/Cpb);
                //    srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata)+Qrad(data,ydata,j))/(rho*Cpb);
                srcTerm = sum / (rho * Cpb);
                //srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata))/(rho*Cpb);
                //    advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
                advTerm = Mdot * (T(j) - T(j - 1)) / dpsim;
                Tdot(j) = data->dPdt / (rho * Cpb)
                          - advTerm
                          + diffTerm
                          - srcTerm;
                if (!data->dirichletOuter) {//Neumann temperature outer BC
                    Tdot(j + 1) = Tdot(j);
                }
            } else {
                sum = ZERO;
                sum1 = ZERO;
                for (size_t k = 1; k <= nsp; k++) {
                    sum = sum + wdot(k) * enthalpy(k);
                    sum1 = sum1 + (Cp(k) / gas->molecularWeight(k - 1))
                                  * HALF * (YVhalf[fluxArrayIndex - 1][k - 1] + YVhalf[fluxArrayIndex][k - 1]);
                }
                sum = sum * Cantera::GasConstant * T(j);
                sum1 = sum1 * Cantera::GasConstant;
                diffTerm = (((rhohalf[fluxArrayIndex] * areahalfsq[fluxArrayIndex] * lambdahalf[fluxArrayIndex] *
                              (T(j + 1) - T(j)) / dpsip)
                             - (rhohalf[fluxArrayIndex - 1] * areahalfsq[fluxArrayIndex - 1] *
                                lambdahalf[fluxArrayIndex - 1] * (T(j) - T(j - 1)) / dpsim))
                            / (dpsiav * Cpb))
                           - (sum1 * area * (cendfp * T(j + 1)
                                             + cendfc * T(j)
                                             + cendfm * T(j - 1)) / Cpb);
                //    srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata)+Qrad(data,ydata,j))/(rho*Cpb);
                //srcTerm   = (sum-Qdot(&t,&R(j),data,j,ydata))/(rho*Cpb);
                srcTerm = sum / (rho * Cpb);
                //    advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
                advTerm = Mdot * (T(j) - T(j - 1)) / dpsim;
                Tdot(j) = data->dPdt / (rho * Cpb)
                          - advTerm
                          + diffTerm
                          - srcTerm;
            }

            //	//energy formulation:
            //	tranTerm = Tdot(j);
            //	sum=ZERO;
            //	sum1=ZERO;
            //	sum2=ZERO;
            //	sum3=ZERO;
            //	for (k = 1; k <=nsp; k++) {
            //		energy(k)=enthalpy(k)-ONE;
            //		sum=sum+wdot(k)*energy(k);
            //		sum1=sum1+(Cp(k)/data->gas->molecularWeight(k-1))*rho
            //		     *HALF*(YVmhalf(k)+YVphalf(k));
            //		sum2=sum2+(YVmhalf(k)/data->gas->molecularWeight(k-1));
            //		sum3=sum3+(YVphalf(k)/data->gas->molecularWeight(k-1));
            //	}
            //	sum=sum*Cantera::GasConstant*T(j);
            //	sum1=sum1*Cantera::GasConstant;
            //	diffTerm  =(( (rhophalf*areaphalfsq*lambdaphalf*(T(j+1)-T(j))/dpsip)
            //		      -(rhomhalf*areamhalfsq*lambdamhalf*(T(j)-T(j-1))/dpsim) )
            //		     /(dpsiav*Cvb) )
            //		  -(sum1*area*(cendfp*T(j+1)
            //		              +cendfc*T(j)
            //		              +cendfm*T(j-1))/Cvb);
            //	srcTerm   = (sum-Qdot(&t,&R(j),&data->ignTime,&data->kernelSize,&data->maxQDot))/(rho*Cvb);
            //	advTerm   = (mdotIn*(T(j)-T(j-1))/dpsim);
            //	advTerm = advTerm + (Cantera::GasConstant*T(j)*area/Cvb)*((sum3-sum2)/dpsiav);
            //	Tres(j)= tranTerm
            //		+advTerm
            //		-diffTerm
            //		+srcTerm;

            /*Species:*/
            /* ∂Yᵢ/∂t = - ṁ(∂Yᵢ/∂ψ)
             * 	    - (∂/∂ψ)(AYᵢVᵢ)
             * 	    + (ώᵢWᵢ/ρ)  */
            if (j == nlpts + 1) {
                sum = ZERO;
                for (size_t k = 1; k <= nsp; k++) {
                    if (k != k_bath) {
                        // diffTerm  = (YVphalf(k)*areaphalf)/dpsiav;
                        diffTerm = (YVhalf[fluxArrayIndex][k - 1] * areahalf[fluxArrayIndex]
                                    - YVhalf[fluxArrayIndex - 1][k - 1] * areahalf[fluxArrayIndex - 1]) / dpsiav;
                        srcTerm = wdot(k)
                                  * (gas->molecularWeight(k - 1)) / rho;
                        // advTerm   = 0.0;
//                    advTerm   = (data->Mdot*(Y(j,k)-Y(j-1,k))/dpsim);
                        advTerm = (Mdot * (Y(j, k) - interfaceYArray[k - 1])) / dpsim;
                        Ydot(j, k) = -advTerm
                                     - diffTerm
                                     + srcTerm;
                        //    Ydot(j-1,k) = Ydot(j,k);

                        //DEBUG
                        //if(k == data->k_e){
                        //   printf("tranTerm= %0.5e\n diffTerm= %0.5e\n srcTerm= %0.5e\n advTerm= %0.5e\n\n",
                        //           tranTerm,         diffTerm,         srcTerm,         advTerm);
                        //}

                        sum = sum + Ydot(j, k);
                    }
                }
                Ydot(j, k_bath) = -sum;

//                if (!data->dirichletInner){ //neumann inner BC
//                    for (size_t k=1;k<=nsp;k++){
//                        if (k!=k_bath){
//                            Ydot(j-1,k) = Ydot(j,k);
//                        }else{
//                            Ydot(j-1,k_bath)=Ydot(j,k_bath);
//                        }
//                    }
//                }

                //    Ydot(j-1,k_bath)=Ydot(j,k_bath);

            } else if (j == nlpts + npts - 1) {
                sum = ZERO;
                for (size_t k = 1; k <= nsp; k++) {
                    if (k != k_bath) {
//		   		 diffTerm  = (-YVmhalf(k)*areamhalf)/dpsiav;
                        diffTerm = (YVhalf[fluxArrayIndex][k - 1] * areahalf[fluxArrayIndex]
                                    - YVhalf[fluxArrayIndex - 1][k - 1] * areahalf[fluxArrayIndex - 1]) / dpsiav;
                        srcTerm = wdot(k)
                                  * (gas->molecularWeight(k - 1)) / rho;
                        // advTerm   = (mdotIn*(Y(j,k)-Y(j-1,k))/dpsim);
                        advTerm = (Mdot * (Y(j, k) - Y(j - 1, k)) / dpsim);
                        Ydot(j, k) = -advTerm
                                     - diffTerm
                                     + srcTerm;
                        //DEBUG
                        //if(k == data->k_e){
                        //   printf("tranTerm= %0.5e\n diffTerm= %0.5e\n srcTerm= %0.5e\n advTerm= %0.5e\n\n",
                        //           tranTerm,         diffTerm,         srcTerm,         advTerm);
                        //}

                        sum = sum + Ydot(j, k);
                    }
                }
                Ydot(j, k_bath) = -sum;

                if (!data->dirichletOuter) { //neumann temperature outer BC
                    for (size_t k = 1; k <= nsp; k++) {
                        Ydot(j + 1, k) = Ydot(j, k);
                    }
                }
            } else {
                sum = ZERO;
                for (size_t k = 1; k <= nsp; k++) {
                    if (k != k_bath) {
                        diffTerm = (YVhalf[fluxArrayIndex][k - 1] * areahalf[fluxArrayIndex]
                                    - YVhalf[fluxArrayIndex - 1][k - 1] * areahalf[fluxArrayIndex - 1]) / dpsiav;
                        srcTerm = wdot(k)
                                  * (gas->molecularWeight(k - 1)) / rho;
                        // advTerm   = (mdotIn*(Y(j,k)-Y(j-1,k))/dpsim);
                        advTerm = (Mdot * (Y(j, k) - Y(j - 1, k)) / dpsim);
                        Ydot(j, k) = -advTerm
                                     - diffTerm
                                     + srcTerm;
                        //DEBUG
                        //if(k == data->k_e){
                        //   printf("tranTerm= %0.5e\n diffTerm= %0.5e\n srcTerm= %0.5e\n advTerm= %0.5e\n\n",
                        //           tranTerm,         diffTerm,         srcTerm,         advTerm);
                        //}

                        sum = sum + Ydot(j, k);
                    }
                }
                Ydot(j, k_bath) = -sum;
            }

            //TODO Account for constant volume case
            /*Pressure:*/
            Pdot(j) = data->dPdt;

            /*Assign values evaluated at p and phalf to m
             * and mhalf to save some cpu cost:****************/
//            areamhalf=areaphalf;
//            areamhalfsq=areaphalfsq;
//            aream=area;
//            rhom=rho;
//            rhomhalf=rhophalf;
//            lambdamhalf=lambdaphalf;
//            for (k = 1; k <=nsp; k++) {
//                YVmhalf(k)=YVphalf(k);
//            }
            /**************************************************/

        }

    }

    /*******************************************************************///


    /*Fill up res with right side (wall) boundary conditions:***********/
    /*We impose zero fluxes at the wall:*/
    //printf("set gas at j=N\n");
    // setGas(data,ydata,npts,GAS);
//    setGas(data,ydata,nlpts+npts);
//    rho=gas->density();
//    area = calc_area(gasR(npts+1),&m);
//
//    /*Mass:*/
//    dpsim=(psi(npts+1)-psi(npts))*mass;

    /*Energy:*/
    if(data->dirichletOuter){
        Tdot(npts+nlpts)=ZERO;
        for (size_t k = 1; k <= nsp; k++){
            Ydot(npts+nlpts, k) = ZERO;
        }
    }

    /*Pressure:*/
    if(data->constantPressure){
        Pdot(npts+nlpts) = data->dPdt;
    }
    else{
        //TODO Account for this and pressure across all points for isochoric case
        //Pdot(npts)=R(npts)-data->domainLength;
        //Pres(npts)=Rdot(npts);
    }


    /*******************************************************************/

    //for (j = 1; j <=npts; j++) {
    //	//for (k = 1; k <=nsp; k++) {
    //	//	Yres(j,k)=Ydot(j,k);
    //	//}
    //	//Tres(j)=Tdot(j);
    //}

    // flux variables of nlpts-1
    double rhoLiquidhalf[nlpts-1],lambdaLiquidhalf[nlpts-1],areaLiquidhalf[nlpts-1],areaLiquidhalfsq[nlpts-1];

#pragma omp parallel for default(none) shared(rhoLiquidhalf,lambdaLiquidhalf,areaLiquidhalf,areaLiquidhalfsq,nlpts,data,ydata,m)
    for (size_t j = 1; j <= nlpts - 1; j++) {
        areaLiquidhalf[j - 1] = calc_area(HALF * (liquidR(j) + liquidR(j + 1)), &m);
        areaLiquidhalfsq[j - 1] = areaLiquidhalf[j - 1] * areaLiquidhalf[j - 1];
        singleLiquidFuel->setLiquidState(HALF * (T(j) + T(j + 1)), P(j));
        lambdaLiquidhalf[j - 1] = singleLiquidFuel->calculateThermalConductivity();
        rhoLiquidhalf[j - 1] = singleLiquidFuel->calculateDensity();
    }

    /***** LIQUID PHASE PART GOVERNING EQUATIONS CODE SECTION *****/
    // -- Precompute initial states for j = 1 and j = 2
//    singleLiquidFuel->setLiquidState(T(1),P(1));
//    rhom = singleLiquidFuel->calculateDensity();
//    lambdamhalf = singleLiquidFuel->calculateThermalConductivity();
//
//    singleLiquidFuel->setLiquidState(T(2),P(2));
//    rho = singleLiquidFuel->calculateDensity();
//
//    // j = 1+1/2
//    singleLiquidFuel->setLiquidState(HALF*(T(1)+ T(2)),HALF* (P(1)+P(2)));
//    rhomhalf = singleLiquidFuel->calculateDensity();
//
//    areamhalf = calc_area(HALF*(liquidR(1)+ liquidR(2)),&m);
//    areamhalfsq = areamhalf * areamhalf;

#pragma omp parallel default(none) shared(rhoLiquidhalf,lambdaLiquidhalf,areaLiquidhalf,areaLiquidhalfsq,nlpts,data,ydata,ydotdata,psidata,m,\
    dropletMass,Mdot,nsp)
    {
        double rho,Cpb, diffTerm, advTerm , area;
        /*Aliases for various grid spacings:*/
    //    double dpsip, dpsiav, dpsipm, dpsim, dpsimm;
        double detap, detaav, detapm, detam, detamm;
        detap=detaav=detapm=detam=detamm=ONE;

#pragma omp for
    for(size_t j =2;j<= (nlpts-1);j++){ // handle inner grid points
        detam = std::abs(eta(j)- eta(j-1));
        detap = std::abs(eta(j+1) - eta(j));
        detaav = HALF * std::abs(eta(j+1) - eta(j-1));
        detapm = std::abs(eta(j+1) - eta(j-1));

        // Zero out Pdot and Ydot at j
        Pdot(j) = data->dPdt;
        for (size_t k = 1; k <= nsp; ++k) {
            Ydot(j, k) = 0.0;
        }

        // set liquid phase and fetch corresponding parameters
        singleLiquidFuel->setLiquidState(T(j),P(j));
        rho = singleLiquidFuel->calculateDensity();
        Cpb = singleLiquidFuel->calculateSpecificHeatCapacity() ;
        area= calc_area(liquidR(j),&m);


        if (j==2) {
            Pdot(j-1) = Pdot(j);
            for (size_t k = 1;k <= nsp;k++){
                Ydot(j-1,k) = Ydot(j,k);
            }
//            advTerm = ZERO; // neumann B.C. at droplet center
            diffTerm = 1/ (Cpb * dropletMass * dropletMass) *
                       (rhoLiquidhalf[j-1] * areaLiquidhalfsq[j-1] * lambdaLiquidhalf[j-1] * (T(j + 1) - T(j)) / detap) / detaav;
            Tdot(j) = diffTerm ;
            Tdot(j-1) = Tdot(j) ;
        }else{
            advTerm = eta(j) * (-Mdot) / dropletMass * (T(j)-T(j-1)) / detam;
            diffTerm = 1.0 / (Cpb * dropletMass*dropletMass) *
                       (rhoLiquidhalf[j-1] * areaLiquidhalfsq[j-1] * lambdaLiquidhalf[j-1] * (T(j + 1) - T(j)) / detap -
                               rhoLiquidhalf[j-2] * areaLiquidhalfsq[j-2] * lambdaLiquidhalf[j-2] * (T(j) - T(j - 1)) / detam) / detaav;
            Tdot(j) = advTerm + diffTerm;
        }
}
//        //re-use some parameters
//        lambdamhalf = lambdaphalf;
//        areamhalf = areaphalf;
//        areamhalfsq = areaphalfsq;
//        rhom = rho;
//        rho = rhop;
//        rhomhalf = rhophalf;
    }

    {    // last internal point handling
        double detam,detap,detaav,
            rho,Cpb,lambdaphalf,
            areaphalf,areaphalfsq,areamhalf,areamhalfsq,
            advTerm,diffTerm,rhomhalf,lambdamhalf;

        detam = std::abs(eta(nlpts)- eta(nlpts-1));
        detap = std::abs(eta(nlpts+1) - eta(nlpts));
        detaav = HALF * std::abs(eta(nlpts+1) - eta(nlpts-1));

        singleLiquidFuel->setLiquidState(T(nlpts),P(nlpts));
        rho = singleLiquidFuel->calculateDensity();
        Cpb = singleLiquidFuel->calculateSpecificHeatCapacity();
        lambdaphalf = singleLiquidFuel->calculateThermalConductivity();

        areaphalf = calc_area(HALF*(liquidR(nlpts)+liquidR(nlpts+1)),&m) ;
        areaphalfsq = areaphalf * areaphalf ;

        areamhalf = calc_area(HALF*(liquidR(nlpts)+liquidR(nlpts-1)),&m) ;
        areamhalfsq = areamhalf * areamhalf ;
        rhomhalf = rhoLiquidhalf[nlpts-2] ;
        lambdamhalf = lambdaLiquidhalf[nlpts-2];

        Pdot(nlpts) = ZERO ;
        for (size_t k = 1;k <= nsp;k++){
            Ydot(nlpts,k) = ZERO;
        }
        advTerm = eta(nlpts) * (-Mdot) / dropletMass * (T(nlpts)-T(nlpts-1)) * detam;
        diffTerm = 1.0 / (Cpb * dropletMass*dropletMass) *
                   (rho * areaphalfsq * lambdaphalf * ( data->interfaceLiquidCellArr[1] - T(nlpts)) / detap -
                    rhomhalf * areamhalfsq * lambdamhalf * (T(nlpts) - T(nlpts - 1)) / detam) / detaav;
        Tdot(nlpts) = advTerm + diffTerm;
    }

    return(0);
}

void printDropletGlobalHeader(UserData data,FILE* output){
    fprintf(output,"%15s\t","#1");
    fprintf((output), "%15s\t%15s\t%15s\t","2","3","4");
    fprintf((output), "\n");

    fprintf((output), "%15s\t","#t");
    fprintf((output), "%15s\t%15s\t%15s\t","Mass (kg)","Radius (m)","surfaceT (K)");
    fprintf((output), "\n");
}

void printDropletGlobalOutput(UserData data,FILE* output,double t){
    double T = data->interfaceGasCellArr[1];
    fprintf(output,"%15.6e\t",t);
    fprintf((output), "%15.6e\t%15.6e\t%15.6e\t",data->dropletMass,interfaceR,T);
    fprintf((output), "\n");
}


void printSpaceTimeHeader(UserData data, FILE* output)
{
	fprintf((output), "%15s\t","#1");
	for (size_t k = 1; k <=data->nvar+2; k++) {
		fprintf((output), "%15lu\t",k+1);
	}
	fprintf((output), "\n");

	fprintf((output), "%15s\t%15s\t%15s\t","#psi","time(s)","dpsi");
	fprintf((output), "%15s\t%15s\t","radius(m)","Temp(K)");
	for (size_t k = 1; k <=data->nsp; k++) {
		fprintf((output), "%15s\t",gas->speciesName(k-1).c_str());
	}
	fprintf((output), "%15s\n","Pressure(Pa)");
}

void printSpaceTimeOutput(double t, N_Vector* y, FILE* output, UserData data)
{
	double *ydata,*psidata;
  	ydata    = N_VGetArrayPointer_OpenMP(*y);

	if(data->adaptiveGrid){
  		psidata  = data->grid->x;
	}else{
  		psidata  = data->uniformGrid;
	}

	for (size_t i = 0; i < data->npts; i++) {
		fprintf(output, "%15.6e\t%15.6e\t",psi(i+1),t);
		if(i==0){
			fprintf(output, "%15.6e\t",psi(2)-psi(1));
		}
		else{
			fprintf(output, "%15.6e\t",psi(i+1)-psi(i));
		}
		fprintf(output, "%15.6e\t",R(i+1));
		for (size_t j = 0; j < data->nvar; j++) {
			fprintf(output, "%15.6e\t",ydata[j+i*data->nvar]);
		}
		fprintf(output, "\n");
	}
	fprintf(output, "\n\n");
}

void printSpaceTimeOutput(double t, double* ydata, FILE* output, UserData data)
{
	//double *ydata,*psidata;
  	//ydata    = N_VGetArrayPointer_OpenMP(*y);
	double *psidata;

	if(data->adaptiveGrid){
  		psidata  = data->grid->x;
	}else{
  		psidata  = data->uniformGrid;
	}

	for (size_t i = 0; i < data->npts; i++) {
		fprintf(output, "%15.6e\t%15.6e\t",psi(i+1),t);
		if(i==0){
			fprintf(output, "%15.6e\t",psi(2)-psi(1));
		}
		else{
			fprintf(output, "%15.6e\t",psi(i+1)-psi(i));
		}
		fprintf(output, "%15.9e\t",R(i+1));
		for (size_t j = 0; j < data->nvar; j++) {
			fprintf(output, "%15.9e\t",ydata[j+i*data->nvar]);
		}
		fprintf(output, "\n");
	}
	fprintf(output, "\n\n");
}

void printSpaceTimeOutputNew(double t,double* ydata,FILE* output,UserData data){
    double *psidata;

    if(data->adaptiveGrid){
        psidata  = data->grid->x;
    }else{
        psidata  = data->uniformGrid;
    }

    for (size_t i=0;i<data->nlpts;i++){  // print liquid phase information
        fprintf(output,"%15.6e\t%15.6e\t", eta(i+1),t);
        if(i==0){
            fprintf(output, "%15.6e\t",std::abs(eta(2)-eta(1)));
        }
        else{
            fprintf(output, "%15.6e\t",std::abs(eta(i+1)-eta(i)));
        }
        fprintf(output,"%15.6e\t", liquidR(i+1));
        for (size_t j = 0; j < data->nvar; j++) {
            fprintf(output, "%15.6e\t",ydata[j+i*data->nvar]);
        }
        fprintf(output,"\n");
    }

    for(size_t i=data->nlpts+1; i<=data->nlpts+data->npts;i++){ // print gas phase information
        size_t gasII = i - data->nlpts+1;
        fprintf(output,"%15.6e\t%15.6e\t", psi(gasII),t);
        fprintf(output,"%15.6e\t", psi(gasII) - psi(gasII-1));
        fprintf(output,"%15.6e\t", gasR(gasII)) ;
        for (size_t j = 0; j < data->nvar; j++) {
            fprintf(output, "%15.6e\t",ydata[j+(i-1)*data->nvar]);
        }
        fprintf(output,"\n");
    }

    fprintf(output, "\n\n");
}


void writeRestart(double t, N_Vector* y, FILE* output, UserData data){
	double *ydata,*psidata;
	ydata = N_VGetArrayPointer_OpenMP(*y);
	if(data->adaptiveGrid){
  		psidata  = data->grid->x;
	}else{
  		psidata  = data->uniformGrid;
	}
 	fwrite(&t, sizeof(t), 1, output);		//write time
 	fwrite(psidata, data->npts*sizeof(psidata), 1, output);	//write grid
 	fwrite(ydata, data->neq*sizeof(ydata), 1, output);	//write solution
}

void readRestart(N_Vector* y, N_Vector* ydot, FILE* input, UserData data){
	double *ydata,*psidata, *ydotdata;
	double t;
	if(data->adaptiveGrid){
  		psidata  = data->grid->x;
	}else{
  		psidata  = data->uniformGrid;
	}
	ydata = N_VGetArrayPointer_OpenMP(*y);
	ydotdata = N_VGetArrayPointer_OpenMP(*ydot);
 	fread(&t, sizeof(t), 1, input);
	data->tNow=t;
 	fread(psidata, data->npts*sizeof(psidata), 1, input);
 	fread(ydata, data->neq*sizeof(ydata), 1, input);
 	fread(ydotdata, data->neq*sizeof(ydotdata), 1, input);
	if(data->adaptiveGrid){
		storeGrid(data->grid->x,data->grid->xOld,data->npts);
	}
}

//void printGlobalHeader(UserData data)
//{
//	fprintf((data->globalOutput), "%8s\t","#Time(s)");
//	//fprintf((data->globalOutput), "%15s","S_u(exp*)(m/s)");
//	fprintf((data->globalOutput), "%15s","bFlux(kg/m^2/s)");
//	fprintf((data->globalOutput), "%16s","  IsothermPos(m)");
//	fprintf((data->globalOutput), "%15s\t","Pressure(Pa)");
//	fprintf((data->globalOutput), "%15s\t","Pdot(Pa/s)");
//	fprintf((data->globalOutput), "%15s\t","gamma");
//	fprintf((data->globalOutput), "%15s\t","S_u(m/s)");
//	fprintf((data->globalOutput), "%15s\t","Tu(K)");
//	fprintf((data->globalOutput), "\n");
//}
//
//void printSpaceTimeRates(double t, N_Vector ydot, UserData data)
//{
//	double *ydotdata,*psidata;
//  	ydotdata    = N_VGetArrayPointer_OpenMP(ydot);
//  	psidata    = N_VGetArrayPointer_OpenMP(data->grid);
//	for (int i = 0; i < data->npts; i++) {
//		fprintf((data->ratesOutput), "%15.6e\t%15.6e\t",psi(i+1),t);
//		for (int j = 0; j < data->nvar; j++) {
//			fprintf((data->ratesOutput), "%15.6e\t",ydotdata[j+i*data->nvar]);
//		}
//		fprintf((data->ratesOutput), "\n");
//	}
//	fprintf((data->ratesOutput), "\n\n");
//}
//
//
//void printGlobalVariables(double t, N_Vector* y, N_Vector* ydot, UserData data)
//{
//	double *ydata,*ydotdata, *innerMassFractionsData, *psidata;
//	innerMassFractionsData =  data->innerMassFractions;
//
//	if(data->adaptiveGrid){
//  		psidata  = data->grid->x;
//	}else{
//  		psidata  = data->uniformGrid;
//	}
//
//	double TAvg, RAvg, YAvg, psiAvg;
//  	ydata    = N_VGetArrayPointer_OpenMP(*y);
//  	ydotdata = N_VGetArrayPointer_OpenMP(*ydot);
//	TAvg=data->isotherm;
//	double sum=ZERO;
//	double dpsim,area,aream,drdt;
//	double Cpb,Cvb,gamma,rho,flameArea,Tu;
//
//
//	/*Find the isotherm chosen by the user*/
//	size_t j=1;
//	size_t jj=1;
//	size_t jjj=1;
//	double wdot[data->nsp];
//	double wdotMax=0.0e0;
//	double advTerm=0.0e0;
//	psiAvg=0.0e0;
//	if(T(data->npts)>T(1)){
//		while (T(j)<TAvg) {
//			j=j+1;
//		}
//		YAvg=innerMassFractionsData[data->k_oxidizer-1]-Y(data->npts,data->k_oxidizer);
//		while (fabs((T(jj+1)-T(jj))/T(jj))>1e-08) {
//			jj=jj+1;
//		}
//
//		setGas(data,ydata,jj,GAS);
//		Tu=T(jj);
//		rho=data->gas->density();
//        	Cpb=data->gas->cp_mass();       	//J/kg/K
//        	Cvb=data->gas->cv_mass();       	//J/kg/K
//		gamma=Cpb/Cvb;
//
//	}
//	else{
//		while (T(j)>TAvg) {
//			j=j+1;
//		}
//		YAvg=innerMassFractionsData[data->k_oxidizer-1]-Y(1,data->k_oxidizer);
//		while (fabs((T(data->npts-jj-1)-T(data->npts-jj))/T(data->npts-jj))>1e-08) {
//			jj=jj+1;
//		}
//
//		setGas(data,ydata,data->npts-jj,GAS);
//		Tu=T(data->npts-jj);
//		rho=data->gas->density();
//        	Cpb=data->gas->cp_mass();       	//J/kg/K
//        	Cvb=data->gas->cv_mass();       	//J/kg/K
//		gamma=Cpb/Cvb;
//	}
//	if(T(j)<TAvg){
//		RAvg=((R(j+1)-R(j))/(T(j+1)-T(j)))*(TAvg-T(j))+R(j);
//	}
//	else{
//		RAvg=((R(j)-R(j-1))/(T(j)-T(j-1)))*(TAvg-T(j-1))+R(j-1);
//	}
//
//	////Experimental burning speed calculation:
//	//int nMax=0;
//	////nMax=maxCurvIndex(ydata, data->nt, data->nvar,
//	////     	data->grid->x, data->npts);
//	//nMax=maxGradIndex(ydata, data->nt, data->nvar,
//	//     	data->grid->x, data->npts);
//	//advTerm=(T(nMax)-T(nMax-1))/(data->mass*(psi(nMax)-psi(nMax-1)));
//	//aream=calc_area(R(nMax),&data->metric);
//	////setGas(data,ydata,nMax);
//	////rho=data->gas->density();
//	//psiAvg=-Tdot(nMax)/(rho*aream*advTerm);
//	////if(t>data->ignTime){
//	////	for(size_t n=2;n<data->npts;n++){
//	////		setGas(data,ydata,n);
//	////		data->gas->getNetProductionRates(wdot); //kmol/m^3
//	////		advTerm=(T(n)-T(n-1))/(data->mass*(psi(n)-psi(n-1)));
//	////		if(fabs(wdot[data->k_oxidizer-1])>=wdotMax){
//	////			aream=calc_area(R(n),&data->metric);
//	////			psiAvg=-Tdot(n)/(rho*aream*advTerm);
//	////			wdotMax=fabs(wdot[data->k_oxidizer-1]);
//	////		}
//	////	}
//	////}
//	////else{
//	////	psiAvg=0.0e0;
//	////}
//
//
//	//drdt=(RAvg-data->flamePosition[1])/(t-data->flameTime[1]);
//	//data->flamePosition[0]=data->flamePosition[1];
//	//data->flamePosition[1]=RAvg;
//	//data->flameTime[0]=data->flameTime[1];
//	//data->flameTime[1]=t;
//	//flameArea=calc_area(RAvg,&data->metric);
//
//	/*Use the Trapezoidal rule to calculate the mass burning rate based on
//	 * the consumption of O2*/
//	aream= calc_area(R(1)+1e-03*data->domainLength,&data->metric);
//	for (j = 2; j <data->npts; j++) {
//        	dpsim=(psi(j)-psi(j-1))*data->mass;
//		area= calc_area(R(j),&data->metric);
//		sum=sum+HALF*dpsim*((Ydot(j-1,data->k_oxidizer)/aream)
//				   +(Ydot(j,data->k_oxidizer)/area));
//		aream=area;
//	}
//
//	//double maxOH,maxHO2;
//	//maxOH=0.0e0;
//	//maxHO2=0.0e0;
//	//for(j=1;j<data->npts;j++){
//	//	if(Y(j,data->k_OH)>maxOH){
//	//		maxOH=Y(j,data->k_OH);
//	//	}
//	//}
//	//for(j=1;j<data->npts;j++){
//	//	if(Y(j,data->k_HO2)>maxHO2){
//	//		maxHO2=Y(j,data->k_HO2);
//	//	}
//	//}
//
//	fprintf((data->globalOutput), "%15.6e\t",t);
//	//fprintf((data->globalOutput), "%15.6e\t",psiAvg);
//	fprintf((data->globalOutput), "%15.6e\t",fabs(sum)/YAvg);
//	fprintf((data->globalOutput), "%15.6e\t",RAvg);
//	fprintf((data->globalOutput), "%15.6e\t",P(data->npts));
//	fprintf((data->globalOutput), "%15.6e\t",Pdot(data->npts));
//	fprintf((data->globalOutput), "%15.6e\t",gamma);
//	fprintf((data->globalOutput), "%15.6e\t",fabs(sum)/(YAvg*rho));
//	fprintf((data->globalOutput), "%15.6e\t",Tu);
//	fprintf((data->globalOutput), "\n");
//}

void getR(N_Vector y, UserData data){

   //TODO: Needs testing
   double *ydata, *psidata;
   
	ydata 	= N_VGetArrayPointer_OpenMP(y); 
	psidata = data->uniformGrid;

   double dpsim, mass, rhom, rho, rhomh;

   mass = data->mass;

   R(1) = data->Rd;
//    setGas(data,ydata,1,GAS);
	setGas(data,ydata,1);
   rhom = gas->density();

   for(size_t gridII = 2; gridII <= data->npts; gridII++){
       setGas(data,ydata,gridII);
      rho = gas->density();
      rhomh = (rho + rhom)/TWO;

      dpsim = (psi(gridII)   - psi(gridII-1));                            

      if(data->metric == 0){
         R(gridII) = dpsim * mass / rhomh + R(gridII-1);
      }else if(data->metric == 1){
         R(gridII) = std::sqrt(TWO * dpsim * mass / rhomh + R(gridII-1)*R(gridII-1));
      }else if(data->metric == 2){
         R(gridII) = std::cbrt(THREE * dpsim * mass / rhomh + R(gridII-1)*R(gridII-1)*R(gridII-1));
      }

      rhom = rho;
   }
}


void getR(double* ydata, UserData data){

    //TODO: Needs testing
    double  *psidata;

//    ydata 	= N_VGetArrayPointer_OpenMP(y);
    psidata = data->uniformGrid;

    //DEBUG
    FILE* outYdataFile = fopen("checkYdataAgain.txt", "w");
    if (outYdataFile == nullptr) {
        printf(("allocate outYdataFile failed ! \n"));
    }
    // Write array data to file using fprintf
    for (int i = 0; i < data->npts; ++i) {
        for (int j = 0; j < data->nvar; ++j) {
            fprintf(outYdataFile, "%15.6e\t", ydata[i*data->nvar+j]);  // Print each element with 2 decimal places
        }
        fprintf(outYdataFile, "\n");                     // Newline after each row
    }
    // Close the file
    fclose(outYdataFile);

    double dpsim, mass, rhom, rho, rhomh;

    mass = data->mass;

    R(1) = data->Rd;
//    setGas(data,ydata,1,GAS);
    setGas(data,ydata,1);
    rhom = gas->density();

    for(size_t gridII = 2; gridII <= data->npts; gridII++){
        setGas(data,ydata,gridII);
        rho = gas->density();
        rhomh = (rho + rhom)/TWO;

        dpsim = (psi(gridII)   - psi(gridII-1));

        if(data->metric == 0){
            R(gridII) = dpsim * mass / rhomh + R(gridII-1);
        }else if(data->metric == 1){
            R(gridII) = std::sqrt(TWO * dpsim * mass / rhomh + R(gridII-1)*R(gridII-1));
        }else if(data->metric == 2){
            R(gridII) = std::cbrt(THREE * dpsim * mass / rhomh + R(gridII-1)*R(gridII-1)*R(gridII-1));
        }

        rhom = rho;
    }
}


void getRNew(double* ydata, UserData data){
    //DEBUG
    debugYdata(ydata,data) ;

    //TODO: Needs testing
    double  *psidata;

//    ydata 	= N_VGetArrayPointer_OpenMP(y);
    psidata = data->uniformGrid;
    size_t npts = data->npts;
    size_t nlpts = data->nlpts;

//    //DEBUG
//    FILE* outYdataFile = fopen("checkYdataAgain.txt", "w");
//    if (outYdataFile == nullptr) {
//        printf(("allocate outYdataFile failed ! \n"));
//    }
//    // Write array data to file using fprintf
//    for (int i = 0; i < data->npts; ++i) {
//        for (int j = 0; j < data->nvar; ++j) {
//            fprintf(outYdataFile, "%15.6e\t", ydata[i*data->nvar+j]);  // Print each element with 2 decimal places
//        }
//        fprintf(outYdataFile, "\n");                     // Newline after each row
//    }
//    // Close the file
//    fclose(outYdataFile);

    double dpsim,detam, mass, dropletMass , rhom, rho, rhomh;

    mass = data->mass;
    dropletMass = data->dropletMass ;

    // Update the liquid phase R
    liquidR(1) = ZERO;
    singleLiquidFuel->setLiquidState(T(1), P(1));
    rhom = singleLiquidFuel->calculateDensity();

    for (size_t gridII = 2; gridII <= nlpts + 1; ++gridII) {
        double T_current = (gridII <= nlpts) ? T(gridII) : data->interfaceLiquidCellArr[1];
        double P_current = (gridII <= nlpts) ? P(gridII) : data->interfaceLiquidCellArr[data->nvar];

        singleLiquidFuel->setLiquidState(T_current, P_current);
        rho = singleLiquidFuel->calculateDensity();
        rhomh = (rho + rhom) / TWO;

        detam = std::abs(eta(gridII) - eta(gridII - 1));
        if (data->metric == 2) {
            liquidR(gridII) = std::cbrt(THREE * detam * dropletMass / rhomh +
                                        liquidR(gridII - 1) * liquidR(gridII - 1) * liquidR(gridII - 1));
        } else {
            printf("Make sure metric equals to 2.\n");
            break; // Exit the loop if metric is not 2
        }
        rhom = rho;
    }

    //update the droplet radius
    data->Rd = liquidR(nlpts+1);

    // update both interface cell array's R
    data->interfaceGasCellArr[0] = data->Rd;
    data->interfaceLiquidCellArr[0] = data->Rd;

    // second step:we update the gas phase R array
//    R(1) = data->Rd;
    gasR(1) = data->Rd;
//    setGas(data,ydata,1,GAS);
//    setGas(data,ydata,1);
    gas->setState_TPY(data->interfaceGasCellArr[1],data->interfaceGasCellArr[data->nvar],&data->interfaceGasCellArr[2]);
    rhom = gas->density();

    for(size_t gridII = 2; gridII <= (npts+1); gridII++){
        size_t gasII = gridII + nlpts-1;
        // DEBUG
//        double massFracSum = 0.0;
//        for (size_t k = 0;k<data->nsp;k++){
//            massFracSum += Y(gasII,k+1);
//        }
//        printf("gasII = %zu , T = %.2f [K] , P = %.2f [Pa] , Sum of mass fraction : %.2f .\n", gasII, T(gasII), P(gasII),massFracSum);
        gas->setState_TPY(T(gasII), P(gasII), &Y(gasII,1));
//        setGas(data,ydata,gridII);
        rho = gas->density();
        rhomh = (rho + rhom)/TWO;

        dpsim = std::abs((psi(gridII) - psi(gridII-1)));

        if(data->metric == 0){
            gasR(gridII) = dpsim * mass / rhomh + gasR(gridII-1);
        }else if(data->metric == 1){
            gasR(gridII) = std::sqrt(TWO * dpsim * mass / rhomh + gasR(gridII-1)*gasR(gridII-1));
        }else if(data->metric == 2){
            gasR(gridII) = std::cbrt(THREE * dpsim * mass / rhomh + gasR(gridII-1)*gasR(gridII-1)*gasR(gridII-1));
//            R(gridII) = std::cbrt(THREE * dpsim * mass / rhomh + R(gridII-1)*R(gridII-1)*R(gridII-1));
        }
        rhom = rho;
    }

}
//
//void printSpaceTimeOutputInterpolated(double t, N_Vector y, UserData data)
//{
//	double *ydata,*psidata;
//  	ydata    = N_VGetArrayPointer_OpenMP(y);
//  	psidata    = N_VGetArrayPointer_OpenMP(data->grid);
//	for (int i = 0; i < data->npts; i++) {
//		fprintf((data->gridOutput), "%15.6e\t%15.6e\t",psi(i+1),t);
//		for (int j = 0; j < data->nvar; j++) {
//			fprintf((data->gridOutput), "%15.6e\t",ydata[j+i*data->nvar]);
//		}
//		fprintf((data->gridOutput), "\n");
//	}
//	fprintf((data->gridOutput), "\n\n");
//}
//
//
////void repairSolution(N_Vector y, N_Vector ydot, UserData data){
////	int npts=data->npts;
////  	double *ydata;
////  	double *ydotdata;
////	ydata    = N_VGetArrayPointer_OpenMP(y); 
////	ydotdata    = N_VGetArrayPointer_OpenMP(ydot); 
////
////	T(2)=T(1);
////	T(npts-1)=T(npts);
////	Tdot(2)=Tdot(1);
////	Tdot(npts-1)=Tdot(npts);
////	for (int k = 1; k <=data->nsp; k++) {
////		    Y(2,k)=Y(1,k);
////		    Y(npts-1,k)=Y(npts,k);
////
////		    Ydot(2,k)=Ydot(1,k);
////		    Ydot(npts-1,k)=Ydot(npts,k);
////	}
////}

//void isobaricAdvance(double* ydata, UserData data){
//
//    size_t nsp = data->nsp ;
//    size_t npts = data->npts;
//    size_t nvar = data->nvar;
//    double dt = data->deltaT/TWO;
//    clock_t start,end;
//    start=clock();
//
//    /*Allocate reactors:*/
////    size_t nReactors=data->nThreads;
////    isobaricReactor_t *reactor[nReactors];
////    reactorNet_t *sim[nReactors];
////    for (size_t i=0;i<nReactors;i++) {
////        reactor[i]=isobaricReactorCreate();
//////        addGasToIsobaricReactor(params->gas[i],reactor[i]);
////        addGasToIsobaricReactor(data->gas,reactor[i]); //pressure is inherited implicitly from gas
////        sim[i] = reactorNetCreate();
////        addIsobaricReactorToNet(reactor[i],sim[i]);
////        reactorNetSetTolerances(sim[i],data->relativeTolerance,data->massFractionTolerance);
////    }
//    data->constPressureReactor = new Cantera::IdealGasConstPressureReactor() ;
//    data->reactorNet = new Cantera::ReactorNet() ;
//
//    addGasToIsobaricReactor(data);
//    addIsobaricReactorToNet(data);
//    reactorNetSetTolerances(data);
//
//    double w[nsp+2];
//
//    /*loop over all grid points*/
//    for (size_t gridII = 1;gridII <= npts; gridII++){
//        w[0] = 1.0e0;
//        w[1] = T(gridII);
//        for (size_t speciesII = 1;speciesII <= nsp; speciesII++){
//            w[1+speciesII] = Y(gridII,speciesII);
//        }
//        setGas(data,ydata,gridII);
//        data->constPressureReactor->syncState();
//        /*since nTreads always equals to 1,
//         * so only 1 isobaric reactor will be used
//         * throughout the loop*/
////        setIsobaricReactorState(data,w);
//        setNetTime(data,0.0);
//        reactorNetAdvance(data,dt) ;
//        getIsobaricReactorState(data,w);
////        setIsobaricReactorState(reactor[0], w);
////        setNetTime(sim[0],0.0);
////        reactorNetAdvance(sim[0],dt);
////        getIsobaricReactorState(reactor[0], w);
//
//        T(gridII) = w[1];
//        for (size_t speciesII = 1;speciesII <= nsp; speciesII++){
//            Y(gridII,speciesII) = w[1+speciesII];
//        }
//    }
//
////    if(T[]>=params->TCutOff){
////        w[0]=1.0e0;
////        w[1]=T[];
////        for (size_t i = 0; i < params->nSpecies; i++)
////            w[i+2]=Yi[i];
////        setIsobaricReactorState(reactor[j], w);
////        setNetTime(sim[j],0.0);
////        reactorNetAdvance(sim[j],dt);
////        getIsobaricReactorState(reactor[j], w);
////        T[]=w[1];
////        k=0;
////        for (scalar s in Y){
////            s[]=w[2+k];
////            Yi[k]=w[2+k];
////            k++;
////        }
////    }
//    delete data->constPressureReactor;
//    delete data->reactorNet;
//    end = clock();
//    printf("isobaricAdvance time: %f seconds\n",(double)(end-start) / CLOCKS_PER_SEC);
//}


// ... [Other necessary includes and definitions]

//void isobaricAdvance2(double* ydata, UserData data) {
//    size_t nsp = data->nsp;
//    size_t npts = data->npts;
//    double dt = data->deltaT / 2.0;
//    clock_t start, end;
//    start = clock();
//
//    // Number of threads
//    int nThreads = omp_get_max_threads();
//
//    // Allocate reactors, reactor networks, and gases per thread
//    std::vector<Cantera::IdealGasConstPressureReactor*> reactors(nThreads);
//    std::vector<Cantera::ReactorNet*> reactorNets(nThreads);
//    std::vector<Cantera::IdealGasMix*> gases(nThreads);
//
//    // Initialize reactors and reactor networks
//#pragma omp parallel for
//    for (int i = 0; i < nThreads; ++i) {
//        gases[i] = new Cantera::IdealGasMix(data->model);  // Clone gas for thread safety
//        reactors[i] = new Cantera::IdealGasConstPressureReactor();
//        reactors[i]->insert(* gases[i]);
//        reactorNets[i] = new Cantera::ReactorNet();
//        reactorNets[i]->addReactor(*reactors[i]);
//        reactorNets[i]->setTolerances(data->relativeTolerance, data->massFractionTolerance);
//    }
//
//    // Prepare storage for state vectors
//    std::vector<std::vector<double>> wThread(nThreads, std::vector<double>(nsp + 2));
//
//    // Parallel loop over grid points
//#pragma omp parallel for
//    for (size_t gridII = 1; gridII <= npts; ++gridII) {
//        int threadID = omp_get_thread_num();
//        auto& reactor = reactors[threadID];
//        auto& reactorNet = reactorNets[threadID];
//        auto& gas = gases[threadID];
//        auto& w = wThread[threadID];
//
//        // Prepare the state vector
//        w[0] = 1.0;
//        w[1] = T(gridII);
//        for (size_t speciesII = 1; speciesII <= nsp; ++speciesII) {
//            w[1 + speciesII] = Y(gridII, speciesII);
//        }
//
//        // Set reactor state efficiently
//        gas->setState_TPY(w[1], data->initialPressure * Cantera::OneAtm, &w[2]);
//        reactor->syncState(); //sync reactor with gas state
//
//        // Advance the reactor
//        reactorNet->setInitialTime(0.0);
//        reactorNet->advance(dt);
//
//        // Get updated state
//        reactor->getState(&w[0]);
////        w[1] = gas->temperature();
////        gas->getMassFractions(&w[2]);
//
//        // Update temperature and species mass fractions
//        T(gridII) = w[1];
//        for (size_t speciesII = 1; speciesII <= nsp; ++speciesII) {
//            Y(gridII, speciesII) = w[1 + speciesII];
//        }
//    }
//
//    // Clean up
//    for (int i = 0; i < nThreads; ++i) {
//        delete reactors[i];
//        delete reactorNets[i];
//        delete gases[i];
//    }
//    end = clock();
//    printf("isobaricAdvance time: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
//}

//
//
//isobaricReactor_t *isobaricReactorCreate(){
//    isobaricReactor_t *r;
//    Cantera::IdealGasConstPressureReactor *obj;
//    r      = (decltype(r))malloc(sizeof(*r));
//    obj = new Cantera::IdealGasConstPressureReactor();
//    r->obj = obj;
//    return r;
//}
//
//reactorNet_t *reactorNetCreate(){
//    reactorNet_t *sim;
//    Cantera::ReactorNet *obj;
//    sim      = (decltype(sim))malloc(sizeof(*sim));
//    obj = new Cantera::ReactorNet();
//    sim->obj = obj;
//    return sim;
//}
//
//
////void idealGasDestroy(idealGas_t *gas){
////    if(gas==NULL) return;
////    delete static_cast<IdealGasMix *>(gas->obj);
////    free(gas);
////    printf("gas object freed!\n");
////}
//
////
////void isochoricReactorDestroy(isochoricReactor_t *r){
////    if(r==NULL) return;
////    delete static_cast<IdealGasReactor *>(r->obj);
////    free(r);
////}
//
//void isobaricReactorDestroy(isobaricReactor_t *r){
//    if(r==NULL) return;
//    delete static_cast<Cantera::IdealGasConstPressureReactor *>(r->obj);
//    free(r);
//}
//
//void reactorNetDestroy(reactorNet_t *sim){
//    if(sim==NULL) return;
//    delete static_cast<Cantera::ReactorNet *>(sim->obj);
//    free(sim);
//}
//
//
////void addGasToIsochoricReactor(idealGas_t *gas, isochoricReactor_t *r)
////{
////    IdealGasMix *gasobj;
////    IdealGasReactor *reactorobj;
////    gasobj = static_cast<IdealGasMix *>(gas->obj);
////    reactorobj = static_cast<IdealGasReactor *>(r->obj);
////    reactorobj->insert(*gasobj);
////}
//
////void addGasToIsobaricReactor(idealGas_t *gas, isobaricReactor_t *r)
////{
////    IdealGasMix *gasobj;
////    IdealGasConstPressureReactor *reactorobj;
////    gasobj = static_cast<IdealGasMix *>(gas->obj);
////    reactorobj = static_cast<IdealGasConstPressureReactor *>(r->obj);
////    reactorobj->insert(*gasobj);
////}
//
///* modified version for LTORC */
//void addGasToIsobaricReactor(Cantera::IdealGasMix *gas, isobaricReactor_t *r)
//{
//    Cantera::IdealGasMix *gasobj;
//    Cantera::IdealGasConstPressureReactor *reactorobj;
////    gasobj = static_cast<IdealGasMix *>(gas->obj);
//    gasobj = gas;
//    reactorobj = static_cast<Cantera::IdealGasConstPressureReactor *>(r->obj);
//    reactorobj->insert(*gasobj);
//}
//
////void setIsochoricReactorState(isochoricReactor_t *r, double* y)
////{
////    IdealGasReactor *reactorobj;
////    reactorobj = static_cast<IdealGasReactor *>(r->obj);
////    reactorobj->updateState(y);
////}
//
//void setIsobaricReactorState(isobaricReactor_t *r, double* y)
//{
//    Cantera::IdealGasConstPressureReactor *reactorobj;
//    reactorobj = static_cast<Cantera::IdealGasConstPressureReactor *>(r->obj);
//    reactorobj->updateState(y);
//}
//
////void getIsochoricReactorState(isochoricReactor_t *r, double* y)
////{
////    IdealGasReactor *reactorobj;
////    reactorobj = static_cast<IdealGasReactor *>(r->obj);
////    reactorobj->getState(y);
////}
//
//void getIsobaricReactorState(isobaricReactor_t *r, double* y)
//{
//    Cantera::IdealGasReactor *reactorobj;
//    reactorobj = static_cast<Cantera::IdealGasReactor *>(r->obj);
//    reactorobj->getState(y);
//}
//
//void setNetTime(reactorNet_t *sim, double t)
//{
//    Cantera::ReactorNet* reactornetobj;
//    reactornetobj = static_cast<Cantera::ReactorNet *>(sim->obj);
//    reactornetobj->setInitialTime(t);
//}
//
////void addIsochoricReactorToNet(isochoricReactor_t *r, reactorNet_t *sim)
////{
////    IdealGasReactor *reactorobj;
////    ReactorNet* reactornetobj;
////    reactorobj = static_cast<IdealGasReactor *>(r->obj);
////    reactornetobj = static_cast<ReactorNet *>(sim->obj);
////    reactornetobj->addReactor(*reactorobj);
////}
//
//void addIsobaricReactorToNet(isobaricReactor_t *r, reactorNet_t *sim)
//{
//    Cantera::IdealGasConstPressureReactor *reactorobj;
//    Cantera::ReactorNet* reactornetobj;
//    reactorobj = static_cast<Cantera::IdealGasConstPressureReactor *>(r->obj);
//    reactornetobj = static_cast<Cantera::ReactorNet *>(sim->obj);
//    reactornetobj->addReactor(*reactorobj);
//}
//
//void reactorNetAdvance(reactorNet_t *sim, double t)
//{
//    Cantera::ReactorNet *obj;
//    obj = static_cast<Cantera::ReactorNet *>(sim->obj);
//    obj->advance(t);
//}
//
//
////void reactorNetReinitialize(reactorNet_t *sim)
////{
////    ReactorNet *obj;
////    obj = static_cast<ReactorNet *>(sim->obj);
////    obj->reinitialize();
////}
//
//
//void reactorNetSetTolerances(reactorNet_t *sim, double rtol, double atol)
//{
//    Cantera::ReactorNet *obj;
//    obj = static_cast<Cantera::ReactorNet *>(sim->obj);
//    obj->setTolerances(rtol,atol);
//}
//


//void addGasToIsobaricReactor(UserData  data){
//    data->constPressureReactor->insert(*data->gas);
//}

//void addIsobaricReactorToNet(UserData data){
//    data->reactorNet->addReactor(*data->constPressureReactor);
//}

//void reactorNetSetTolerances(UserData data){
//    double rtol = data->relativeTolerance;
//    double atol = data->massFractionTolerance;
//    data->reactorNet->setTolerances(rtol,atol);
//}

//void setIsobaricReactorState(UserData data, double* y){
//    data->constPressureReactor->updateState(y);
//}
//
//void setNetTime(UserData data, double t){
//    data->reactorNet->setInitialTime(t);
//}
//
//void reactorNetAdvance(UserData data,double t){
//    data->reactorNet->advance(t);
//}
//
//void getIsobaricReactorState(UserData data,double *y){
//    data->constPressureReactor->getState(y);
//}

void debugFuncReadInitialFile(double* ydata, UserData data){
    FILE* file1 = fopen("initYdata.dat","w"); // print ydata to debug
    for(size_t i = 0; i< (data->npts+data->nlpts);i++){
        for(size_t j =0;j<data->nvar;j++){
            fprintf(file1,"%15.6e\t",ydata[j+i*data->nvar]);
        }
        fprintf(file1,"\n");
    }
    fclose(file1) ;

    FILE * file2= fopen("initRdata.dat","w"); // print Rdata to debug
    for (size_t i=0;i<(data->npts + data->nlpts+1) ; i++){
        fprintf(file2,"%15.6e\n",data->Rdata[i]);
    }
    fprintf(file2,"\n");
    fclose(file2);

    FILE* file3 = fopen("interfaceGasCellArraydata.dat","w");
    for(size_t j=0;j< (data->nvar+1);j++){
        fprintf(file3,"%15.6e\n",data->interfaceGasCellArr[j]);
    }
    fprintf(file3,"\n");
    fclose(file3);

    FILE* file4 = fopen("interfaceLiquidCellArraydata.dat","w");
    for(size_t j=0;j< (data->nvar+1);j++){
        fprintf(file4,"%15.6e\n",data->interfaceLiquidCellArr[j]);
    }
    fprintf(file4,"\n");
    fclose(file4);
}

void debugFuncInitialPsiEtaGrid(double* psidata, UserData data){
    double nlpts = data->nlpts;
    double npts = data->npts;
    FILE* file = fopen("initialPsiEta.dat","w");
    for(size_t i=0;i< nlpts+npts+1; i++){
        fprintf(file,"%15.6e\n",data->uniformGrid[i]);
    }
    fprintf(file,"\n");
    fclose(file);
}

void debugYdata(double* ydata, UserData data){
    FILE * file = fopen("debugYdata.dat","w");
    for(size_t i = 0; i< (data->npts+data->nlpts);i++){
        for(size_t j =0;j<data->nvar;j++){
            fprintf(file,"%15.6e\t",ydata[j+i*data->nvar]);
        }
        fprintf(file,"\n");
    }
    fclose(file) ;
}

void debugRdata(UserData data){
    FILE* file = fopen("debugRdata.dat","w");
    for (size_t i=0;i<(data->npts + data->nlpts+1) ; i++){
        fprintf(file,"%15.6e\n",data->Rdata[i]);
    }
    fprintf(file,"\n");
    fclose(file);
}

double heptaneVaporPressure(double temperature){
    double A = 4.02832;
    double B = 1268.636;
    double C = -56.199;
    double P= ZERO;
    P = 1.0e5 * std::pow(10.0,A-(B/(temperature+C))) ;
    return P;
}

double octaneVaporPressure(double temperature){
    double A = 3.93679;
    double B = 1257.84;
    double C = -52.415;
    double P= ZERO;
    P = 1.0e5 * std::pow(10.0,A-(B/(temperature+C))) ;
    return P;
}

//double heptaneLatentHeat(double temperature){
//    double l_ref = 3.17e5 ;
//    double T_c = 540.2 ;
//    double l = l_ref * std::pow((1-temperature/T_c),0.38) ;
//    return l;
//}

double heptaneLatentHeat(double T) {
    // Get enthalpy of saturated vapor (J/kg)
    // double Hvap = CoolProp::PropsSI("H", "T", T, "Q", 1, "n-Heptane");

    // // Get enthalpy of saturated liquid (J/kg)
    // double Hliq = CoolProp::PropsSI("H", "T", T, "Q", 0, "n-Heptane");

    // Compute latent heat of vaporization
    return 312400.00; // J/kg

    // return (Hvap - Hliq);  // J/kg
}

double octaneLatentHeat(double T) {
    // Get enthalpy of saturated vapor (J/kg)
    double Hvap = CoolProp::PropsSI("H", "T", T, "Q", 1, "nOctane");

    // Get enthalpy of saturated liquid (J/kg)
    double Hliq = CoolProp::PropsSI("H", "T", T, "Q", 0, "nOctane");

    // Compute latent heat of vaporization
    return (Hvap - Hliq);  // J/kg
}

void updateInterfaceMassFracArray(double temperature, double PAmbience, double* massFracArrayOld,
                                  double* massFracArrayNew, double* MWArray, const size_t nsp, size_t fuelIndex) {
    size_t fuel_II = fuelIndex - 1;
    double fuelMoleFrac = heptaneVaporPressure(temperature) / PAmbience;
    double moleFracNoNorm[nsp], moleFracNorm[nsp];
    double sumMoleFrac = 0.0, totalMass = 0.0;

    // Calculate non-normalized mole fractions
    for (size_t i = 0; i < nsp; ++i) {
        if (i != fuel_II) {
            moleFracNoNorm[i] = massFracArrayOld[i] / MWArray[i];
            sumMoleFrac += moleFracNoNorm[i];
        }else{ // i = fuel_II
            moleFracNoNorm[i] = 0.0 ;
        }
    }

    // Normalize mole fractions and include fuel mole fraction
    for (size_t i = 0; i < nsp; ++i) {
        moleFracNorm[i] = (i == fuel_II) ? fuelMoleFrac : (moleFracNoNorm[i] / sumMoleFrac) *(1-fuelMoleFrac);
        totalMass += moleFracNorm[i] * MWArray[i];
    }

    // Compute normalized mass fractions
    for (size_t i = 0; i < nsp; ++i) {
        massFracArrayNew[i] = moleFracNorm[i] * MWArray[i] / totalMass;
    }
}

double interfaceIterationResidue(double x, void* params) {
    auto* parameters = static_cast<interfaceProblemPara>(params);
    double factor = 1.0 ;

    // Extract parameters for clarity
//    const double latentHeat = 316887; // Unit: J/kg
//    const double latentHeat = octaneLatentHeat(x) *0.9 ; // units: J/kg
    const double latentHeat = heptaneLatentHeat(x) ; // units: J/kg
    const double lambdaGas = parameters->lambda_gas; // Unit: kg·m/s³·K (W/m/K)
    const double lambdaLiquid = parameters->lambda_liquid;
    const double leftR = parameters->left_R;
    const double rightR = parameters->right_R;
    const double dropletRadius = parameters->interface_R;
    const double leftT = parameters->left_T;
    const double rightT = parameters->right_T;
    const double D = parameters->gasDiffCoeff;
    const double rhoGas = parameters->gasDensity;
    double YV[parameters->nsp] ;
    double deltaR_R = rightR -  dropletRadius;

    std::vector<double> x1 = {parameters->left_R1,parameters->left_R,parameters->interface_R};
    std::vector<double> x2 = {parameters->interface_R,parameters->right_R,parameters->right_R1};
    std::vector<double> y1 = {parameters->left_T1,parameters->left_T,x};
    std::vector<double> y2 = {x,parameters->right_T,parameters->right_T1};

    double gradT_l = computeDerivative(x1,y1,dropletRadius);
    double gradT_g = computeDerivative(x2,y2,dropletRadius);

    double YArrayInterfaceNew[parameters->nsp];
    updateInterfaceMassFracArray(
            x, parameters->P, parameters->YArray_Right, YArrayInterfaceNew,
            parameters->MWArray, parameters->nsp, parameters->fuelIndex
    );

    getInterfaceMassFlux(leftT,rightT,parameters->P,YArrayInterfaceNew,parameters->YArray_Right,deltaR_R,YV) ;

    // Calculate mdot and associated terms
    const double YInterfaceFuel = YArrayInterfaceNew[parameters->fuelIndex - 1];
    const double YRightFuel = parameters->YArray_Right[parameters->fuelIndex - 1];
    const double mdot = rhoGas * (-D) * (YRightFuel - YInterfaceFuel) / (rightR - dropletRadius) /
            (1 - HALF*(YInterfaceFuel+YRightFuel)); // mdot here is the flux
//    const double mdot = YV[parameters->fuelIndex - 1]/ (rightR - dropletRadius) /
//            (1 - HALF*(YInterfaceFuel+YRightFuel )); // mdot here is the flux
//////            (1 - YInterfaceFuel); // mdot here is the flux

    const double vaporHeatTerm = factor* mdot * latentHeat ; //units: Watt

    // Calculate heat fluxes
//    const double gasHeatFlux = lambdaGas * (rightT - x) / (rightR - dropletRadius);
//    const double liquidHeatFlux = lambdaLiquid * (x - leftT) / (dropletRadius - leftR);
    const double gasHeatFlux = lambdaGas* gradT_g ;
    const double liquidHeatFlux = lambdaLiquid * gradT_l;
    const double res = (vaporHeatTerm+liquidHeatFlux-gasHeatFlux) ;

    // Calculate and return the residue
//    return std::abs(vaporHeatTerm) + std::abs(liquidHeatFlux) - std::abs(gasHeatFlux);
    return res;
}


//void updateInterfaceState(double* ydata, UserData data, double delta_t){
//    // allocate the required global parameters first
//    size_t nsp = data->nsp ;
//    size_t fuelIndex = data->dropII ; // Note: "1" based index
//    size_t P = data->initialPressure* Cantera::OneAtm ;
//
//    // allocate interface cell related parameters
//    double interfaceT = data->interfaceGasCellArr[1];
//    double interfaceYArray_Old[nsp], MWArray[nsp],gasDiffCoeffs[nsp],interfaceYArray_New[nsp];
//    double fuelDiffCoeff, lambdaGas,rhoGas;
//    for (size_t i = 0; i<nsp;i++) {
//        MWArray[i] = data->gas->molecularWeight(i) ;
//        interfaceYArray_Old[i] = data->interfaceGasCellArr[2+i] ;
//    }
//    data->gas->setState_TPY(interfaceT,P,&interfaceYArray_Old[0]) ;
//    lambdaGas = data->trmix->thermalConductivity() ; //units: W/(m*K)
//    data->trmix->getMixDiffCoeffs(&gasDiffCoeffs[0]) ; //units: m^2/s
//    fuelDiffCoeff = gasDiffCoeffs[fuelIndex -1];
//    rhoGas = data->gas->density() ; //units: kg/m^3
//    double dropletRadius = data->interfaceGasCellArr[0] ;
//
//    // allocate liquid phase related parameters
//    double leftT = T(data->nlpts) ;
//    data->singleLiquidFuel->setLiquidState(leftT,P);
//    double lambdaLiquid = data->singleLiquidFuel->calculateThermalConductivity(); //units: W/(m*K)
//    double leftR = liquidR(data->nlpts);
//
//    // allocate gas phase related parameters
//    double rightT = T(data->nlpts+1) ;
//    double rightR = gasR(2);
//    double YArrayRight[nsp] ;
//    for (size_t k =0;k<nsp;k++){
//        YArrayRight[k] =  Y(data->nlpts+1,k+1) ;
//    }
//
//    if (data->flagSolveInterfaceProblem){
//        // allocate data structure to solve interface problem
//        interfaceProblemPara nP = new interfaceProblemParameters;
//        nP->nsp = nsp ;
//        nP->fuelIndex = fuelIndex ;
//        nP->P = P;
//        nP->gasDensity = rhoGas ;
//        nP->gasDiffCoeff = fuelDiffCoeff;
//        nP->lambda_gas = lambdaGas;
//        nP->lambda_liquid = lambdaLiquid;
//        nP->left_R = leftR;
//        nP->interface_R = dropletRadius;
//        nP->right_R = rightR;
//        nP->left_T = leftT ;
//        nP->right_T = rightT ;
//        nP->YArrayMiddle_Old = new double[nsp];
//        nP->YArray_Right = new double[nsp];
//        nP->MWArray = new double[nsp] ;
//        for(size_t k=0;k<nsp;k++){
//            nP->YArrayMiddle_Old[k] = interfaceYArray_Old[k] ;
//            nP->MWArray[k] = MWArray[k];
//            nP->YArray_Right[k] = YArrayRight[k] ;
//        }
//
//        // gsl brent solver
//        double x_lo = interfaceT;
//        double x_hi = rightT*0.99;
//        const gsl_root_fsolver_type *fT;
//        gsl_root_fsolver *fS;
//        gsl_function F;
//        F.function = &interfaceIterationResidue;
//        F.params = nP;
//        fT = gsl_root_fsolver_brent ;
//        fS = gsl_root_fsolver_alloc(fT);
//        gsl_root_fsolver_set(fS,&F,x_lo,x_hi);
//
//        int status;
//        int iter=0,max_iter = 20;
//        double r = 0.00;
//        do{
//            iter++;
//            status = gsl_root_fsolver_iterate(fS);
//            r = gsl_root_fsolver_root (fS);
//            x_lo = gsl_root_fsolver_x_lower (fS);
//            x_hi = gsl_root_fsolver_x_upper (fS);
//            status = gsl_root_test_interval (x_lo, x_hi,
//                                             0, 1e-09);
//        }  while (status == GSL_CONTINUE && iter < max_iter);
//
//        // vapor pressore should not exceed ambient pressure
//        double temp_interfaceT = r;
//        double vaporPressure = heptaneVaporPressure(temp_interfaceT);
//        double area = calc_area(dropletRadius,&data->metric);
//        if (vaporPressure <= P){ // temperature lower than saturated temperature
//            updateInterfaceMassFracArray(temp_interfaceT,P,&interfaceYArray_Old[0],&interfaceYArray_New[0],&MWArray[0],nsp,fuelIndex);
//            double mdot = rhoGas*(-fuelDiffCoeff)*(interfaceYArray_New[fuelIndex-1] - YArrayRight[fuelIndex-1])/(rightR - dropletRadius)/(1-interfaceYArray_New[fuelIndex-1]) * area;
//
//            // update relevant parameters in data structure
//            data->dropletMass = data->dropletMass - std::abs(mdot) * delta_t ;
//            data->interfaceGasCellArr[1] = temp_interfaceT ; // update cell temperature
//            data->interfaceLiquidCellArr[1] = temp_interfaceT ;
//            for (size_t k = 0;k<nsp;k++){
//                data->interfaceGasCellArr[2+k] = interfaceYArray_New[k];
//            }
//        }else{
//            data->flagSolveInterfaceProblem = false;
//            double mdot = rhoGas* (-fuelDiffCoeff)*(interfaceYArray_Old[fuelIndex-1] - YArrayRight[fuelIndex-1])/(rightR - dropletRadius)/(1-interfaceYArray_Old[fuelIndex-1]) *area ;
//            data->dropletMass = data->dropletMass - std::abs(mdot) * delta_t ;
//        }
//        delete nP->YArrayMiddle_Old;
//        delete nP->YArray_Right;
//        delete nP->MWArray ;
//        delete nP;
//    }else{ // flag = false, no need to solve iteration problem
//        double area = calc_area(dropletRadius,&data->metric);
//        double mdot = rhoGas*(-fuelDiffCoeff) *(interfaceYArray_Old[fuelIndex-1] -YArrayRight[fuelIndex-1])/(rightR - dropletRadius)/(1-interfaceYArray_Old[fuelIndex-1]) *area;
//        data->dropletMass = data->dropletMass - std::abs(mdot) * delta_t ;
//    }
//}


void updateInterfaceState(double* ydata, UserData data, double delta_t) {
    // Extract reusable parameters
    size_t nsp = data->nsp;
    size_t fuelIndex = data->dropII;
    double P = data->initialPressure * Cantera::OneAtm;
    double interfaceT = data->interfaceGasCellArr[1];
    double dropletRadius = interfaceR;

    // Liquid and gas phase properties
    double leftT = T(data->nlpts);
    double rightT = T(data->nlpts + 1);
    double leftR = liquidR(data->nlpts);
    double rightR = gasR(2);
    double delta_T = 10.00;

    double T_low = interfaceT ;
    double T_high = rightT ;

    // Initialize arrays and gas phase properties
    std::vector<double> MWArray(nsp), interfaceYArray_Old(nsp), gasDiffCoeffs(nsp), YArrayRight(nsp),YArrayAvg(nsp);
    for (size_t i = 0; i < nsp; ++i) {
        MWArray[i] = gas->molecularWeight(i);
        interfaceYArray_Old[i] = data->interfaceGasCellArr[2 + i];
        if (i < nsp) YArrayRight[i] = Y(data->nlpts + 1, i + 1);
        YArrayAvg[i] = HALF * (interfaceYArray_Old[i] + YArrayRight[i]);
    }
    gas->setState_TPY(HALF*(interfaceT+rightT), P, YArrayAvg.data());
    double lambdaGas = trmix->thermalConductivity();
    trmix->getMixDiffCoeffs(gasDiffCoeffs.data());
    double fuelDiffCoeff = gasDiffCoeffs[fuelIndex - 1];
    double rhoGas = gas->density();


    singleLiquidFuel->setLiquidState(HALF*(interfaceT+leftT), P);
    double lambdaLiquid = singleLiquidFuel->calculateThermalConductivity();

    if (data->flagSolveInterfaceProblem) {
        // New solver: enforce non-penetration for non-fuel species; fuel only crosses
        if(solveInterfaceNonPenetration(data, ydata,delta_t) == 0){
            // Keep liquid-side interface temperature consistent
            data->interfaceLiquidCellArr[1] = data->interfaceGasCellArr[1];
            return;
        }
        // Prepare for interface problem solving
        interfaceProblemPara nP = new interfaceProblemParameters;
//        auto nP = std:: <interfaceProblemParameters>();
        *nP = {nsp, fuelIndex, P, rhoGas, fuelDiffCoeff, lambdaGas, lambdaLiquid,
               liquidR(data->nlpts-1),leftR, dropletRadius, rightR, gasR(3),leftT, rightT,T(data->nlpts-1), T(data->nlpts + 2),
               interfaceYArray_Old.data(), YArrayRight.data(), MWArray.data()};

        // GSL root solver
        gsl_function F = {&interfaceIterationResidue, nP};
        gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // Evaluate the function at the endpoints
        double f_low = F.function(T_low, F.params);
        double f_high = F.function(T_high, F.params);

        // DEBUG
        printf("f_low  = %.6f .\n", f_low) ;
        printf("f_high = %.6f .\n", f_high) ;

        // --- NEW: Add checks for NaN and Infinity FIRST ---
        if (std::isnan(f_low) || std::isinf(f_low) || std::isnan(f_high) || std::isinf(f_high)) {
            printf("Function evaluated to NaN or Infinity at the boundaries. Cannot solve.\n");
            printf("f_low: %f, f_high: %f. Retaining old interface state.\n", f_low, f_high);
            
            // Clean up and skip the solver, just like in the other failure cases
            gsl_root_fsolver_free(solver);
            delete nP;

            double area = calc_area(dropletRadius, &data->metric);
            double mdot = rhoGas * (-fuelDiffCoeff) *
                          (interfaceYArray_Old[fuelIndex - 1] - YArrayRight[fuelIndex - 1]) /
                          (rightR - dropletRadius) / (1 - interfaceYArray_Old[fuelIndex - 1]) * area;
            data->dropletMass -= mdot * delta_t;

            return;
        }

        if (f_low * f_high > 0) {
            // If the endpoints do not straddle y = 0, retain old data and skip the solver
            printf("Brent solver endpoints do not straddle y=0. Retaining old interface state.\n");
            gsl_root_fsolver_free(solver);
            delete nP;

            // Skip solving iteration problem
            double area = calc_area(dropletRadius, &data->metric);
            double mdot = rhoGas * (-fuelDiffCoeff) *
                          (interfaceYArray_Old[fuelIndex - 1] - YArrayRight[fuelIndex - 1]) /
                          (rightR - dropletRadius) / (1 - interfaceYArray_Old[fuelIndex - 1]) * area;

            data->dropletMass -= mdot * delta_t;

            return;
        }

        // Set up and solve using the Brent solver
        gsl_root_fsolver_set(solver, &F, T_low, T_high);
        double root = 0.0;
        int iter = 0, status;
        const int max_iter = 20;

        do {
            gsl_root_fsolver_iterate(solver);
            root = gsl_root_fsolver_root(solver);
            status = gsl_root_test_interval(gsl_root_fsolver_x_lower(solver),
                                            gsl_root_fsolver_x_upper(solver),
                                            0, 1e-9);
        } while (status == GSL_CONTINUE && ++iter < max_iter);

        gsl_root_fsolver_free(solver);

        // Update based on vapor pressure
//        double vaporPressure = octaneVaporPressure(root);
        double vaporPressure = heptaneVaporPressure(root) ;
        double area = calc_area(dropletRadius, &data->metric);
        if (vaporPressure <= P) {
            std::vector<double> interfaceYArray_New(nsp);
            updateInterfaceMassFracArray(root, P, YArrayRight.data(),
                                         interfaceYArray_New.data(), MWArray.data(), nsp, fuelIndex);

            double mdot = rhoGas * (-fuelDiffCoeff) *
                          (interfaceYArray_New[fuelIndex - 1] - YArrayRight[fuelIndex - 1]) /
                          (rightR - dropletRadius) / (1 - interfaceYArray_New[fuelIndex - 1]) * area;

            data->dropletMass -= std::abs(mdot) * delta_t;
            data->interfaceGasCellArr[1] = root;
            data->interfaceLiquidCellArr[1] = root;
            std::copy(interfaceYArray_New.begin(), interfaceYArray_New.end(),
                      data->interfaceGasCellArr + 2);
        } else {
            data->flagSolveInterfaceProblem = false;
            double mdot = rhoGas * (-fuelDiffCoeff) *
                          (interfaceYArray_Old[fuelIndex - 1] - YArrayRight[fuelIndex - 1]) /
                          (rightR - dropletRadius) / (1 - interfaceYArray_Old[fuelIndex - 1]) * area;

            data->dropletMass -= mdot * delta_t;
        }
        delete nP;
    } else {
        // Skip solving iteration problem
        double area = calc_area(dropletRadius, &data->metric);
        double mdot = rhoGas * (-fuelDiffCoeff) *
                      (interfaceYArray_Old[fuelIndex - 1] - YArrayRight[fuelIndex - 1]) /
                      (rightR - dropletRadius) / (1 - interfaceYArray_Old[fuelIndex - 1]) * area;

        data->dropletMass -= mdot * delta_t;
    }
}

void updateDropletMass(double* ydata, UserData data, double delta_t) {
    // Extract reusable parameters
    size_t nsp = data->nsp;
    size_t fuelIndex = data->dropII - 1; // Convert to 0-based index
    double P = data->initialPressure * Cantera::OneAtm;
    double dropletRadius = interfaceR;
    double interfaceT = data->interfaceGasCellArr[1];

    // Initialize molecular weights and mass fractions
    std::vector<double>  interfaceYArray_Old(nsp), gasDiffCoeffs(nsp), YArrayRight(nsp);
    for (size_t i = 0; i < nsp; ++i) {
        interfaceYArray_Old[i] = data->interfaceGasCellArr[2 + i];
        YArrayRight[i] = Y(data->nlpts + 1, i + 1);
    }

    // Set gas state and extract properties
    gas->setState_TPY(interfaceT, P, interfaceYArray_Old.data());
    double lambdaGas = trmix->thermalConductivity();
    trmix->getMixDiffCoeffs(gasDiffCoeffs.data());
    double rhoGas = gas->density();
    double fuelDiffCoeff = gasDiffCoeffs[fuelIndex];

    // Calculate mass flux and update droplet mass
    double area = calc_area(dropletRadius, &data->metric);
    double mdot = rhoGas * (-fuelDiffCoeff) *
                  (interfaceYArray_Old[fuelIndex] - YArrayRight[fuelIndex]) /
                  (gasR(2) - dropletRadius) / (1 - interfaceYArray_Old[fuelIndex]) * area;
    data->dropletMass -= std::abs(mdot) * delta_t;
}

// Function to compute the derivative at a given x using Lagrange interpolation
double computeDerivative(const std::vector<double>& x, const std::vector<double>& y, double targetX) {
    // Number of points (3 in this case)
    int n = x.size();

    // Derivative at targetX
    double derivative = 0.0;

    // Loop to compute the derivative using Lagrange polynomial
    for (int i = 0; i < n; ++i) {
        double term = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue; // Skip the current point
            double product = 1.0;
            for (int k = 0; k < n; ++k) {
                if (k == i || k == j) continue; // Skip the current and j-th points
                product *= (targetX - x[k]) / (x[i] - x[k]);
            }
            term += product / (x[i] - x[j]);
        }
        derivative += y[i] * term;
    }

    return derivative;
}