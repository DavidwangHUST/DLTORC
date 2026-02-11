/*
 _____ ___  ____   ____ 
|_   _/ _ \|  _ \ / ___|
  | || | | | |_) | |    
  | || |_| |  _ <| |___ 
  |_| \___/|_| \_\\____|
                        
*/

#include "UserData.h"
#include "solution.h"
#include "residue.h"
#include "macros.h"
#include "timing.h"

//DEBUG
//#include <math.h>

//#include <ida/ida.h>
//#include <ida/ida_direct.h>
//#include <ida/ida_spils.h>
//#include <sunmatrix/sunmatrix_band.h>
//#include <sunlinsol/sunlinsol_lapackband.h>
//#include <sunlinsol/sunlinsol_spgmr.h>
//#include <sunlinsol/sunlinsol_band.h>
//#include <ida/ida_band.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts. */
#include <kinsol/kinsol.h>		
#include <nvector/nvector_openmp.h>    /* serial N_Vector types, fcts., macros */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix */
#include <sunlinsol/sunlinsol_lapackband.h>
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver */
#include <cvode/cvode_direct.h>        /* access to CVDls interface */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>    /* definition of ABS and EXP */
#include <ctime>
#include <omp.h>


static int check_flag(void *flagvalue, 
		      const char *funcname, 
		      int opt);

void freeAtLast(void* memCVODE, N_Vector *y,
		N_Vector *ydot, 
		N_Vector *atolv,
		N_Vector *constraints,UserData data);

N_Vector make_vec(size_t n);

int main(){

  	FILE *input;input=fopen("input.dat","r");
  	UserData data;
      data=NULL;
      data=allocateUserData(input);
	fclose(input);
	data->clockStart=get_wall_time();

	if(data==NULL){
		printf("check input file!\n");
  		freeUserData(data);
		return(-1);
	}

    omp_set_dynamic(0);
    if (data->nThreads > 0) {
        omp_set_num_threads(data->nThreads);
    }
	long int ier,mu,ml,count,netf,ncfn,njevals,nrevals;
	realtype tNow,*atolvdata,*constraintsdata,finalTime,tolsfac;

  	N_Vector y,ydot,atolv,constraints;
  	y=ydot=atolv=constraints=NULL;
 	ier=allocateSolution(data->neq,data->nThreads,&y,&ydot,&atolv,&constraints);
	ier=setInitialCondition(&y,&ydot,data);
	if(ier==-1)return(-1);
	tNow=data->tNow;
	finalTime=data->finalTime;

	double* ydata;
	double* ydotdata;
  	ydata    = N_VGetArrayPointer_OpenMP(y);
  	ydotdata = N_VGetArrayPointer_OpenMP(ydot);
      double deltaT = data->deltaT;

//	void *mem;mem=NULL;mem = CVodeCreate(CV_BDF, CV_NEWTON);
//    void *mem;mem=NULL;mem = CVodeCreate(CV_BDF,CV_NEWTON);
    void *memCVODE;memCVODE=NULL;memCVODE = CVodeCreate(CV_BDF,CV_NEWTON);
	// void *memCVODE;mem=NULL;mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);

   ier = CVodeSetUserData(memCVODE, data);
   if(check_flag(&ier, "CVodeSetUserData", 1)) return(1);

//   ier = CVodeInit(memCVODE, fun, tNow, y);
    ier = CVodeInit(memCVODE,funNew,tNow,y);
   if(check_flag(&ier, "CVodeInit", 1)) return(1);
   // Atol array

   atolvdata = N_VGetArrayPointer_OpenMP(atolv);
   for (size_t i = 1; i <=(data->nlpts+data->npts); i++) {
     atolT(i) = data->temperatureTolerance;
     for (size_t k = 1; k <=data->nsp; k++) {
		   if(k!=data->k_bath){
                		atolY(i,k) = data->massFractionTolerance;
			}
			else if(k == data->k_bath){
                		atolY(i,k) = data->bathGasTolerance;
         }
			// else if(k == data->k_e){
            //     		atolY(i,k) = data->electronTolerance;
			// }

     }
     //atolR(i) = data->radiusTolerance;
     atolP(i) = data->pressureTolerance;
	}
   //ier = CVodeSStolerances(mem, data->relativeTolerance, atolv);
   ier = CVodeSVtolerances(memCVODE, data->relativeTolerance, atolv);
   if (check_flag(&ier, "CVodeSStolerances", 1)) return(1);

  	mu = 2*data->nvar; ml = mu;
	SUNMatrix A; A=NULL;
	A=SUNBandMatrix(data->neq,mu,ml,mu+ml);
//    A = SUNBandMatrix(data->neq,mu,ml);
	SUNLinearSolver LS; LS=NULL;
	//LS=SUNBandLinearSolver(y,A);
	LS=SUNLapackBand(y,A);

//--------------------------------create KINSOL solver--------------------------------
	const size_t n_if = 2 + data->nsp -1 ; // Ti,mdot , Y1...Y_{nsp-1}
	N_Vector u_iface = make_vec(n_if);
	N_Vector s_scale = make_vec(n_if);
    N_Vector f_scale = make_vec(n_if);
    
	NV_Ith_S(u_iface,0) = data->interfaceGasCellArr[1];
	NV_Ith_S(u_iface,1) = data->Mdot;
	for(size_t k =0 ;k<data->nsp-1;k++) NV_Ith_S(u_iface,k+2) = data->interfaceGasCellArr[k+2];
	N_VConst(1.0,s_scale);
    N_VConst(1.0e-3,f_scale);



//	printf("CVODE memory: %p\n", memCVODE);
//	printf("SUNMatrix A: %p\n", A);
//	printf("SUNLinearSolver LS: %p\n", LS);

	ier=CVDlsSetLinearSolver(memCVODE,LS,A);
//    ier = CVodeSetLinearSolver(memCVODE,LS,A);

   ier=CVodeSetStabLimDet(memCVODE, SUNTRUE);

   //constraintsdata = N_VGetArrayPointer_OpenMP(constraints);
	//if(data->setConstraints){
	//	for (size_t i = 1; i <=data->npts; i++) {
   //          		for (size_t k = 1; k <=data->nsp; k++) {
   //          		   constraintsY(i,k) = ONE;
   //          		}
	//	}
  	//	ier=IDASetConstraints(mem, constraints);
	//}

    if(data->printYdot){
        for(size_t ii = 0; ii < data->neq; ii++){
            ydotdata[ii] = 0.0;
        }

        fun(0.0, y, ydot, data);

        printf("Printing ydot.dat\n");
        FILE* ydotOutput = fopen("ydot.dat","w");
        printSpaceTimeHeader(data, ydotOutput);
        printSpaceTimeOutput(0.0, &ydot, ydotOutput, data);
        fclose(ydotOutput);
        printf("done!\n");
        return(0);
    }

    if(data->printY){
        printf("Printing y.dat\n");
        FILE* yOutput = fopen("y.dat","w");
        printSpaceTimeHeader(data, yOutput);
        printSpaceTimeOutput(0.0, &y, yOutput, data);
        fclose(yOutput);
        printf("done!\n");
        return(0);
    }
	
	printSpaceTimeHeader(data,data->output);
    printDropletGlobalHeader(data,data->dropletGlobalOutput);
//	printGlobalHeader(data);
	printSpaceTimeOutput(tNow, &y, data->output, data);
	printSpaceTimeOutput(tNow, &y, data->gridOutput, data);
    printDropletGlobalOutput(data,data->dropletGlobalOutput,tNow);

	if(!data->dryRun){
		count=0;
		double dt=1e-08;
		double t1=0.0e0;
		double xOld=0.0e0;
		double x=0.0e0;
		double dx=0.0e0;
		double dxMin=1.0e0;
		double dxRatio=dx/dxMin;
		double writeTime=tNow;
        double maxT;
		int move=0;
		int kcur=0;
		if(data->adaptiveGrid){
			dxMin=data->grid->leastMove;
			//xOld=maxCurvPosition(ydata, data->nt, data->nvar,
			//		     data->grid->x, data->npts);
			xOld=isothermPosition(ydata, data->isotherm, data->nt, 
	                		   data->nvar, data->grid->x, data->npts);
		}

      //DEBUG
      //double maxJac = SM_ELEMENT_B(A,0,0);
      //double curJac;
      //size_t maxII = 0, maxJJ = 0;

        while (tNow <= finalTime && data->Rd >= (data->initialRd * 0.1)) {
            maxT = maxTemperature(ydata,data);
            printf("Max temperature in the domain : %.3f [K] \n",maxT);
			t1=tNow;
//            isobaricAdvance(ydata,data) ;
//            isobaricAdvance2(ydata,data);

            //Compute ydot
//            fun(tNow, y, ydot, data);
            //Compute spatial coordinates
//            getR(y, data);
            getRNew(ydata,data) ;
            clock_t start_1,end_1;
            start_1 = clock();

//         ier = CVode(memCVODE, finalTime, y, &tNow, CV_ONE_STEP);
            ier = CVode(memCVODE, finalTime, y, &tNow, CV_ONE_STEP);

            end_1 = clock();
  			if(check_flag(&ier, "CVode", 1)){
				freeAtLast(memCVODE,&y,&ydot,&atolv,&constraints,data);
		       		return(-1);
			}
            printf("Time elapsed for CVode: %f\n",(double)(end_1-start_1)/CLOCKS_PER_SEC);

            dt=tNow-t1;

            //TODO: update the interface gas and liquid array (except R), mdot and dropletMass
//            updateInterfaceState(ydata,data,dt) ;
            updateInterfaceCell(ydata,data,dt) ;

            //TODO: only update the droplet mass and radius
//            updateDropletMass(ydata,data,dt) ;

            // DEBUG
            printf("Interface temperature : %.3f [K] \n",data->interfaceGasCellArr[1]);
            printf("Droplet radius : %.3e [m] \n",data->interfaceGasCellArr[0]);
            printDropletGlobalOutput(data,data->dropletGlobalOutput,tNow);

            ier = CVodeGetCurrentOrder(memCVODE, &kcur);
//            isobaricAdvance(ydata,data);
//            isobaricAdvance2(ydata,data);


         //DEBUG
         //maxJac = SM_ELEMENT_B(A,0,0);
         //curJac = maxJac;
         //maxII = 0;
         //maxJJ = 0;
         //for(int JJ = 0; JJ < data->neq; JJ++){
         //   for(int II = std::max(0,JJ - (int)mu); II <= std::min((int)data->neq-1, JJ+(int)ml); II++){
         //      curJac = SM_ELEMENT_B(A,II,JJ);
         //      if(abs(curJac) > maxJac){
         //         maxJac = curJac;
         //         maxII = II;
         //         maxJJ = JJ;
         //      }
         //   }
         //}
         //printf("maxJac=%0.15e, II=%d, JJ=%d\n",maxJac,maxII,maxJJ);

         //DEBUG
         //if(dt < 1e-11 && tNow > 1e-9){
         //   FILE* jacFile = fopen("jac.dat","w");
         //   SUNBandMatrix_Print(A,jacFile);
         //   fclose(jacFile);
         //   return(-1);
         //}




			//if(data->adaptiveGrid==1 && data->moveGrid==1){
			//	//x=maxCurvPosition(ydata, data->nt, data->nvar,
			//	//		  data->grid->x, data->npts);

			//	x=isothermPosition(ydata, data->isotherm, data->nt, 
	      //          			   data->nvar, data->grid->x, data->npts);

			//	//x=maxGradPosition(ydata, data->nt, data->nvar,
			//	//		  data->grid->x, data->npts);
			//	dx=x-xOld;

			//	if(dx*dxMin>0.0e0){
			//		move=1;
			//	}else{
			//		move=-1;
			//	}

			//	//if(fabs(dx)>=dxMin && x+(double)(move)*0.5e0*dxMin<=1.0e0){
			//	dxRatio=fabs(dx/dxMin);
			//	if(dxRatio>=1.0e0 && dxRatio<=2.0e0){
			//		printf("Regridding!\n");

			//		data->regrid=1;
			//		//printSpaceTimeOutput(tNow, &y, data->gridOutput, data);

			//		ier=reGrid(data->grid, x+(double)(move)*0.5e0*dxMin);
			//		if(ier==-1){
			//			freeAtLast(mem,&y,&ydot,&atolv,&constraints,data);
			//			return(-1);
			//		}

			//		updateSolution(ydata, ydotdata, data->nvar,
			//		               data->grid->xOld,data->grid->x,data->npts);
			//		storeGrid(data->grid->x,data->grid->xOld,data->npts);
			//		xOld=x;

			//		printf("Regrid Complete! Restarting Problem at %15.6e s\n",tNow);
			//		ier = IDAReInit(mem, tNow, y, ydot);
  			//		if(check_flag(&ier, "IDAReInit", 1)){
			//			freeAtLast(mem,&y,&ydot,&atolv,&constraints,data);
		   //    				return(-1);
			//		}
			//		ier = IDASetInitStep(mem,1e-01*dt);
			//		printf("Reinitialized! Calculating Initial Conditions:\n");
			//		printf("Cross your fingers...\n");
  			//		ier = IDACalcIC(mem, IDA_YA_YDP_INIT, tNow+1e-01*dt);
  			//		if(check_flag(&ier, "IDACalcIC", 1)){
  			//			ier = IDACalcIC(mem, IDA_YA_YDP_INIT, tNow+1e+01*dt);
			//		}
			//		//Every once in a while, for reasons
			//		//that befuddle this mathematically
			//		//lesser author, IDACalcIC fails. Here,
			//		//I desperately try to make it converge
			//		//again by sampling increasingly larger
			//		//time-steps:
			//		for (int i = 0; i < 10; i++) {
  			//			ier = IDACalcIC(mem, IDA_YA_YDP_INIT, tNow+(1e-01+pow(10,i)*dt));
			//			if(ier==0){
			//				break;
			//			}
			//		}
			//		//Failure :( Back to the drawing board:
  			//		if(check_flag(&ier, "IDACalcIC", 1)){
			//			freeAtLast(mem,&y,&ydot,&id,&atolv,&constraints,data);
		   //    				return(-1);
			//		}
			//		printf("Initial (Consistent) Conditions Calculated!\n");
			//		printf("------------------------------------------\n\n");
			//		//if(data->writeEveryRegrid){

			//	}
			//}

            if(count%data->nSaves==0 && !data->writeEveryRegrid){
                printSpaceTimeOutputNew(tNow,ydata,data->output,data);
//				printSpaceTimeOutput(tNow, &y, data->output, data);
				//if(data->writeRates){
				//	printSpaceTimeRates(tNow, ydot, data);
				//}
			}

			/*Print global variables only if time-step is of high order!*/
			if(data->nTimeSteps==0){
				data->flamePosition[0]=0.0e0;
				data->flamePosition[1]=0.0e0;
				data->flameTime[0]=tNow;
				data->flameTime[1]=tNow;
			}
//			if(kcur>=2 ){
//			      printGlobalVariables(tNow,&y,&ydot,data);
//			}

//			if(tNow>writeTime){
//				writeTime+=data->writeDeltaT;
//				printSpaceTimeOutput(tNow, &y, data->output, data);
//				FILE* fp;
//				fp=fopen("restart.bin","w");
//				writeRestart(tNow, &y, fp, data);
//				fclose(fp);
//			}

            //TODO: Redo this for CVODE
            ier = CVodeGetNumErrTestFails(memCVODE, &netf);
            ier = CVodeGetNumNonlinSolvConvFails(memCVODE, &ncfn);
            ier = CVDlsGetNumJacEvals(memCVODE, &njevals);
            ier = CVodeGetNumRhsEvals(memCVODE, &nrevals);
            printf("etf = %ld ,"
                   "nlsf= %ld ,"
                   "J=%ld ,"
                   "R=%ld ,"
                   "o=%d ,",netf, ncfn, njevals, nrevals, kcur);
            printf("Time=%15.6e s,",tNow);
            printf("dt=%15.6e s,",dt);
            printf("frac: %15.6e\n",dxRatio);

			count++;
			data->nTimeSteps=count;

			printf("\nElapsed Wall Clock Time: %15.6e s\n",get_wall_time()-data->clockStart);
     		printf("Elapsed CPU Time:        %15.6e s\n",get_cpu_time());
			printf("---------------------------------------------------------\n\n\n");
		}
	}



	SUNLinSolFree(LS);
	SUNMatDestroy(A);
	freeAtLast(memCVODE,&y,&ydot,&atolv,&constraints,data);

	return(0);
}

void freeAtLast(void* memCVODE,
		N_Vector *y, 
		N_Vector *ydot, 
		N_Vector *atolv,
		N_Vector *constraints,UserData data){

  	CVodeFree(&memCVODE);
	freeSolution(y,ydot,atolv,constraints);
  	freeUserData(data);
}

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); 
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}

N_Vector make_vec(size_t n){return N_VNew_Serial(static_cast<long>(n));}