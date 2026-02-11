#include "UserData.h"
#include "parse.h"

#ifndef GSL_DEF
#define GSL_DEF
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#endif

#include <cstring>
#include "omp.h"

thread_local Cantera::IdealGasMix* gas = nullptr;
thread_local Cantera::Transport* trmix = nullptr;
thread_local pSingleLiquidFuel singleLiquidFuel= nullptr;

void freeUserData(UserData data){
    // free global vars
#pragma omp parallel default(none)
    {
        if(trmix) {
            delete trmix;
            trmix = nullptr;
            printf("Transport on thread %d deleted!\n", omp_get_thread_num());
        }
        if(gas) {
            delete gas;
            gas = nullptr;
            printf("Gas on thread %d deleted!\n", omp_get_thread_num());
        }
        if(singleLiquidFuel){
            delete singleLiquidFuel;
            singleLiquidFuel= nullptr;
            printf("SingleLiquidFuel1 on thread %d deleted!\n", omp_get_thread_num());
        }
    }
	if(data!=NULL){
//		if(data->trmix!=NULL){
//			delete data->trmix;
//			printf("Transport Deleted!\n");
//
//		}
//		if(data->gas!=NULL){
//			delete data->gas;
//			printf("Gas Deleted!\n");
//
//		}
//        if(data->singleLiquidFuel!=NULL){
//            delete data->singleLiquidFuel;
//            printf("Single component liquid fuel structure deleted from Memory!\n");
//        }
//        if(data->constPressureReactor!=NULL){
//            delete data->constPressureReactor;
//            printf("Constant pressure reactor Deleted!\n");
//
//        }
//        if(data->reactorNet!=NULL){
//            delete data->reactorNet;
//            printf("Reactor net Deleted!\n");
//        }
		if(data->adaptiveGrid){
			if(data->grid->xOld!=NULL){
				delete[] data->grid->xOld;
				printf("old grid array Deleted!\n");
			}
			if(data->grid->x!=NULL){
				delete[] data->grid->x;
				printf("current grid array Deleted!\n");
			}
			if(data->grid!=NULL){
				free(data->grid);
				printf("grid object Freed!\n");
			}
		}
		else{
			if(data->uniformGrid!=NULL){
				delete[] data->uniformGrid;
				printf("uniformGrid deleted!\n");
			}
		}
		if(data->innerMassFractions!=NULL){
			delete[] data->innerMassFractions;
			printf("innerMassFractions array Deleted!\n");
		}
		if(data->output!=NULL){
  			fclose(data->output);
			printf("Output File Cleared from Memory!\n");
		}
		if(data->gridOutput!=NULL){
  			fclose(data->gridOutput);
			printf("Grid Output File Cleared from Memory!\n");
		}
		//if(data->ratesOutput!=NULL){
  		//	fclose(data->ratesOutput);
		//	printf("Rates Output File Cleared from Memory!\n");
		//}
		if(data->globalOutput!=NULL){
  			fclose(data->globalOutput);
			printf("Global Output File Cleared from Memory!\n");
		}
        if(data->dropletGlobalOutput!= NULL){
            fclose(data->dropletGlobalOutput);
            printf("Droplet Global Output File Cleared from Memory!\n");
        }

//        if(data->leftGhostCellArr!=NULL){
//            delete data->leftGhostCellArr;
//            printf("Left Ghost Cell Array deleted from Memory!\n");
//        }

        if(data->interfaceLiquidCellArr!=NULL){
            delete data->interfaceLiquidCellArr;
            printf("Interface Liquid Cell Array deleted from Memory!\n");
        }

        if(data->interfaceGasCellArr!=NULL){
            delete data->interfaceGasCellArr;
            printf("Interface Gas Cell Array deleted from Memory!\n");
        }

//        if(data->rightGhostCellArr!=NULL){
//            delete data->rightGhostCellArr;
//            printf("Right Ghost Cell Array deleted from Memory!\n");
//        }
//        if(data->heatSpline!=NULL){
//	        gsl_spline_free(data->heatSpline);
//            data->heatSpline=NULL;
//			printf("Heat Spline Cleared from Memory!\n");
//        }
//        if(data->heatAcc!=NULL){
//	        gsl_interp_accel_free(data->heatAcc);
//            data->heatAcc=NULL;
//			printf("Heat Acc Cleared from Memory!\n");
//        }
	}
  	free(data);             /* Free the user data */
	printf("\n\n");
}

UserData allocateUserData(FILE *input){

	UserData data;
  	data = (UserData) malloc(sizeof *data);

	if(!data){
		printf("Allocation Failed!\n");
		return(NULL);
	}
	setSaneDefaults(data);

	int ier;

	ier=parseNumber<size_t>(input, "gasBasePts"	, MAXBUFLEN, &data->npts);
	if(ier==-1 || data->npts<=0){
		printf("Enter non-zero basePts!\n");
		return(NULL);
	}

    ier=parseNumber<size_t>(input, "liquidBasePts"	, MAXBUFLEN, &data->nlpts);
    if(ier==-1 || data->nlpts<=0){
        printf("Enter non-zero liquidBasePts!\n");
        return(NULL);
    }

	ier=parseNumber<double>(input, "domainLength", MAXBUFLEN, &data->domainLength);
	if(ier==-1 || data->domainLength<=0.0e0){
		printf("domainLength error!\n");
		return(NULL);
	}

	ier=parseNumber<double>(input, "Rd", MAXBUFLEN, &data->Rd);
	if(ier==-1 || data->Rd<=0.0e0){
		printf("Rd error!\n");
		return(NULL);
	}
    data->initialRd = data->Rd ;

	ier=parseNumber<int>(input, "constantPressure"	, MAXBUFLEN, &data->constantPressure);
	if(ier==-1 || (data->constantPressure!=0 && data->constantPressure!=1)){
		printf("constantPressure error!\n");
		return(NULL);
	}

	ier=parseNumber<int>(input, "problemType"	, MAXBUFLEN, &data->problemType);
	if(ier==-1 || (data->problemType!=0 
		       && data->problemType!=1
			&& data->problemType!=2
			&& data->problemType!=3)){
		printf("include valid problemType!\n");
		printf("0: premixed combustion with NO equilibrated ignition kernel\n");
		printf("1: premixed combustion WITH equilibrated ignition kernel\n");
		printf("2: arbitrary initial conditions\n");
		printf("3: Restart\n");
		return(NULL);
	}

	ier=parseNumber<double>(input, "dPdt", MAXBUFLEN, &data->dPdt);

	ier=parseNumber<double>(input, "Rg", MAXBUFLEN, &data->Rg);
    if(data->Rg < 0.0){
        printf("Rg must be greater than 0");
        return(NULL);
    }

	ier=parseNumber<int>   (input, "reflectProblem"	, MAXBUFLEN, &data->reflectProblem);
	if(data->reflectProblem!=0 && data->reflectProblem!=1){
		printf("Invalid entry for reflectProblem! Can be only 1 or 0.\n");
		return(NULL);
	}
	
	ier=parseNumber<double>(input, "mdot"		, MAXBUFLEN, &data->mdot);
	
	ier=parseNumber<double>(input, "initialTemperature", MAXBUFLEN,
			&data->initialTemperature);
	if(ier==-1 || data->initialTemperature<=0.0e0){
		printf("Enter positive initialTemperature in K!\n");
		return(NULL);
	}


	ier=parseNumber<double>(input, "initialPressure", MAXBUFLEN,
			&data->initialPressure);
	if(ier==-1 || data->initialTemperature<=0.0e0){
		printf("Enter positive initialPressure in atm!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "metric"	, MAXBUFLEN, &data->metric);
	if(data->metric!=0 && data->metric!=1 && data->metric!=2){
		printf("Invalid entry for metric!\n");
		printf("0: Cartesian\n");
		printf("1: Cylindrical\n");
		printf("2: Spherical\n");
		return(NULL);
	}
    
	// ier=parseNumber<int>   (input, "heatType", MAXBUFLEN, &data->heatType);
	// if(data->heatType!=0 && data->heatType!=1){
	// 	printf("Invalid entry for heatType!\n");
	// 	printf("0: Constant heat source\n");
	// 	printf("1: Custom heat source\n");
	// 	return(NULL);
	// }

	ier=parseNumber<double>(input, "QDot", MAXBUFLEN, &data->maxQDot);
	if(ier==-1 && data->heatType==0){
		printf("Need to specify QDot for heatType 0!\n");
		return(NULL);
	}

	// ier=parseNumber<double>(input, "kernelSize", MAXBUFLEN, &data->kernelSize);
	// if(ier==-1){
	// 	printf("Need to specify kernelSize!\n");
	// 	return(NULL);
	// }

	ier=parseNumber<double>(input, "ignTime", MAXBUFLEN, &data->ignTime);
	if(ier==-1 && data->problemType==0){
		printf("Need to specify ignTime for problemType 0!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "rxn"	, MAXBUFLEN,
			&data->rxn);
	if(ier==-1 && data->rxn!=0 && data->rxn!=1){
		printf("rxn can either be 0 or 1!\n");
		return(NULL);
	}
    
    // if(data->heatType == 0){
	//     ier=parseNumber<double>(input, "ignTime", MAXBUFLEN, &data->ignTime);
	//     if(ier==-1){
	//     	printf("Need to specify ignTime for problemType 0!\n");
	//     	return(NULL);
	//     }
    // }else if(data->heatType == 1){
	//    	printf("heatType 1 selected, using last time in heatFile as ignTime\n");
    // }

	// ier=parseNumber<int>   (input, "electrodeLoss", MAXBUFLEN, &data->electrodeLoss);
	// if(data->electrodeLoss!=0 && data->electrodeLoss!=1){
	// 	printf("Invalid entry for electrodeLoss!\n");
	// 	printf("0: No electrode heat loss\n");
	// 	printf("1: Electrode heat loss\n");
	// 	return(NULL);
	// }

	// ier=parseNumber<double>(input, "eLossDelta", MAXBUFLEN,
	// 		&data->eLossDelta);
	// if(ier==-1 && data->electrodeLoss == 1){
	// 	printf("Need to specify eLossDelta when electodeLoss is being used!\n");
	// 	return(NULL);
	// }
   
	// ier=parseNumber<double>(input, "electrodeRadius", MAXBUFLEN,
	// 		&data->electrodeRadius);
	// if(ier==-1 && data->electrodeLoss == 1){
	// 	printf("Need to specify electrodeRadius when electodeLoss is being used!\n");
	// 	return(NULL);
	// }


	ier=parseNumber<double>(input, "mixingWidth", MAXBUFLEN,
			&data->mixingWidth);
	if(ier==-1){
		printf("Need to specify mixingWidth!\n");
		return(NULL);
	}
	
	ier=parseNumber<double>(input, "shift", MAXBUFLEN, &data->shift);

	ier=parseNumber<double>(input, "firstRadius", MAXBUFLEN, &data->firstRadius);
	
	ier=parseNumber<double>(input, "wallTemperature", MAXBUFLEN, &data->wallTemperature);

	ier=parseNumber<int>   (input, "dirichletInner"	, MAXBUFLEN,
			&data->dirichletInner);
	if(data->dirichletInner!=0 && data->dirichletInner!=1){
		printf("dirichletInner can either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "dirichletOuter"	, MAXBUFLEN,
			&data->dirichletOuter);
	if(data->dirichletOuter!=0 && data->dirichletOuter!=1){
		printf("dirichletOuter can either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "adaptiveGrid"	, MAXBUFLEN,
			&data->adaptiveGrid);
	if(ier==-1 || (data->adaptiveGrid!=0 && data->adaptiveGrid!=1)){
		printf("specify adaptiveGrid as 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "moveGrid"	, MAXBUFLEN,
			&data->moveGrid);
	if(ier==-1 || (data->moveGrid!=0 && data->moveGrid!=1)){
		printf("specify moveGrid as 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<double>   (input, "isotherm"	, MAXBUFLEN,
			&data->isotherm);
	if(ier==-1){
		printf("specify temperature of isotherm!\n");
		return(NULL);
	}

	ier=parseNumber<double>(input, "gridOffset", MAXBUFLEN, &data->gridOffset);

	ier=parseNumber<int>   (input, "nSaves"	, MAXBUFLEN, &data->nSaves);
	if(data->nSaves<0 ){
		printf("nSaves must be greater than 0!\n");
		return(NULL);
	}
	
	ier=parseNumber<int>   (input, "writeRates"	, MAXBUFLEN,
			&data->writeRates);
	if(data->writeRates!=0 && data->writeRates!=1){
		printf("writeRates must either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "writeEveryRegrid", MAXBUFLEN,
			&data->writeEveryRegrid);
	if(data->writeEveryRegrid!=0 && data->writeEveryRegrid!=1){
		printf("writeEveryRegrid must either be 0 or 1!\n");
		return(NULL);
	}

//	 ier=parseNumber<double>   (input, "writeDeltaT", MAXBUFLEN,
//	 		&data->writeDeltaT);

	ier=parseNumber<int>   (input, "setConstraints"	, MAXBUFLEN,
			&data->setConstraints);
	if(data->setConstraints!=0 && data->setConstraints!=1){
		printf("setConstraints must either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "suppressAlg"	, MAXBUFLEN,
			&data->suppressAlg);
	if(data->suppressAlg!=0 && data->suppressAlg!=1){
		printf("suppressAlg must either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "dryRun"	, MAXBUFLEN,
			&data->dryRun);
	if(data->dryRun!=0 && data->dryRun!=1){
		printf("dryRun must either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "printYdot"	, MAXBUFLEN,
			&data->printYdot);
	if(data->printYdot!=0 && data->printYdot!=1){
		printf("printYdot must either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<int>   (input, "printY"	, MAXBUFLEN,
			&data->printY);
	if(data->printY!=0 && data->printY!=1){
		printf("printY must either be 0 or 1!\n");
		return(NULL);
	}

	ier=parseNumber<double>   (input, "finalTime"	, MAXBUFLEN,
			&data->finalTime);
    ier=parseNumber<double>   (input, "deltaT"	, MAXBUFLEN,
                               &data->deltaT);
	ier=parseNumber<double>   (input, "relativeTolerance"	, MAXBUFLEN,
			&data->relativeTolerance);
	ier=parseNumber<double>   (input, "radiusTolerance"	, MAXBUFLEN,
			&data->radiusTolerance);
	ier=parseNumber<double>   (input, "temperatureTolerance", MAXBUFLEN,
			&data->temperatureTolerance);
	ier=parseNumber<double>   (input, "pressureTolerance", MAXBUFLEN,
			&data->pressureTolerance);
	ier=parseNumber<double>   (input, "massFractionTolerance", MAXBUFLEN,
			&data->massFractionTolerance);
	ier=parseNumber<double>   (input, "bathGasTolerance", MAXBUFLEN,
			&data->bathGasTolerance);
	// ier=parseNumber<double>   (input, "electronTolerance", MAXBUFLEN,
	// 		&data->electronTolerance);

	char gasFile[MAXBUFLEN],gasPhase[MAXBUFLEN],transportPhase[MAXBUFLEN],chem[MAXBUFLEN],dropSpecies[MAXBUFLEN],
         mix[MAXBUFLEN],tran[MAXBUFLEN],heat[MAXBUFLEN],cross[MAXBUFLEN],rad[MAXBUFLEN],therCond[MAXBUFLEN];

	// ier=parseNumber<char>(input, "gasFile"	, MAXBUFLEN, gasFile);
	// if(ier==-1){
	// 	printf("Enter gasFile!\n");
	// 	return(NULL);
	// }//else{
	 //  try{
	 //  	data->gas = new Cantera::IdealGasMix(chem);
	 //  	data->ngsp=data->gas->nSpecies(); //assign no: of species

    //		} catch (Cantera::CanteraError& err) {
	 //      printf("Error:\n");
	 //      printf("%s\n",err.what());
	 //      return(NULL);
    //		}
	 //}

// 	ier=parseNumber<char>(input, "gasPhase", MAXBUFLEN, gasPhase);
// 	if(ier==-1){
// 		printf("Enter gasPhase!\n");
// 		return(NULL);
// 	}else{
// 		try{
// 			data->gas = new Cantera::IdealGasMix(gasFile,gasPhase);
// 			data->nsp=data->gas->nSpecies(); //assign no: of species

//     		} catch (Cantera::CanteraError& err) {
// 		    printf("Error:\n");
// 		    printf("%s\n",err.what());
// 		    return(NULL);
//     		}
// 	}

// 	ier=parseNumber<char>(input, "transportPhase", MAXBUFLEN, transportPhase);
// 	if(ier==-1){
// 		printf("Enter transportPhase!\n");
// 		return(NULL);
// 	}else{
// 		try{
// 			data->transport = new Cantera::IdealGasMix(gasFile,transportPhase);

//     		} catch (Cantera::CanteraError& err) {
// 		    printf("Error:\n");
// 		    printf("%s\n",err.what());
// 		    return(NULL);
//     		}
// 	}
//    if(data->nsp-1 != data->transport->nSpecies()){
//       printf("The transport phase should contain all the species of the gas phase excluding the electron\n");
//    }

//     data->k_transport_ = (size_t*) malloc((data->nsp-1) * sizeof(size_t));
//     for(size_t tii = 0; tii < (data->nsp-1); tii++) {
//         std::string tSpecName = data->transport->speciesName(tii);
//         for(size_t gii = 0; gii < data->nsp; gii++) {
//             std::string gSpecName = data->gas->speciesName(gii);
//             if(tSpecName == gSpecName){
//                 data->k_transport_[tii] = gii;
//                 //printf("%d\n",gii);
//                 break;
//             }else if(gii == data->nsp-1){
// 		        printf("%s does not exist in combined file!\n",tSpecName);
// 		        return(NULL);
//             }
//         }
//     }

//	ier=parseNumber<char>(input, "chemistryFile"	, MAXBUFLEN, chem);
//	if(ier==-1){
//		printf("Enter chemistryFile!\n");
//		return(NULL);
//	}else{
//		try{
//            strcpy(data->model,chem);
//			data->gas = new Cantera::IdealGasMix(chem);
//			data->nsp=data->gas->nSpecies(); //assign no: of species
//
//    		} catch (Cantera::CanteraError& err) {
//		    printf("Error:\n");
//		    printf("%s\n",err.what());
//		    return(NULL);
//    		}
//	}

    ier=parseNumber<char>(input, "chemistryFile"	, MAXBUFLEN, chem);
    if(ier==-1){
        printf("Enter chemistryFile!\n");
        return(NULL);
    }else{
        strcpy(data->model,chem);
        int succeed = 0;
        size_t nsp_temp = 0;
#pragma omp parallel default(none) shared(chem, nsp_temp) reduction(+:succeed)
        {
            try {
                if (!gas) {
                    gas = new Cantera::IdealGasMix(chem);
                }
                #pragma omp single
                {
                    nsp_temp = gas->nSpecies();
                }
            } catch (Cantera::CanteraError &err) {
                printf("Error:\n");
                printf("%s\n", err.what());
                succeed +=1;
            }
        }
        if(succeed != 0) {
            return(NULL);
        }
        data->nsp = nsp_temp; //assign no: of species (outside parallel region)
    }



//    data->constPressureReactor = new Cantera::IdealGasConstPressureReactor() ;
//    if(data->constPressureReactor != NULL){
//        printf("Ideal gas constant pressure reactor object created! \n");
//    }else{
//        printf("Ideal gas constant pressure reactor object create failed! \n");
//        return(NULL);
//    }
//
//    data->reactorNet = new Cantera::ReactorNet() ;
//    if(data->reactorNet != NULL){
//        printf("Reactor net object created! \n");
//    }else{
//        printf("Reactor net object create failed! \n");
//        return(NULL);
//    }

	// ier=parseNumber<char>(input, "transportModel", MAXBUFLEN, tran);
	// if(ier==-1){
	// 	printf("Enter transportModel!\n");
	// 	return(NULL);
	// }else{
	// 	try{
  	// 		data->trmix = Cantera::newTransportMgr(tran, data->transport);
	// 	}catch (Cantera::CanteraError& err) {
	// 		printf("Error:\n");
	// 		printf("%s\n",err.what());
	// 		return(NULL);
    // 		}
	// }

//	ier=parseNumber<char>(input, "transportModel", MAXBUFLEN, tran);
//	if(ier==-1){
//		printf("Enter transportModel!\n");
//		return(NULL);
//	}else{
//		try{
//  			data->trmix = Cantera::newTransportMgr(tran, data->gas);
//		}catch (Cantera::CanteraError& err) {
//			printf("Error:\n");
//			printf("%s\n",err.what());
//			return(NULL);
//    		}
//	}

    ier=parseNumber<char>(input, "transportModel", MAXBUFLEN, tran);
    if(ier==-1){
        printf("Enter transportModel!\n");
        return(NULL);
    }else{
        int succeed = 0;
#pragma omp parallel default(none) shared(tran, chem) reduction(+:succeed)
        {
            try {
                if (!gas) {
                    gas = new Cantera::IdealGasMix(chem);
                }
                if (!trmix) {
                    trmix = Cantera::newTransportMgr(tran, gas);
                }
            } catch (Cantera::CanteraError &err) {
                printf("Error:\n");
                printf("%s\n", err.what());
                succeed ++;
            }
        }
        if(succeed != 0) {
            return (NULL);
        }
    }

//	ier=parseNumber<char>(input, "mixtureComposition", MAXBUFLEN, mix);
//	if(ier==-1){
//		printf("Enter mixtureComposition!\n");
//		return(NULL);
//	}else{
//		if(data->gas!=NULL){
//			try{
//				data->gas->setState_TPX(data->initialTemperature,
//							data->initialPressure*Cantera::OneAtm,
//							mix);
//			}catch (Cantera::CanteraError& err) {
//		    		printf("Error:\n");
//		    		printf("%s\n",err.what());
//		    		return(NULL);
//    			}
//		}
//	}

    ier=parseNumber<char>(input, "mixtureComposition", MAXBUFLEN, mix);
    if(ier==-1){
        printf("Enter mixtureComposition!\n");
        return(NULL);
    }else{
        int succeed = 0;
#pragma omp parallel default(none) shared(data, mix,Cantera::OneAtm) reduction(+:succeed)
        {
            if (gas != NULL) {
                try {
                    #pragma omp critical
                    {
                        gas->setState_TPX(data->initialTemperature,
                                        data->initialPressure * Cantera::OneAtm,
                                        mix);
                    }
                } catch (Cantera::CanteraError &err) {
                    printf("Error:\n");
                    printf("%s\n", err.what());
                    succeed++;
                }
            }
        }
        if(succeed != 0) {
            return NULL;
        }
    }


	ier=parseNumber<char>(input, "dropletComposition", MAXBUFLEN, dropSpecies);
	if(ier==-1){
		printf("Enter dropletComposition!\n");
		return(NULL);
	}else{
		if(gas!=NULL){
			try{
				data->dropII = dropletSpeciesIndex(data,dropSpecies);
			}catch (Cantera::CanteraError& err) {
		    		printf("Error:\n");
		    		printf("%s\n",err.what());
		    		return(NULL);
    			}
		}
	}


	// ier=parseNumber<char>(input, "heatFile", MAXBUFLEN, heat);
	// if(ier==-1){
	// 	printf("No heatFile entered!\n");
    //     if(data->heatType==1){
	// 	    printf("Enter heatFile when using heatType 1!\n");
	// 	    return(NULL);
    //     }
	// }else{
	// 	FILE* input;
	// 	if(input=fopen(heat,"r")){
    //         readHeatFile(input, data->heatSpline, data->heatAcc, &data->ignTime);
    //         //printf("before interp!\n");
   	//         //double test = gsl_spline_eval(data->heatSpline,1.0e-5,data->heatAcc);
    //         //printf("%0.5e: end interp!\n",test);
	// 		fclose(input);
	// 	}
	// 	else{
	// 		printf("heatFile file not found!\n");
	// 		return(NULL);
	// 	}
    // }
   
	// ier=parseNumber<char>(input, "crossFile", MAXBUFLEN, cross);
	// if(ier==-1){
	// 	printf("No crossFile entered!\n");
	//     return(NULL);
	// }else{
	// 	FILE* input;
	// 	if(input=fopen(cross,"r")){
    //         readCrossFile(input, data->crossSpline, data->crossAcc);
    //         //printf("before interp!\n");
   	//         //double test = gsl_spline_eval(data->crossSpline,1.0,data->crossAcc);
    //         //printf("%0.5e: end interp!\n",test);
	// 		fclose(input);
	// 	}
	// 	else{
	// 		printf("crossFile file not found!\n");
	// 		return(NULL);
	// 	}
    // }

	// ier=parseNumber<char>(input, "radFile", MAXBUFLEN, rad);
	// if(ier==-1){
	// 	printf("No radFile entered!\n");
	//     return(NULL);
	// }else{
	// 	FILE* input;
	// 	if(input=fopen(rad,"r")){
    //         readRadFile(input, data->radSpline, data->radAcc);
    //         //printf("before interp!\n");
   	//         //double test = gsl_spline_eval(data->radSpline,1.0188e4,data->radAcc);
    //         //printf("%0.5e: end interp!\n",test);
	// 		fclose(input);
	// 	}
	// 	else{
	// 		printf("radFile file not found!\n");
	// 		return(NULL);
	// 	}
    // }

    // ier=parseNumber<char>(input, "therCondFile", MAXBUFLEN, therCond);
	// if(ier==-1){
	// 	printf("No therCondFile entered!\n");
	//     return(NULL);
	// }else{
	// 	FILE* input;
	// 	if(input=fopen(therCond,"r")){
    //         readTherCondFile(input, data->therCondSpline, data->therCondAcc);
    //         //printf("before interp!\n");
   	//         //double test = gsl_spline_eval(data->radSpline,1.0188e4,data->radAcc);
    //         //printf("%0.5e: end interp!\n",test);
	// 		fclose(input);
	// 	}
	// 	else{
	// 		printf("therCondFile file not found!\n");
	// 		return(NULL);
	// 	}
    // }

	ier=parseNumber<int>   (input, "nThreads"	, MAXBUFLEN, &data->nThreads);
	if(data->nThreads<0 ){
		printf("nThreads must be greater than 0!\n");
		return(NULL);
	}

	data->k_bath=BathGasIndex(data);
	data->k_oxidizer=oxidizerIndex(data);
	data->k_OH=OHIndex(data);
	data->k_HO2=HO2Index(data);
	// data->k_e=eIndex(data);
	printf("Oxidizer index: %lu\n",data->k_oxidizer);
	printf("Bath gas index:%lu\n",data->k_bath);
	// printf("Electron index:%lu\n",data->k_e);

	//data->nr=0;
	data->nt=0;
	data->ny=data->nt+1;
	data->np=data->ny+data->nsp;
				
	data->nvar=data->nsp+2;		 //assign no: of variables (T,P,nsp species)
//    data->leftGhostCellArr = new double[data->nvar+1];  //"1" here represents the "R" (radius)
    data->interfaceLiquidCellArr = new double[data->nvar+1];
    data->interfaceGasCellArr = new double [data->nvar+1] ;
//    data->rightGhostCellArr = new double[data->nvar+1];

    //Allocate memory for spatial coordinate array
//    data->Rdata = new double[data->npts];
    data->Rdata = new double[data->nlpts+data->npts+1]; //"1": interface ghost cell


	if(!data->adaptiveGrid){
//		data->uniformGrid = new double[data->npts];
//		data->neq=data->nvar*data->npts;
        data->uniformGrid = new double[data->nlpts+data->npts+1]; //"1": interface ghost cell
        data->neq = data->nvar*(data->npts+data->nlpts);
	}
	else{
		data->grid=(UserGrid) malloc(sizeof *data->grid);
		ier=getGridSettings(input,data->grid);
		if(ier==-1)return(NULL);

		ier=initializeGrid(data->grid);
		if(ier==-1)return(NULL);

		ier=reGrid(data->grid, data->grid->position);
		if(ier==-1)return(NULL);

		data->npts=data->grid->nPts;
		data->neq=data->nvar*data->npts;
	}

	data->output=fopen("output.dat","w");
	data->globalOutput=fopen("globalOutput.dat","w");
	data->gridOutput=fopen("grid.dat","w");
    data->dropletGlobalOutput = fopen("dropletGlobalOuput.dat","w");

    //data->ratesOutput=fopen("rates.dat","w");
	
	data->innerMassFractions = new double [data->nsp];
//    data->singleLiquidFuel = new LiquidFuelProperties;
#pragma omp parallel default(none)
    {
        if (!singleLiquidFuel) {
            singleLiquidFuel = new LiquidFuelProperties;
        }
    }
	return(data);
}


void readHeatFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc, double* ignTime){

	size_t bufLen=10000;
	size_t nRows=0;
	size_t nColumns=2;

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

	double times[nRows], heat[nRows];
	size_t i=0;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
		}
		else{
			ret=strtok(buf," \t");
			times[i]=(double)(atof(ret));
			ret=strtok(NULL," \t");
			heat[i]=(double)(atof(ret));
         //DEBUG
         //printf("time: %15.6e, heat: %15.6e\n",times[i],heat[i]);
			i++;
		}
	}


    *ignTime = times[nRows-1];

	//gsl_interp_accel* acc;
	//gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, nRows);

	gsl_spline_init(spline,times,heat,nRows);


   	//ydata[j+k*nColumns]=gsl_spline_eval(spline,xNew[k],acc);
	
	//gsl_interp_accel_free(acc);
	//gsl_spline_free(spline);

}

void readCrossFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc){

	size_t bufLen=10000;
	size_t nRows=0;
	size_t nColumns=2;

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

	double energy[nRows], crossSec[nRows];
	size_t i=0;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
		}
		else{
			ret=strtok(buf," \t");
			energy[i]=(double)(atof(ret));
			ret=strtok(NULL," \t");
			crossSec[i]=(double)(atof(ret));
			i++;
		}
	}

	//gsl_interp_accel* acc;
	//gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, nRows);

	gsl_spline_init(spline,energy,crossSec,nRows);

   	//ydata[j+k*nColumns]=gsl_spline_eval(spline,xNew[k],acc);
	
	//gsl_interp_accel_free(acc);
	//gsl_spline_free(spline);

}

void readRadFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc){

	size_t bufLen=10000;
	size_t nRows=0;
	size_t nColumns=2;

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

	double temperature[nRows], radiation[nRows];
	size_t i=0;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
		}
		else{
			ret=strtok(buf," \t");
			temperature[i]=(double)(atof(ret));
			ret=strtok(NULL," \t");
			radiation[i]=(double)(atof(ret));
			i++;
		}
	}

    //*ignTime = times[nRows-1];

	//gsl_interp_accel* acc;
	//gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, nRows);

	gsl_spline_init(spline,temperature,radiation,nRows);

   	//ydata[j+k*nColumns]=gsl_spline_eval(spline,xNew[k],acc);
	
	//gsl_interp_accel_free(acc);
	//gsl_spline_free(spline);

}

void readTherCondFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc){

	size_t bufLen=10000;
	size_t nRows=0;
	size_t nColumns=2;

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

	double temperature[nRows], therCond[nRows];
	size_t i=0;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
		}
		else{
			ret=strtok(buf," \t");
			temperature[i]=(double)(atof(ret));
			ret=strtok(NULL," \t");
			therCond[i]=(double)(atof(ret));
			i++;
		}
	}

    //*ignTime = times[nRows-1];

	//gsl_interp_accel* acc;
	//gsl_spline* spline;
	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_steffen, nRows);

	gsl_spline_init(spline,temperature,therCond,nRows);

   	//ydata[j+k*nColumns]=gsl_spline_eval(spline,xNew[k],acc);
	
	//gsl_interp_accel_free(acc);
	//gsl_spline_free(spline);

}
void setSaneDefaults(UserData data){
	data->domainLength=1.0e-02;
	data->constantPressure=1;
	data->problemType=1;
	data->dPdt=0.0e0;
    data->Rg=1.0;
	data->reflectProblem=0;
	data->mdot=0.0;
	data->initialTemperature=300.0;
	data->initialPressure=1.0;
	data->metric=0;
    data->heatType=0;
	data->ignTime=1e-09;
	data->maxQDot=0.0e0;
	data->maxTemperature=300.0e0;
	data->mixingWidth=1e-04;
	data->shift=3.0e-03;
	data->firstRadius=1e-04;
	data->wallTemperature=298.0e0;
	data->dirichletInner=0;
	data->dirichletOuter=0;
	data->adaptiveGrid=0;
	data->moveGrid=0;
	data->gridOffset=0.0e0;
	data->isotherm=1000.0;
	data->nSaves=30;
	data->writeRates=0;
	data->writeEveryRegrid=0;
	data->writeDeltaT=1e-04;
	data->relativeTolerance=1e-06;
	data->radiusTolerance=1e-08;
	data->temperatureTolerance=1e-06;
	data->pressureTolerance=1e-06;
	data->massFractionTolerance=1e-09;
	data->bathGasTolerance=1e-06;
	data->finalTime=1e-02;
	data->tNow=0.0e0;
	data->setConstraints=0;
	data->suppressAlg=1;
	data->regrid=0;
	data->gridDirection=1;
	data->dryRun=0;
	data->nThreads=1;

	data->flamePosition[0]=0.0e0;
	data->flamePosition[1]=0.0e0;
	data->flameTime[0]=0.0e0;
	data->flameTime[1]=0.0e0;
	data->nTimeSteps=0;

	data->Mdot=0.0;
	data->Rd=2.0e-4;
	data->dropII=0;
	data->rxn=0;
    data->deltaT=1.0e-4;
    data->flagSolveInterfaceProblem= true;
}

//Locate bath gas:
size_t BathGasIndex(UserData data){
	size_t index=0;
	double max1;
	double max=gas->massFraction(0);
	for(size_t k=1;k<data->nsp;k++){
		max1=gas->massFraction(k);
		if(max1>=max){
			max=max1;
			index=k;
		}
	}
	return(index+1);
}

//Locate Oxidizer:
size_t oxidizerIndex(UserData data){
	size_t index=0;
	for(size_t k=1;k<=data->nsp;k++){
		if(gas->speciesName(k-1)=="O2"){
			index=k;
		}
	}
	return(index);
}

//Locate OH:
size_t OHIndex(UserData data){
	size_t index=0;
	for(size_t k=1;k<=data->nsp;k++){
		if(gas->speciesName(k-1)=="OH"){
			index=k;
		}
	}
	return(index);
}

//Locate HO2:
size_t HO2Index(UserData data){
	size_t index=0;
	for(size_t k=1;k<=data->nsp;k++){
		if(gas->speciesName(k-1)=="HO2"){
			index=k;
		}
	}
	return(index);
}

//Locate dropletFuelSpecies:
size_t dropletSpeciesIndex(UserData data, char* dropletSpec){
	size_t index=0;
	for(size_t k=1;k<=data->nsp;k++){
		if(gas->speciesName(k-1) == dropletSpec){
			index=k;
		}
	}
	//DEBUG
	printf("Droplet species index : %d  (function:allocateUserData)\n",index);
	return(index);
}

//Locate Electron:
size_t eIndex(UserData data){
	size_t index=0;
	for(size_t k=1;k<=data->nsp;k++){
		if(gas->speciesName(k-1)=="e-"){
			index=k;
		}
	}
	return(index);
}
