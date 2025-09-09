//
// Created by weiye on 11/19/24.
//

#ifndef SINGLEPHASELTORC_NEWFORSOLVER_H
#define SINGLEPHASELTORC_NEWFORSOLVER_H

#include <cantera/IdealGasMix.h>
#include <cantera/transport.h>
#include <cantera/zerodim.h>

//struct idealGas {
//    void *obj;
//};
//
//struct transport {
//    void *obj;
//};
//
//struct isochoricReactor {
//    void *obj;
//};

struct isobaricReactor {
    void *obj;
};

struct reactorNet {
    void *obj;
};

//typedef struct idealGas idealGas_t;
//typedef struct transport transport_t;
//typedef struct isochoricReactor isochoricReactor_t;
typedef struct isobaricReactor isobaricReactor_t;
typedef struct reactorNet reactorNet_t;

//idealGas_t *idealGasCreate(const char* mech);
//transport_t *transportCreate(idealGas_t *gas);
//isochoricReactor_t *isochoricReactorCreate();
isobaricReactor_t *isobaricReactorCreate();
reactorNet_t *reactorNetCreate();

//void idealGasDestroy(idealGas_t *gas);
//void isochoricReactorDestroy(isochoricReactor_t *r);
void isobaricReactorDestroy(isobaricReactor_t *r);
void reactorNetDestroy(reactorNet_t *sim);

//void addGasToIsochoricReactor(idealGas_t *gas, isochoricReactor_t *r);
//void addGasToIsobaricReactor(idealGas_t *gas, isobaricReactor_t *r);
//void setIsochoricReactorState(isochoricReactor_t *r, double* y);
void setIsobaricReactorState(isobaricReactor_t *r, double* y);
//void getIsochoricReactorState(isochoricReactor_t *r, double* y);
void getIsobaricReactorState(isobaricReactor_t *r, double* y);
//void addIsochoricReactorToNet(isochoricReactor_t *r, reactorNet_t *sim);
void addIsobaricReactorToNet(isobaricReactor_t *r, reactorNet_t *sim);

void reactorNetAdvance(reactorNet_t *sim, double t);
//void reactorNetReinitialize(reactorNet_t *sim);
void reactorNetSetTolerances(reactorNet_t *sim, double rtol, double atol);
void setNetTime(reactorNet_t *sim, double t);

/* modified version for LTORC */
void addGasToIsobaricReactor(Cantera::IdealGasMix *gas, isobaricReactor_t *r);


#endif //SINGLEPHASELTORC_NEWFORSOLVER_H
