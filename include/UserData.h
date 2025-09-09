#ifndef CANTERA_DEF
#define CANTERA_DEF
#include <cantera/IdealGasMix.h>
#include <cantera/transport.h>
#include <cantera/zerodim.h>
#endif
#include <stdexcept>

#ifndef COOLPROP_DEF
#define COOLPROP_DEF
#include <CoolProp.h>
#include <CoolPropLib.h>
#endif

#include "gridRoutines.h"

#ifndef LIQUID_DEF
#define LIQUID_DEF
#include <iostream>
#include <stdexcept>
#include "CoolProp.h"  // Include CoolProp for property calculations

typedef struct LiquidFuelProperties {
    double temperature; // Temperature in Kelvin
    double pressure;    // Pressure in Pascal

    /**
     * @brief Default constructor initializes to a reasonable state.
     */
    LiquidFuelProperties()
            : temperature(298), pressure(1e5) {}  // Default to 25°C and 1 atm

    /**
     * @brief Set the state of the liquid with temperature and pressure.
     * @param T Temperature in Kelvin.
     * @param P Pressure in Pascal.
     */
    void setLiquidState(double T, double P) {
        if (T < 182 || T > 600) { // Valid range for n-heptane (boiling point ~ 371K)
            throw std::out_of_range("Temperature is out of valid range (182–600 K) for n-heptane.");
        }
        if (P < 1e5 || P > 1e7) { // Typical valid pressure range
            throw std::out_of_range("Pressure is out of valid range (1e5–1e7 Pa).");
        }
        temperature = T;
        pressure = P;
    }

    /**
     * @brief Calculate the density of liquid n-heptane.
     * @return Density in kg/m³.
     */
    double calculateDensity() const {
        return CoolProp::PropsSI("D", "T", temperature, "P", pressure, "nHeptane");
    }

    /**
     * @brief Calculate the thermal conductivity of liquid n-heptane.
     * @return Thermal conductivity in W/(m·K).
     */
    double calculateThermalConductivity() const {
        return CoolProp::PropsSI("L", "T", temperature, "P", pressure, "nHeptane");
    }

    /**
     * @brief Calculate the specific heat capacity of liquid n-heptane.
     * @return Specific heat capacity in J/(kg·K).
     */
    double calculateSpecificHeatCapacity() const {
        return CoolProp::PropsSI("C", "T", temperature, "P", pressure, "nHeptane");
    }

    /**
     * @brief Print all properties of n-heptane.
     */
    void print() const {
        try {
            std::cout << "n-Heptane Properties:\n";
            std::cout << "Temperature (K): " << temperature << "\n";
            std::cout << "Pressure (Pa): " << pressure << "\n";
            std::cout << "Density (kg/m³): " << calculateDensity() << "\n";
            std::cout << "Thermal Conductivity (W/(m·K)): " << calculateThermalConductivity() << "\n";
            std::cout << "Specific Heat Capacity (J/(kg·K)): " << calculateSpecificHeatCapacity() << "\n";
        } catch (const std::exception &e) {
            std::cerr << "Error while calculating properties: " << e.what() << '\n';
        }
    }
}* pSingleLiquidFuel;


#endif



#ifndef USER_DEF
#define USER_DEF
typedef struct UserDataTag{
	/*An ideal gas object from Cantera. Contains thermodynamic and kinetic
	 *info of all species.*/
//	Cantera::IdealGasMix* gas;
//   /*An ideal gas object that contains transport data for all species except the electron*/
//	Cantera::IdealGasMix* transport;
//	/*A Transport object from Cantera. Contains all transport info of all
//	 * species.*/
//	Cantera::Transport* trmix;
//    /*An isobaric 0D batch reactor */
//    Cantera::IdealGasConstPressureReactor* constPressureReactor;
//    /*A reactor net used to do time integration*/
//    Cantera::ReactorNet* reactorNet;
	/*Length of the domain (in meters):*/			  
	double domainLength;	  
	/*Mass of gas in domain (in kg):*/			  
	double mass;		  
	/*Parameter that indicates the symmetry of the problem;*/			  
	/*metric=0:Planar*/			  
	/*metric=1:Cylindrical*/			  
	/*metric=2:Spherical*/			  
	int metric;		  
	/*No: of species including both gas and plasma phase:*/
	size_t nsp;	  
	/*No: of gas phase species:*/
	size_t ngsp;	  
	/*No: of plasma phase species:*/
	size_t npsp;	  
	/*No: of equations:*/
	size_t neq;
	/*No: of variables:*/
	size_t nvar;
	/*Pointer indices (see "macros.h" for aliases that use these):*/
	/*Pointer index for temperature:*/
	size_t nt;	  
	/*Pointer index for species:*/
	size_t ny;
	/*Pointer index for spatial coordinate:*/
	//size_t nr;
   double* Rdata;
	/*Pointer index for pressure:*/
	size_t np;
	/*Species index of bath gas:*/			  
	size_t k_bath;		  
	/*Species index of oxidizer, OH, HO2, and the electron:*/			  
	size_t k_oxidizer;
	size_t k_OH;
	size_t k_HO2;
   size_t k_e;
	/*Map species associated with the transport phase (all species excluding the electron) 
    *to the phase including all species:*/			  
   size_t* k_transport_;
	/*User-defined mass flux (kg/m^2/s):*/	
	double mdot;
	/*Flag to solve isobaric/isochoric problem;*/
	/*constantPressure=1: isobaric*/
	/*constantPressure=0: isochoric*/
	int constantPressure;
	/*User-defined dPdt (Pa/s), activates when problem is "isobaric":*/			  
	double dPdt;	
	/*Initial temperature of the gas (K):*/
	double initialTemperature;
	/*Initial Pressure of the gas (atm):*/
	double initialPressure;
	/*Classification of problem type;*/
	/*problemType=0: Mixture is premixed and spatially uniform initially.
	 * In order for mixture to ignite, an external heat source (finite
	 * maxQDot) must be used.*/
	/*problemType=1: Mixture is premixed but spatially non-uniform
	 * initially. Equilibrium products are contained within a hot kernel of
	 * size given by "shift" and a mixing length scale given by
	 * "mixingWidth".*/
	/*problemType=2: User specified initial condition. Use file
	 * "initialCondition.dat".*/
	int problemType;

    /*Classification of Heat Source Type;*/ 
    /*heatType=0: Constant heat source whose magnitude is set by "Qdot"*/
    /*heatType=1: Custom heat source whose value is set by "dischargeFile"*/
    int heatType;
	/*Maximum External heat source (K/s):*/
	double maxQDot;
    
	/*Ignition kernel size:*/
	double kernelSize;

    // /*Objects to interpolate heat in time for custom heat source*/
	// gsl_interp_accel* heatAcc;
	// gsl_spline* heatSpline;

    // /*Objects to interpolate cross-sections of electron collisions with temperature*/
	// gsl_interp_accel* crossAcc;
	// gsl_spline* crossSpline;

	// /*Objects to interpolate radiation in temperature */
	// gsl_interp_accel* radAcc;
	// gsl_spline* radSpline;

	// /*Objects to interpolate thermal conductivity of air plasma with temperature*/
	// gsl_interp_accel* therCondAcc;
	// gsl_spline* therCondSpline;

    /*Set to 1 if electrode heat loss should be included, 0 otherwise*/
    int electrodeLoss;

    /*Length scale of electrode loss*/
    double eLossDelta;

    /*Length scale of electrode loss*/
    double electrodeRadius;

	double maxTemperature;
	/*Maximum time for which the external heat source is applied (s):*/
	double ignTime;
	/*Vector of Mass Fractions used to impose Robin Boundary Condition for
	 * species at the domain origin:*/
	double* innerMassFractions;
	/*Value of temperature to be used if Dirichlet Boundary Conditions are
	 * imposed for temperature:*/
	double innerTemperature;
	double wallTemperature;
	/*Isotherm chosen to find the location of a "burning" front (K):*/
	double isotherm;
	/*Interval of time integration:*/
	double finalTime;
	/*Current time:*/
	double tNow;
	/*Flag to reflect initial conditions across center of the domain:*/
	int reflectProblem;
	/*Parameters for initial conditions in setting up profiles:
	increasing function of x: g=0.5*(erf(x-3*w-shift)/w)+1)
	decreasing function of x: f=1-g*/
	double mixingWidth;
	double shift;
	double firstRadius;
	/*Flag to run program without time-integration i.e. simply layout the
	 * initial conditions and quit:*/
	int dryRun;
	/*Flag to run program and print ydot associated with the
	 * initial conditions and quit:*/
	int printYdot;
	/*Flag to run program and print y associated with the
	 * initial conditions and quit:*/
	int printY;
	/*Relative Tolerance:*/
	double relativeTolerance;
	/*Absolute Tolerance for spatial coordinate:*/
	double radiusTolerance;
	/*Absolute Tolerance for Temperature:*/
	double temperatureTolerance;
	/*Absolute Tolerance for Pressure:*/
	double pressureTolerance;
	/*Absolute Tolerance for Mass Fractions:*/
	double massFractionTolerance;
	/*Absolute Tolerance for bath gas mass fraction:*/
	double bathGasTolerance;
	/*Absolute Tolerance for electron mass fraction:*/
	double electronTolerance;
	/*Flag to set constraints on Mass fractions so they don't acquire
	 * negative values:*/
	int setConstraints;
	/*Flag to suppress error checking on algebraic variables:*/
	int suppressAlg;
	/*Number of time-steps elapsed before saving the solution:*/
	int nSaves;		  
	/*Flag to set write for every regrid:*/
	int writeEveryRegrid;
	double writeDeltaT;
	/*Solution output file:*/
	FILE* output;
	/*Flag to write the rates (ydot) of solution components into the
	 * "ratesOutput" file:*/
	int writeRates;
	/*Grid output file:*/
	FILE* gridOutput;
	///*Rate of change (ydot) output file (see "writeRates"):*/
	//FILE* ratesOutput;
	/*Global properties (mdot, radius of flame, etc.) output file:*/
	FILE* globalOutput;

	/*Flag to adapt grid:*/
	int adaptiveGrid;
	/*Flag to move grid:*/
	int moveGrid;
	/*Flag to initiate regrid:*/
	int regrid;
	/*Integer that specifies grid movement direction:
	 * gridDirection = -1: Move Left
	 * gridDirection = +1: Move Right*/ 
	int gridDirection;

	/*Grid Ratio: This replaces the uniform grid. dX0 and dXf are the grid initial
	and final grid spacing, respectivly. The Grid Ratio (Rg) is equal to dXf/dX0. A
	Rg>1 focuses grid points on the droplet surface while a Rg<1 focuses grid
	points at the right boundary. A Rg of 1 is a uniform grid.*/ 
    double Rg;

	/*Total number of points for grid:*/
	size_t npts;

	double gridOffset;

	UserGrid grid;
	double* uniformGrid;

	int dirichletInner,dirichletOuter;

	int nThreads;
	double clockStart;
	
	/*These arrays are used to compute dr/dt, which in turn is used to
	 * compute the flame speed S_u:*/
	double flamePosition[2];
	double flameTime[2];
	size_t nTimeSteps;

	//droplet relevant parameters 
	double Mdot;
	
	/*initial droplet radius*/
	double Rd;
	size_t dropII;

	/*switch for rxn term*/
	int rxn;
    double deltaT;

    char model[MAXBUFLEN];
    size_t nlpts;
//    double* leftGhostCellArr;
    double* interfaceLiquidCellArr;
    double* interfaceGasCellArr;  // arrays only keep state variables
//    double* rightGhostCellArr;
    double dropletMass;

//    pSingleLiquidFuel singleLiquidFuel;
    bool flagSolveInterfaceProblem;
    FILE* dropletGlobalOutput;
    double initialRd;
} *UserData;


UserData allocateUserData(FILE *input);
void setSaneDefaults(UserData data);
void freeUserData(UserData data);
void readHeatFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc, double* ignTime);
void readCrossFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc);
void readRadFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc);
void readTherCondFile(FILE* input, gsl_spline*& spline, gsl_interp_accel*& acc);

size_t BathGasIndex(UserData data);
size_t oxidizerIndex(UserData data);
size_t OHIndex(UserData data);
size_t HO2Index(UserData data);
size_t dropletSpeciesIndex(UserData data, char* dropletSpec);
size_t eIndex(UserData data);
#endif

#ifndef INTERFACE_DEF
#define INTERFACE_DEF
typedef struct interfaceProblemParameters {
    size_t nsp , fuelIndex;
    double P,gasDensity;
    double gasDiffCoeff;
    double lambda_gas , lambda_liquid;
    double left_R1,left_R , interface_R, right_R,right_R1;
    double left_T , right_T;
    double left_T1,right_T1;
    double* YArrayMiddle_Old;
    double* YArray_Right;
    double* MWArray;
}* interfaceProblemPara;

#endif

extern thread_local Cantera::IdealGasMix* gas;
extern thread_local Cantera::Transport* trmix;
extern thread_local pSingleLiquidFuel singleLiquidFuel;
