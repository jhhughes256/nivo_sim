# Nivolumab Population Pharmacokinetic Model (20170803_Nivolumab18)

# ------------------------------------------------------------------------------
# Define the model parameters and equations

library(dplyr)	    #New plyr - required for mrgsolve
library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics


code <- '

$INIT    // Initial Conditions for Compartments

CMT1      =  0,         // Central Compartment
CMT2      =  0,         // Peripheral Compartment
CMT3      =  60,         // Tumor
CMT4      =  0,         // Survival
CMT5      =  0,         // Censor
CMT6      =  0,        // Area under the curve

$SET      // Set Differential Equation Solver Options			
atol          = 1e-8, rtol = 1e-8
maxsteps      = 100000

$PARAM   
// Population parameters

AVCLH     = 0.06,      // clearance(L/h) at health status
VC        = 5.0,       // Central Volume (L)
Q         = 0.5,       // Intercompartmental CL (L/h)
VP        = 5.0,       // Peripheral Volume (L)

TVEMAX    = 0.02,
TVEC50    = 20,
TVTG        = 0.005,    // Tumor gowth rate
TVLAMDA     = 0.01     // Resistance

LAMBS     = 0.0001,     //Scale parameter in the Weibull probability density function for the survival model
ALPHS     = 1,          //Shape parameter in the Weibull probability density function for the survival model
LAMBC     = 0.001,      //Scale parameter in the Weibull probability density function for the drop out model
ALPHC     = 1,          //Shape parameter in the Weibull probability density function for the drop out model

TUMLIM = 1000

//Population values(allocated in population)
ECOG      = 1,
UCENSOR   = 1,
UEVENT    = 1,


// Covariate Effects
TUMSLD_CL = 0.3,       // Power of TUMSLD effect on CL
ECOG_CL   = 1.5,       // ECOG effect on CL
TSHAZ     = 0.02,       // Tumor size on hazard
ECOGHAZ   = 1,          // ECOG on hazard

// Default ETA Values for Simulation (allocated in creation of population)
ETA1      =  0,         // ZCL
ETA2      =  0,         // ZVC
ETA3      =  0,        // ZVP
ETA4      =  0,        // ZEMAX
ETA5      =  0,        // ZEC50
ETA6      =  0,        // ZTG
ETA7      =  0,        // ZR
ETA8      =  0,        // ZHZ


$OMEGA    // Population Parameter Variability
name  = "omega"
block = FALSE
labels  = s(ZCL, ZVC, ZVP, ZEMAX, ZEC50, ZTG, ZR, ZHZ)
0.1 //ZCL
0.1 //ZVC
0.1 //ZVP
1 //ZEMAX
0.01 //ZEC50
0.1 //ZTG
0.5 //ZR
0.1 //ZHZ



$MAIN // Individual Parameter Values

                      //---DRUG EXPOSURE---
double DELT = 1 ;
double CLH = AVCLH * exp(ETA1) ;
double TUM = CMT3 + DELT ;
double CL=CLH * pow(ECOG_CL,ECOG) * pow(TUM,TUMSLD_CL);

double V1 = VC * exp(ETA2) ;
double V2= VP * exp(ETA3) ;

double K10=CL/V1 ;
double K12=Q/V1 ;
double K21=Q/V2 ;
double S1=V1 ;
                     //---TUMOR GROWTH---
double EMAX   =   TVEMAX * exp(ETA4) ;
double EC50   =   TVEC50 * exp(ETA5) ;

double TG = TVTG * exp(ETA6);
double LAMDA = TVLAMDA* exp(ETA7);

double CT  = CMT1/V1 ;

                     //---SURVIVAL WEIBULL---


double IIVHAZ   =   (ETA8);


$ODE

                      //---DRUG EXPOSURE---
double DEL= pow(10,-6);
double C1 =  CMT1/V1;
double C2 =  CMT2/V2;
dxdt_CMT1 = -C1*Q + C2*Q - C1*CL ;
dxdt_CMT2 =  C1*Q - C2*Q ;

                      //---TUMOR GROWTH---
double EFF = EMAX*C1*exp(-LAMDA*SOLVERTIME)/(EC50+C1);
dxdt_CMT3  =  TG*CMT3*log(TUMLIM/CMT3)-EFF*CMT3;
double TUMSLD = CMT3;

                     //---SURVIVAL WEIBULL---
double HAZRATEBASE =  LAMBS*ALPHS*pow((SOLVERTIME+DEL),(ALPHS-1));
double HAZRATECOV  =  TUMSLD*TSHAZ + ECOG*ECOGHAZ;
dxdt_CMT4  = HAZRATEBASE*exp(HAZRATECOV+IIVHAZ);

                     //---DROPOUT WEIBULL---
dxdt_CMT5  = LAMBC*ALPHC*pow((SOLVERTIME+DEL),(ALPHC-1));
dxdt_CMT6  = CT;

$SIGMA   // Residual Unexplained Variability	
         block     =  FALSE
labels    =  s(RESERR)
1	      // Proportional error

$TABLE	 
                    //---Time course of drug concentration---
double CONC    =   CMT1/V1;

                    //---Time course of tumor size---


                    //---TIme to death---

double CHAZS   =   CMT4;
double CHAZC   =   CMT5;

double SURS    =   exp(-CHAZS);
double SURC    =   exp(-CHAZC);

double HASDRP    =   0;
double HASEVT    =   0;
double CENSOR    =   0 ; 

if ((HASDRP == 0) &  (SURC < UCENSOR)){
HASDRP = 1 ;
CENSOR  =  0 ; } 
if ((HASDRP == 0) &  (HASEVT ==0) & (SURS < UEVENT)){
HASEVT =  1 ;
CENSOR = 1;}
if ((HASDRP == 0) & (HASEVT == 0) & (TIME >= 730)){
HASDRP = 1;
CENSOR = 0;}
 

double AUC   =   CMT6;

//Determines Values and Includes in Output	
         double IPRED =  C1;                     //real concentration
         double DV   =  (IPRED)*(1+RESERR); // observed concentration


$CAPTURE 
IPRED,VC,V1,C1,DV,ETA1, ETA2,ETA3,ETA4,ETA5,ETA6,ETA7,ETA8     

'

# ------------------------------------------------------------------------------
# Compile the model code
mod <- mcode("NivolumabmodelSurvival",code)






