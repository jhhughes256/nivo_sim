# Nivolumab Population Pharmacokinetic Model (20170803_Nivolumab18)

# ------------------------------------------------------------------------------
# Define the model parameters and equations

library(dplyr)	    #New plyr - required for mrgsolve
library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics


code <- '

$INIT    // Initial Conditions for Compartments
          CMT1      =  0,         // Central Compartment
          CMT2      =  0,         // Peripheral Compartment
          AUC       =  0,         // Area under the curve


$SET     // Set Differential Equation Solver Options			
atol      =  1e-8, rtol = 1e-8
maxsteps  =  100000


$PARAM   // Population parameters
TVCL     =  0.0095 ,        // Typical value of Clearance
TVVC     =  3.87,  	        // Typical value of Central Volume
TVVP     =  3.01,  	        // Typical value of Peripheral Volume
Q        =  0.0331	 ,        // Intercompartmental Clearance
AVBWT   = 80,               // Population bodyeight
AVGFR   = 80,               // Population GFR
AVALB   =  4,               // Population albumin
AVTS    =  60,

// Covariate Effects
CL_BWT          =  0.738,  // Effect of Body weight on Clearance
CL_GFR          =  0.189,  // Effect of Renal function on Clearance
CL_ALB          =  -0.723, // Effect of Albumin on Clearance
CL_PS           =  0.092,  // Effect of Performance Status on Clearance
CL_ADApos       =  1.11,   // Effect of ADA positive on Clearance
CL_ADAunk       = 1.04,    // Effect of ADA unknown on Clearance
CL_TUMORRCC     = 0.071,   // Effect of tumor type = RCC on Clearance
CL_TUMOROTH     =-0.0411,   // Effect of tumor type = OTHER on Clearance
CL_TS           = 0.111,    // Effect of tumor size on Clearance

VC_BWT   =  0.582,       // Effect of Weight on Central Volume
VC_MALE =   0.11,        // Effect of Sex on Central Volume
VC_CELL = -0.123,        // Effect of cell type on VC (SQ or NSQ)

CLH       = 0.06,       //   Typical value of disease-severity-irrelevant clearance ??
TVTmax    = 0.218  ,    // Typical value of the maximal change of clearance relative to baseline
T50       = 66    ,     // Time for 50% of maximal clearance change
HILL      = 7.82,       // hill coefficient
   

// Additive and proportional errors
AERR = 0 ,
PERR = 0.194,


// Default Covariate Values for Simulation (allocated in creation of population)
GFR       =  80,        // 80mL/min
BWT       =  80,        // 80kg
ALB       =   4,        // serum albumin = 4mg/dL
AGE       =  60, 
SEX      =   1,
PS        =   1,
RCC       =   1,
OTHERC    =   1,
ADApos    =   1,
ADAunk    =   1,
SQNSQ     =   1,
TS        =   60,

// Default ETA Values for Simulation (allocated in creation of population)
         ETA1      =  0,         // ZCL
         ETA2      =  0,         // ZVC
         ETA3      =  0,        // ZVP
         ETA4      =  0,        // ZTMAX


$OMEGA   // Population parameter Variability
name      = "omega1"
block     =  FALSE
labels    =  s(ZCL, ZVC, ZVP, ZTMAX)
0.096721           // ZCL
0.099225           // ZVC
0.185761           // ZVP
0.044521           // ZTMAX




$SIGMA   // Residual Unexplained Variability	
block     =  FALSE
labels    =  s(RESERR)
1	      // Proportional error


$MAIN    // Individual Parameter Values

double COVca = pow ( exp ( CL_PS ),PS ) * pow ( exp (CL_TUMORRCC) , RCC  ) *pow ( exp (CL_TUMOROTH) , OTHERC  ) * pow ( exp ( CL_ADApos ) , ADApos ) * pow ( exp ( CL_ADAunk), ADAunk);
double COVco = pow ( BWT / AVBWT , CL_BWT)  * pow ( GFR / AVGFR , CL_GFR)  * pow ( ALB / AVALB , CL_ALB);
double CLTSPK = TVCL*COVco*COVca*exp(ETA1);
double VCTSPK = TVVC * pow( BWT, VC_BWT)*pow(exp(VC_MALE),SEX) * pow(exp(VC_CELL),SQNSQ) * exp(ETA2);
double CLTSPKtumcov = CLTSPK * pow (TS/AVTS , CL_TS);
double Tmax = TVTmax + (ETA4) ;
double CLtime = exp (Tmax*pow(TIME,HILL) / (pow(T50,HILL)+pow(TIME,HILL)));
double CLTDPKtumcov = CLTSPKtumcov * CLtime;



double CL=CLTDPKtumcov ; 
double V1=VCTSPK;
double V2=TVVP*exp(ETA3) ;
double S1=V1;


D_CMT1    =  1;       // Infusion duration


$ODE     // Differential Equations
double C1 =  CMT1/V1;
double C2 =  CMT2/V2;

dxdt_CMT1 = -C1*Q + C2*Q - C1*CL ;
dxdt_CMT2 =  C1*Q - C2*Q;
dxdt_AUC  =  60*C1/1000;


$TABLE	 //Determines Values and Includes in Output	
double IPRED =  C1;                   //real concentration
double DV   =  IPRED*(1+RESERR); // observed concentration

$CAPTURE 
SEX AGE BWT  IPRED DV CL V1 V2 Q C1 C2 COVca COVco CLTSPK VCTSPK CLTSPKtumcov Tmax CLtime CLTDPKtumcov ETA1 ETA2 ETA3 ETA4

'

# ------------------------------------------------------------------------------
# Compile the model code
mod <- mcode("Nivolumabmodel1",code)