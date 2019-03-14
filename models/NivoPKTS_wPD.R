# Nivolumab Population Pharmacokinetic Model - No Tumour Size
# ------------------------------------------------------------------------------
# Original PopPK Model from manuscript with the addition of a tumour size model 
#   emailed from the model creator.

# Model sourced from:
#   C Liu, J Yu, H Li et al. (2017) Association of Time‐Varying Clearance of 
#   Nivolumab With Disease Dynamics and Its Implications on Exposure Response 
#   Analysis. Clin. Pharmacol. Ther., 101: 657-666. doi:10.1002/cpt.656

# Reference individual according to supplementary material:

# A typical subject (reference) is female, weighing 80 kg, 
#   eGFR of 80 mL/min/1.73 m2 , serum albumin of 4 mg/dL, tumor size of 60 mm, 
#   cell type/histology of other (i.e., not squamous or non-squamous), 
#   performance status of 0, tumor type of NSCLC, and ADA assay negative. 
#   The reference values for continuous covariates were selected to approximate 
#   the median values.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load libraries
  # library(dplyr)
  # library(mrgsolve)

# Define model code

  code <- '
$INIT  // Initial Conditions for Compartments
  CMT1 =  0,   // Central Compartment
  CMT2 =  0,   // Peripheral Compartment
  TUM =  54.6,  // Tumor size
  SRV =  0,   // Survival
  DRP =  0,   // Dropout
  AUC  =  0,   // Area under the curve

$SET     // Set Differential Equation Solver Options			
  atol      =  1e-8, rtol = 1e-8
  maxsteps  =  100000

$PARAM  // Population Parameters
  // Pharmacokinetic Population parameters
  TVCL = 0.0095*24,  // Typical value of Clearance (L/day)
  TVVC = 3.87,       // Typical value of Central Volume (L)
  TVVP = 3.01,       // Typical value of Peripheral Volume (L)
  Q = 0.0331*24,     // Intercompartmental Clearance (L/day)

  AVBWT = 80,      // Population bodyeight
  AVGFR = 80,      // Population GFR
  AVALB = 4,       // Population albumin
  AVTS = 54.6,       // Population Tumour Size
  TVTmax = 0.218,  // Typical value of the maximal change of clearance relative to baseline
  T50 = 66,        // Time for 50% of maximal clearance change
  HILL = 7.82,     // hill coefficient

  // Tumour Growth population parameters 
  TVEMAX = 0.02,     // typical value of Emax (day-1)
  TVEC50 = 20,       // typical value of EC50 (day-1)
  TVTG = 0.005,     // tumor growth rate (day-1)
  TVLAMDA = 0.001,   // resistance (day-1)
  TUMLIM = 1000,     // tumour limit (mm)

  // Weibull probability density function parameters
  LAMBS = 0.0001,  // Scale parameter for survival model
  ALPHS = 1,       // Shape parameter for survival model
  LAMBC = 0.001,   // Scale parameter for drop out model
  ALPHC = 1,       // Shape parameter for drop out model

  // Covariate Effects
  CL_BWT = 0.738,         // Effect of Body weight on Clearance
  CL_GFR = 0.189,         // Effect of Renal function on Clearance
  CL_ALB = -0.723,        // Effect of Albumin on Clearance
  CL_ECOG = 0.092,          // Effect of Performance Status on Clearance
  CL_ADApos = 1.11,       // Effect of ADA positive on Clearance
  CL_ADAunk = 1.04,       // Effect of unknown ADA status on Clearance
  CL_TUMORRCC = 0.071,    // Effect of tumor type RCC on Clearance
  CL_TUMOROTH = -0.0411,  // Effect of tumor type not RCC on Clearance
  CL_TS = 0.111,          // Effect of tumor size on Clearance

  VC_BWT = 0.582,         // Effect of Weight on Central Volume
  VC_MALE = 0.11,         // Effect of Sex on Central Volume
  VC_CELL = -0.123,       // Effect of cell type on VC (SQ or NSQ)

  TSHAZ = 0.02,           // Tumor size on hazard
  ECOGHAZ = 1,              // Performance status on hazard

  // Additive and proportional errors
  AERR = 0,
  PERR = 0.194,

  // Default Covariate Values for Simulation
  GFR = 80,     // glomerular filtration (mL/min)
  BWT = 80,     // body weight (kg)
  ALB = 4,      // serum albumin (mg/dL)
  AGE = 60,     // Age (years)
  SEX = 0,      // Sex (Male = 1, Female = 0)
  ECOG = 0,     // Performance Status
  RCC = 0,      // Renal Cell Carcinoma (positive = 1, negative = 0)
  OTHERC = 1,   // Other Cancer (positive = 1, negative = 0)
  ADApos = 0,   // Anti-Drug Antibodies (pos = 1, neg = 0)
  ADAunk = 1,   // unknown Anti-Drug Antibody Status (pos = 1, neg = 0)
  SQNSQ = 0,    // cell histology Squamous or non-squamous (pos = 1, neg = 0)
  UCENSOR = 1,  // chance of dropout
  UEVENT = 1,   // chance of event

  // Default ETA Values for Simulation
  // Allocated in population so set to zero
  ETA1 = 0,  // ZCL
  ETA2 = 0,  // ZVC
  ETA3 = 0,  // ZVP
  ETA4 = 0,  // ZEMAX
  ETA5 = 0,  // ZEC50
  ETA6 = 0,  // ZTG
  ETA7 = 0,  // ZR
  ETA8 = 0,  // ZHZ
  ETA9 = 0,  // ZTMAX

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  labels = s(ZCL, ZVC, ZVP, ZEMAX, ZEC50, ZTG, ZR, ZHZ, ZTMAX)
  0.096721  // ZCL
  0.099225  // ZVC
  0.185761  // ZVP
  1.000000  // ZEMAX
  0.010000  // ZEC50
  0.100000  // ZTG
  0.500000  // ZR
  0.100000  // ZHZ
  0.044521  // ZTMAX

$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  label = s(RESERR)
  1  // Error defined as THETA in $PARAM

$MAIN  // Drug Exposure
  // Covariate Values
  double COVca = pow(exp(CL_ECOG), ECOG)*pow(exp(CL_TUMORRCC), RCC)*
    pow(exp(CL_TUMOROTH), OTHERC)*pow(exp(CL_ADApos), ADApos)*
    pow(exp(CL_ADAunk), ADAunk);
  double COVco = pow(BWT/AVBWT, CL_BWT)*pow(GFR/AVGFR, CL_GFR)*
    pow(ALB/AVALB, CL_ALB);
  double CLTSPK = TVCL*COVco*COVca;
  double VCTSPK = TVVC*pow(BWT/AVBWT, VC_BWT)*pow(exp(VC_MALE), SEX)*
    pow(exp(VC_CELL), SQNSQ);
  double CLtime = exp(Tmax*pow(TIME, HILL)/(pow(T50, HILL) + pow(TIME, HILL)));

  // Individual Parameter Values
  double CL = CLTSPK*exp(ETA1); 
  double V1 = VCTSPK*exp(ETA2);
  double V2 = TVVP*exp(ETA3);
  double Tmax = TVTmax + ETA9;

  // Tumour Growth
  // Individual Parameter Values
  double EMAX = TVEMAX*exp(ETA4);
  double EC50 = TVEC50*exp(ETA5);

  double TG = TVTG*exp(ETA6);
  double LAMDA = TVLAMDA*exp(ETA7);

  // Survival Weibull
  double IIVHAZ = (ETA8);

$ODE  // Differential Equations
  // Drug Exposure
  double DEL = pow(10, -6);
  double C1 = CMT1/V1;
  double C2 = CMT2/V2;
  double CLTSPKtumcov = CL*pow(TUMSLD/AVTS, CL_TS);
  double CLTDPKtumcov = CLTSPKtumcov*CLtime;

  dxdt_CMT1 = -C1*Q + C2*Q - C1*CLTDPKtumcov ;
  dxdt_CMT2 =  C1*Q - C2*Q;
  dxdt_AUC = C1;

  // Tumour Growth
  double TUMSLD = TUM;
  double EFF = EMAX*C1*exp(-LAMDA*(SOLVERTIME))/(EC50 + C1);
  dxdt_TUM = TG*TUMSLD*log(TUMLIM/TUMSLD) - EFF*TUMSLD;

  // Survival Weibull
  double HAZRATEBASE = LAMBS*ALPHS*pow((SOLVERTIME + DEL), (ALPHS - 1)); 
  double HAZRATECOV = TUMSLD*TSHAZ + ECOG*ECOGHAZ;
  double HAZRATE = HAZRATEBASE*exp(HAZRATECOV + IIVHAZ);
  dxdt_SRV = HAZRATE;

  // Dropout Weibull
  dxdt_DRP = LAMBC*ALPHC*pow((SOLVERTIME + DEL), (ALPHC - 1));

$TABLE  // Determines Values and Includes in Output	
  // Drug Exposure
  double IPRED = C1;               // real concentration
  double DV = IPRED*(1 + RESERR);  // observed concentration

  // Time to Death
  double CHAZS = SRV;
  double CHAZC = DRP;
  
  double SURS = exp(-CHAZS);
  double SURC = exp(-CHAZC);
  
  double HASDRP = 0;
  double HASEVT = 0;
  double CENSOR = 0; 

  // If patient hasnt dropped out
  // Patient drops out if non-dropout chance too low
  if ((HASDRP == 0) &  (SURC < UCENSOR)) {
  HASDRP = 1;
  CENSOR = 0;} 
  // Patient dies if survival chance too low (if not already dead)
  if ((HASDRP == 0) &  (HASEVT == 0) & (SURS < UEVENT)) {
  HASEVT = 1;
  CENSOR = 1;}
  // Patient drops out if they are in study for 2 years (if not dead)
  if ((HASDRP == 0) & (HASEVT == 0) & (TIME >= 730)) {
  HASDRP = 1;
  CENSOR = 0;}

$CAPTURE 
  SEX AGE BWT ECOG IPRED DV AUC CL V1 V2 Q Tmax EMAX EC50 C1 C2 
  TG LAMDA IIVHAZ HASDRP HASEVT UEVENT CENSOR UCENSOR
  C1 C2 // COVca COVco CLTSPK VCTSPK CLTSPKtumcov CLtime CLTDPKtumcov  // Debug
  EFF TUMSLD HAZRATE // HAZRATEBASE HASRATECOV CHAZS CHAZC SURS SURC  // Debug
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 ETA9
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mcode("NivoPKTS", code)