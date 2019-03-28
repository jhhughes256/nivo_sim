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
  TVCL = 0.0082*24,  // Typical value of Clearance (L/day)
  TVVC = 3.86,       // Typical value of Central Volume (L)
  TVVP = 3.73,       // Typical value of Peripheral Volume (L)
  Q = 0.0307*24,     // Intercompartmental Clearance (L/day)

  AVBWT = 80,      // Population body weight (kg) 
  AVGFR = 80,      // Population GFR (mL/min/1.73m2)
  AVALB = 4,       // Population albumin (paper: mg/dl; reality?: g/dl)
  AVTS = 54.6,     // Population Tumour Size (mm)

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
  CL_BWT = 0.724,         // Effect of Body weight on Clearance
  CL_GFR = 0.187,         // Effect of Renal function on Clearance
  CL_ALB = -0.72,        // Effect of Albumin on Clearance
  CL_ECOG = 0.093,          // Effect of Performance Status on Clearance
  CL_ADApos = 1.11,       // Effect of ADA positive on Clearance
  CL_ADAunk = 0.988,       // Effect of unknown ADA status on Clearance
  CL_TUMORRCC = 0.0551,    // Effect of tumor type RCC on Clearance
  CL_TUMOROTH = -0.0463,  // Effect of tumor type not RCC on Clearance
  CL_TS = 0.126,          // Effect of tumor size on Clearance

  VC_BWT = 0.582,         // Effect of Weight on Central Volume
  VC_MALE = 0.11,         // Effect of Sex on Central Volume
  VC_CELL = -0.123,       // Effect of cell type on VC (SQ or NSQ)

  TSHAZ = 0.02,           // Tumor size on hazard
  ECOGHAZ = 1,              // Performance status on hazard

  // Additive and proportional errors
  AERR = 0,
  PERR = 0.199,

  // Default Covariate Values for Simulation
  GFR = 80,     // glomerular filtration (mL/min)
  BWT = 80,     // body weight (kg)
  ALB = 4,      // serum albumin (mg/dL)
  AGE = 60,     // Age (years)
  SEX = 0,      // Sex (Male = 1, Female = 0)
  ECOG = 0,     // Performance Status
  RCC = 1,      // Renal Cell Carcinoma (positive = 1, negative = 0)
  OTHERC = 0,   // Other Cancer (positive = 1, negative = 0)
  ADApos = 0,   // Anti-Drug Antibodies (pos = 1, neg = 0)
  ADAunk = 0,   // unknown Anti-Drug Antibody Status (pos = 1, neg = 0)
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

  // Default EPS values for simulation
  // Allocated in population so set to zero
  EPS1 = 0,  // EPROP

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  0.119025  // ZCL
  0.099225  // ZVC
  0.199809  // ZVP
  1.000000  // ZEMAX
  0.010000  // ZEC50
  0.100000  // ZTG
  0.500000  // ZR
  0.100000  // ZHZ

$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  1  // Error defined as THETA in $PARAM

$MAIN  // Drug Exposure
  D_CMT1 = 0.5/24;    // Infusion Duration (days)

  // Covariate Values
  double COVca = pow(exp(CL_ECOG), ECOG)*pow(exp(CL_TUMORRCC), RCC)*
    pow(exp(CL_TUMOROTH), OTHERC)*pow(exp(CL_ADApos), ADApos)*
    pow(exp(CL_ADAunk), ADAunk);
  double COVco = pow(BWT/AVBWT, CL_BWT)*pow(GFR/AVGFR, CL_GFR)*
    pow(ALB/AVALB, CL_ALB);
  double CLTSPK = TVCL*COVco*COVca;
  double VCTSPK = TVVC*pow(BWT/AVBWT, VC_BWT)*pow(exp(VC_MALE), SEX)*
    pow(exp(VC_CELL), SQNSQ);

  // Individual Parameter Values
  double CLi = CLTSPK*exp(ETA1); 
  double V1 = VCTSPK*exp(ETA2);
  double V2 = TVVP*exp(ETA3);

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
  double CL = CLi*pow(TUMSLD/AVTS, CL_TS);

  dxdt_CMT1 = -C1*Q + C2*Q - C1*CL ;
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
  double DV = IPRED*(1 + PERR*EPS1);  // observed concentration

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
  SEX AGE BWT GFR ALB ECOG RCC OTHERC ADApos ADAunk SQNSQ  // Covariates
  IPRED DV AUC EFF TUM HAZRATE HASDRP HASEVT CENSOR  // Outputs
  CL CLi V1 V2 Q EMAX EC50 TG LAMDA IIVHAZ  // Individual Parameters
  // C1 C2 SRV DRP HAZRATEBASE HASRATECOV CHAZS CHAZC SURS SURC // Debug
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 EPS1 UEVENT UCENSOR  // Variability
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mcode("NivoPKTS", code)