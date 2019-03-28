# Nivolumab Population Pharmacokinetic/Pharmacodynamic Model
# ------------------------------------------------------------------------------
# Second PK/PD model based on the email sent from model creator.

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
$INIT    // Initial Conditions for Compartments
  CMT1 = 0,  // Central Compartment
  CMT2 = 0,  // Peripheral Compartment
  TUM = 54.6,  // Tumor
  SRV = 0,  // Survival
  DRP = 0,  // Dropout
  AUC = 0,  // Area under the curve

$SET  // Set Differential Equation Solver Options			
  atol = 1e-8, rtol = 1e-8
  maxsteps = 100000

$PARAM  // Population parameters
  AVCLH = 0.06,     // clearance at health status (L/day)
  VC = 5.0,            // Central Volume (L)
  Q = 0.5,          // Intercompartmental CL (L/day)
  VP = 5.0,            // Peripheral Volume (L)
  TVEMAX = 0.02,    // Maximal effect on tumor suppression (day-1)
  TVEC50 = 20,         // EC50 that gives half-maximal effect (ug/mL)
  TVTG = 0.005,     // Tumor growth rate (day-1)
  TVLAMDA = 0.001,  // Resistance parameter of tumor (day-1)
  TUMLIM = 1000,       // Greatest possible size of tumour (mm)

  // Weibull probability density function parameters
  LAMBS = 0.0001,  // Scale parameter for survival model
  ALPHS = 1,       // Shape parameter for survival model
  LAMBC = 0.001,   // Scale parameter for drop out model
  ALPHC = 1,       // Shape parameter for drop out model

  // Covariate Effects
  TUMSLD_CL = 0.3,  // Power of TUMSLD effect on CL
  ECOG_CL = 1.5,    // ECOG effect on CL
  TSHAZ = 0.02,     // Tumor size on hazard
  ECOGHAZ = 1,      // ECOG on hazard

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

  // Default EPS values for simulation
  // Allocated in population so set to zero
  EPS1 = 0,  // EPROP

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  0.10  // ZCL
  0.10  // ZVC
  0.10  // ZVP
  1.00  // ZEMAX
  0.01  // ZEC50
  0.10  // ZTG
  0.50  // ZR
  0.10  // ZHZ

$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  1  // Error defined as THETA in $PARAM

$MAIN // Drug Exposure
  D_CMT1 = 0.5/24;    // Infusion Duration (days)

  // Covariate Values
  double DELT = 1;
  double TUMS = TUM + DELT;
  double CLH = AVCLH*pow(ECOG_CL, ECOG)*pow(TUMS, TUMSLD_CL);

  // Individual Parameter Values
  double CL = CLH*exp(ETA1); 
  double CLi = CL;
  double V1 = VC*exp(ETA2);
  double V2 = VP*exp(ETA3);
  double TMAX = 0;
  
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
  double C1 =  CMT1/V1;
  double C2 =  CMT2/V2;

  dxdt_CMT1 = -C1*Q + C2*Q - C1*CL;
  dxdt_CMT2 =  C1*Q - C2*Q;
  dxdt_AUC  =  C1;

  // Tumour Growth
  double TUMSLD = TUM;
  double EFF = EMAX*C1*exp(-LAMDA*SOLVERTIME)/(EC50 + C1); 
  dxdt_TUM = TG*TUMSLD*log(TUMLIM/TUMSLD) - EFF*TUMSLD;

  // Survival Weibull
  double HAZRATEBASE = LAMBS*ALPHS*pow((SOLVERTIME + DEL), (ALPHS - 1)); 
  double HAZRATECOV = TUMSLD*TSHAZ + ECOG*ECOGHAZ;
  double HAZRATE = HAZRATEBASE*exp(HAZRATECOV + IIVHAZ);
  dxdt_SRV = HAZRATE;

  // Dropout Weibull
  dxdt_DRP = LAMBC*ALPHC*pow((SOLVERTIME + DEL), (ALPHC - 1)); 

$TABLE  // Determine values and output
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
  CL CLi V1 V2 Q EMAX EC50 TG LAMDA IIVHAZ TMAX  // Individual Parameters
  // C1 C2 SRV DRP HAZRATEBASE HASRATECOV CHAZS CHAZC SURS SURC // Debug
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 ETA9 UEVENT UCENSOR  // Variability
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mcode("NivoPKPD", code)
