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
  TUM = 60,  // Tumor
  SRV = 0,  // Survival
  DRP = 0,  // Dropout
  AUC = 0,  // Area under the curve

$SET  // Set Differential Equation Solver Options			
  atol = 1e-8, rtol = 1e-8
  maxsteps = 100000

$PARAM  // Population parameters
  AVCLH = 0.06,  // clearance at health status (L/h)
  VC = 5.0,         // Central Volume (L)
  Q = 0.5,       // Intercompartmental CL (L/h)
  VP = 5.0,         // Peripheral Volume (L)
  TVEMAX = 0.02,    // 
  TVEC50 = 20,      //
  TVTG = 0.005,     // Tumor gowth rate
  TVLAMDA = 0.01,    // Resistance
  TUMLIM = 1000,

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

  // Default Covariate Values for Simulation
  ECOG = 1,
  UCENSOR = 1,
  UEVENT = 1,

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

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  labels = s(ZCL, ZVC, ZVP, ZEMAX, ZEC50, ZTG, ZR, ZHZ)
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
  label = s(RESERR)
  1  // Error defined as THETA in $PARAM

$MAIN // Drug Exposure
  // Covariate Values
  double DELT = 1;
  double TUMS = TUM + DELT;
  double CLH = AVCLH*pow(ECOG_CL, ECOG)*pow(TUMS, TUMSLD_CL);

  // Individual Parameter Values
  double CL = CLH*exp(ETA1); 
  double V1 = VC*exp(ETA2);
  double V2 = VP*exp(ETA3);
  
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
  double EFF = EMAX*C1*exp(-LAMDA*SOLVERTIME)/(EC50 + C1); // Divided by 24?
  dxdt_TUM = TG*TUMSLD*log(TUMLIM/TUMSLD) - EFF*TUMSLD;

  // Survival Weibull
  double HAZRATEBASE = LAMBS*ALPHS*pow((SOLVERTIME + DEL), (ALPHS - 1)); // Divided by 24?
  double HAZRATECOV = TUMSLD*TSHAZ + ECOG*ECOGHAZ;
  dxdt_SRV  = HAZRATEBASE*exp(HAZRATECOV + IIVHAZ);

  // Dropout Weibull
  dxdt_DRP  = LAMBC*ALPHC*pow((SOLVERTIME + DEL), (ALPHC - 1)); // Divided by 24?

$TABLE  // Determine values and output
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
  if ((HASDRP == 0) & (HASEVT == 0) & (TIME >= 730*24)) {
  HASDRP = 1;
  CENSOR = 0;}

$CAPTURE 
  ECOG IPRED DV AUC CL V1 V2 Q EMAX EC50 TG LAMDA IIVHAZ 
  HASDRP HASEVT UEVENT CENSOR UCENSOR
  C1 C2 CLH EFF TUMSLD TUMS  // Debug Tumour Growth
  // HAZRATEBASE HASRATECOV CHAZS CHAZC SURS SURC  // Debug Time to Death
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8
'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mcode("NivoPKPD", code)
