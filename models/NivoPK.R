# Nivolumab Population Pharmacokinetic Model - No Tumour Size
# ------------------------------------------------------------------------------
# Original PopPK Model based on the model written in the manuscript mentioned
#   below.

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
  CMT1 =  0,  // Central Compartment
  CMT2 =  0,  // Peripheral Compartment
  AUC =  0,   // Area under the curve

$SET  // Set Differential Equation Solver Options			
  atol = 1e-8, rtol = 1e-8
  maxsteps = 100000

$PARAM  // Population parameters
  TVCL = 0.0095,   // Typical value of Clearance
  TVVC = 3.87,     // Typical value of Central Volume
  TVVP = 3.01,     // Typical value of Peripheral Volume
  Q = 0.0331,      // Intercompartmental Clearance

  AVBWT = 80,      // Population bodyeight
  AVGFR = 80,      // Population GFR
  AVALB = 4,       // Population albumin
  AVTS = 60,       // Population Tumour Size
  CLH = 0.06,      // Typical value of disease-severity-irrelevant CL
  TVTmax = 0.218,  // Typical value of the maximal change of clearance relative to baseline
  T50 = 66,        // Time for 50% of maximal clearance change
  HILL = 7.82,     // hill coefficient

  // Covariate Effects
  CL_BWT = 0.738,         // Effect of Body weight on Clearance
  CL_GFR = 0.189,         // Effect of Renal function on Clearance
  CL_ALB = -0.723,        // Effect of Albumin on Clearance
  CL_PS = 0.092,          // Effect of Performance Status on Clearance
  CL_ADApos = 1.11,       // Effect of ADA positive on Clearance
  CL_ADAunk = 1.04,       // Effect of unknown ADA status on Clearance
  CL_TUMORRCC = 0.071,    // Effect of tumor type RCC on Clearance
  CL_TUMOROTH = -0.0411,  // Effect of tumor type not RCC on Clearance
  CL_TS = 0.111,          // Effect of tumor size on Clearance

  VC_BWT = 0.582,         // Effect of Weight on Central Volume
  VC_MALE = 0.11,         // Effect of Sex on Central Volume
  VC_CELL = -0.123,       // Effect of cell type on VC (SQ or NSQ)

  // Additive and proportional errors
  AERR = 0,
  PERR = 0.194,

  // Default Covariate Values for Simulation
  GFR = 80,    // glomerular filtration (mL/min)
  BWT = 80,    // body weight (kg)
  ALB = 4,     // serum albumin (mg/dL)
  AGE = 60,    // Age (years)
  SEX = 0,     // Sex (Male = 1, Female = 0)
  PS = 0,      // Performance Status
  RCC = 0,     // Renal Cell Carcinoma (positive = 1, negative = 0)
  OTHERC = 1,  // Other Cancer (positive = 1, negative = 0)
  ADApos = 0,  // Anti-Drug Antibodies (positive = 1, negative = 0)
  ADAunk = 1,  // unknown Anti-Drug Antibody Status (positive = 1, negative = 0)
  SQNSQ = 0,   // cell type/histology Squamous or non-squamous (pos = 1, neg = 0)
  TS = 60,     // Tumour Size

  // Default ETA Values for Simulation
  // Allocated in population so set to zero
  ETA1 = 0,  // ZCL
  ETA2 = 0,  // ZVC
  ETA3 = 0,  // ZVP
  ETA4 = 0,  // ZTMAX

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  labels = s(ZCL, ZVC, ZVP, ZTMAX)
  0.096721  // ZCL
  0.099225  // ZVC
  0.185761  // ZVP
  0.044521  // ZTMAX

$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  label = s(RESERR)
  1  // Error defined as THETA in $PARAM

$MAIN  // Determine covariate values
  double COVca = pow(exp(CL_PS), PS)*pow(exp(CL_TUMORRCC), RCC)*
    pow(exp(CL_TUMOROTH), OTHERC)*pow(exp(CL_ADApos), ADApos)*
    pow(exp(CL_ADAunk), ADAunk);
  double COVco = pow(BWT/AVBWT, CL_BWT)*pow(GFR/AVGFR, CL_GFR)*
    pow(ALB/AVALB, CL_ALB);
  double CLTSPK = TVCL*COVco*COVca;
  double VCTSPK = TVVC*pow(BWT, VC_BWT)*pow(exp(VC_MALE), SEX)*
    pow(exp(VC_CELL), SQNSQ);
  double CLTSPKtumcov = CLTSPK*pow(TS/AVTS, CL_TS);
  double CLtime = exp(Tmax*pow(TIME, HILL)/(pow(T50, HILL) + pow(TIME, HILL)));
  double CLTDPKtumcov = CLTSPKtumcov*CLtime;

  // Individual Parameter Values
  double CL = CLTDPKtumcov*exp(ETA1); 
  double V1 = VCTSPK*exp(ETA2);
  double V2 = TVVP*exp(ETA3);
  double Tmax = TVTmax + ETA4;

$ODE  // Differential Equations
  double C1 = CMT1/V1;
  double C2 = CMT2/V2;

  dxdt_CMT1 = -C1*Q + C2*Q - C1*CL ;
  dxdt_CMT2 =  C1*Q - C2*Q;
  dxdt_AUC = 60*C1/1000;

$TABLE  // Determines Values and Includes in Output	
  double IPRED = C1;               // real concentration
  double DV = IPRED*(1 + RESERR);  // observed concentration

$CAPTURE 
  SEX AGE BWT IPRED DV CL V1 V2 Q Tmax C1 C2
  // COVca COVco CLTSPK VCTSPK CLTSPKtumcov CLtime CLTDPKtumcov  // Debug
  ETA1 ETA2 ETA3 ETA4
'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mcode("NivoPK", code)