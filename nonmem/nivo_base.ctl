;; 1. Based on: 005
;; 2. Description: sim_FINAL
;; x1. Author: CL
$SIZES MAXFCN=100000000
$PROBLEM    sim_PK_TUM
;Overall Survival Simulations
$INPUT      ID TIME DV BASESLD ECOG AMT II ADDL MDV EVID FLAG DOSEGRP


$DATA       nonmem_pop.csv IGNORE=# 
$SUBROUTINE ADVAN9 TOL=6

$MODEL      NCOMP=6

COMP  =  (CENTRAL)
COMP  =  (PERIPH)
COMP  =  (TUMOR)
COMP  =  (SURVIVAL)
COMP  =  (CENSOR)
COMP  =  (AUC)

$THETA
(0.06) FIX      ;--th1-CLH: Clearance (L/h) at health status
(0.3) FIX      ;--th2-CL: Power of TUMSLD effect on CL
(1.5) FIX      ;--th3-CL: ECOG effect on CL
(5.0) FIX      ;--th4-VC: Central Volume (L)
(0.5) FIX      ;--th5-Q: Intercompartmental CL (L/h)
(5.0) FIX     ;--th6-VP: Peripheral Volume (L)

(0.02) FIX    ;--th7-TVEMAX
(20) FIX       ;--th8-TVEC50
(0.005) FIX   ;--th9-TG: Tumor gowth rate
(0.001) FIX   ;--th10-LAMDA: Resistance

(0.0001) FIX   ;--th11-LAMBS: Scale parameter in the Weibull probability density function for the survival model
(1) FIX        ;--th12-ALPHS: Shape parameter in the Weibull probability density function for the survival model
(0.001) FIX    ;--th13-LAMBC: Scale parameter in the Weibull probability density function for the drop out model
(1) FIX        ;--th14-ALPHC: Shape parameter in the Weibull probability density function for the drop out model

(0.02) FIX     ;--th15-TSHAZ: Tumor size on hazard
(1) FIX      ;--th16-ECOGHAZ: ECOG on hazard

$OMEGA
0.1 FIX  ;--eta1-IIV in CL [exponential]
0.1 FIX  ;--eta2-IIV in VC [exponential]
0.1 FIX  ;--eta3-IIV in VP [exponential]
1 FIX  ;--eta4-IIV in EMAX [exponential]
0.01 FIX  ;--eta5-IIV in EC50 [exponential]
0.1 FIX  ;--eta6-IIV in Tumor growth rate [exponential]
0.5 FIX  ;--eta7-IIV in drug resistance [exponential]
0.1 FIX  ;--eta8-IIV in hazard rate [additive]

$PK
;-----DRUG EXPOSURE------------------
DELT = 1
CLH = THETA(1)*EXP(ETA(1))
TUM = A(3) + DELT
CL  =  CLH * THETA(3)**ECOG * TUM**THETA(2)

TVV1 = THETA(4)
V1  =  TVV1*EXP(ETA(2))
TVQ=THETA(5)
Q=TVQ
TVV2=THETA(6)
V2=TVV2*EXP(ETA(3))

K10=CL/V1
K12=Q/V1
K21=Q/V2
S1=V1 ;dose=mg,volume in L, conc=ug/mL

A_0(1) = 0
A_0(2) = 0
A_0(3) = BASESLD
A_0(4) = 0
A_0(5) = 0


;----- TUMOR GROWTH ------------
TVEMAX =   THETA(7)
TVEC50 =   THETA(8)

EMAX   =   TVEMAX*EXP(ETA(4))
EC50   =   TVEC50*EXP(ETA(5))

TG = THETA(9)*EXP(ETA(6))
LAMDA = THETA(10)*EXP(ETA(7))

;----- Weibull model for survival -------

LAMBS    =   THETA(11) ; Scale parameter in the Weibull probability density function for the survival model
ALPHS    =   THETA(12) ; Shape parameter in the Weibull probability density function for the survival model
LAMBC    =   THETA(13) ; Scale parameter in the Weibull probability density function for the drop out model
ALPHC    =   THETA(14) ; Shape parameter in the Weibull probability density function for the drop out model

TSHAZ    =   THETA(15) ; Tumor size on hazard
ECOGHAZ  =   THETA(16) ; ECOG on hazard

IIVHAZ   =   ETA(8)

;----- Generate Random Number -------

IF (ICALL.EQ.4) THEN
IF (NEWIND.LE.1) THEN  ; first record

CALL RANDOM(2,R)
UEVENT=R                 ; Uniform random number for event

CALL RANDOM(3,R)
UCENSOR=R                ; Uniform random number for dropout

ENDIF
ENDIF

$DES

;-----DRUG EXPOSURE------------------
DEL=1E-6
DADT(1) =  K21*A(2) - A(1)*K12 - A(1)*K10
DADT(2) =  K12*A(1) - K21*A(2)
CT = A(1)/V1

;-----TUMOR GROWTH ------------------

EFF = EMAX*CT*EXP(-LAMDA*T)/(EC50+CT)
TUMLIM = 1000

DADT(3)   =   TG*A(3)*LOG(TUMLIM/A(3)) - EFF*A(3)
TUMSLD    =   A(3)    ;tumor size time course

;-----Weibull model for survival-------


HAZRATEBASE =  LAMBS*ALPHS*(T+DEL)**(ALPHS-1)
HAZRATECOV  =  TUMSLD*TSHAZ + ECOG*ECOGHAZ
DADT(4)     =  HAZRATEBASE*EXP(HAZRATECOV+IIVHAZ)

;-----Weibull model for dropout

DADT(5) = LAMBC*ALPHC*(T+DEL)**(ALPHC-1)

DADT(6) = CT

$SIGMA 0 FIX

$ERROR
;-----------Time course of drug concentration------
CONC  =  A(1)/V1
;-----------Time course of tumor size--------------
TUM_SLD    =   A(3)    ;tumor size time course

;-----------TTD (Time to death)--------------
CHAZS     =   A(4)            ; cumulative hazard for survival
CHAZC     =   A(5)            ; cumulative hazard for drop out

SURS      =   EXP(-CHAZS)      ; survival probability
SURC      =   EXP(-CHAZC)      ; drop out probability

HASDRP    =   0
HASEVT    =   0
CENSOR    =   0

IF (ICALL.EQ.4) THEN
IF (FLAG.EQ.1.AND.HASDRP.EQ.0.AND.SURC.LT.UCENSOR) THEN
HASDRP     =   1 ; dropout event
CENSOR     =   0 ; censored event
ELSE
IF (FLAG.EQ.1.AND.HASDRP.EQ.0.AND.HASEVT.EQ.0.AND.SURS.LT.UEVENT) THEN
HASEVT     =   1 ; Event
CENSOR     =   1
ELSE
IF(FLAG.EQ.1.AND.HASDRP.EQ.0.AND.HASEVT.EQ.0.AND.TIME.GE.730) THEN
HASDRP     =   1
CENSOR     =   0
ENDIF
ENDIF
ENDIF
ENDIF

AUC = A(6)

$SIM (12345) (13579 UNIFORM) (24680 UNIFORM) ONLYSIM SUBPROB=1000

$TABLE ID TIME BASESLD DOSEGRP CL CLH ECOG EVID CONC TUM_SLD EMAX EC50 TG CL LAMDA EFF CENSOR
HASEVT HASDRP SURS SURC UCENSOR UEVENT AUC 

ONEHEADER NOPRINT FILE=nonmem_pop.fit
