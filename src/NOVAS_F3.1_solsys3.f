*** SOLSYS VERSION 3 PACKAGE: SOLSYS, SUN, IDSS ***

      SUBROUTINE SOLSYS (TJD,M,K,POS,VEL,IERR)
*
*     SUBROUTINE SOLSYS VERSION 3.
*     THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF THE
*     EARTH AT EPOCH TJD BY EVALUATING A CLOSED-FORM THEORY WITHOUT
*     REFERENCE TO AN EXTERNAL FILE.  THIS ROUTINE CAN ALSO PROVIDE
*     THE POSITION AND VELOCITY OF THE SUN.
*
*          TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
*          M    = BODY IDENTIFICATION NUMBER (IN)
*                 SET M=0 OR M=1 FOR THE SUN
*                 SET M=2 OR M=3 FOR THE EARTH
*          K    = ORIGIN SELECTION CODE (IN)
*                 SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
*                 SET K=1 FOR ORIGIN AT CENTER OF SUN
*          POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
*                 OF J2000.0, COMPONENTS IN AU (OUT)
*          VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
*                 OF J2000.0, COMPONENTS IN AU/DAY (OUT)
*          IERR = ERROR INDICATOR (OUT)
*                 IERR=0 MEANS EVERYTHING OK
*                 IERR=1 MEANS TJD BEFORE FIRST ALLOWED DATE
*                 IERR=2 MEANS TJD AFTER LAST ALLOWED DATE
*
*
       DOUBLE PRECISION TJD,POS,VEL,PI,TWOPI,T0,OBL,EL,C,P,TLAST,
     .     PM,PA,PE,PJ,PO,PW,PL,PN,       
     .     TMASS,SE,CE,SI,CI,SN,CN,SW,CW,P1,P2,P3,Q1,Q2,Q3,ROOTE,A,B,
     .     QJD,E,MLON,MA,U,SINU,COSU,ANR,PPLAN,VPLAN,F,PBARY,VBARY,
     .     DFLOAT,DABS,DMOD,DSIN,DCOS,DSQRT
      DIMENSION POS(3), VEL(3), EL(21), C(13), P(3,3),
     .     PM(4), PA(4), PE(4), PJ(4), PO(4), PW(4), PL(4), PN(4), 
     .     A(3,4), B(3,4), PPLAN(3), VPLAN(3), PBARY(3), VBARY(3)
      SAVE
      
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( TWOPI  = 2.D0 * PI             )
      PARAMETER ( T0     = 2451545.0D0           )
      PARAMETER ( OBL    = 23.43927944D0         )
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0
*     OBL = OBLIQUITY OF ECLIPTIC AT EPOCH J2000.0
      
      DATA EL, C, P / 43*0.D0 /,   TLAST / 0.D0 /

*     ARRAYS BELOW CONTAIN MASSES AND ORBITAL ELEMENTS OF THE FOUR
*     LARGEST PLANETS (SEE EXPLANATORY SUPPLEMENT (1992), P. 316)
*     WITH ANGLES IN RADIANS
*     THIS DATA USED FOR BARYCENTER COMPUTATIONS ONLY
*                 JUPITER        SATURN        URANUS       NEPTUNE
      DATA PM /  1047.349D 0,  3497.898D 0,   22903.0D 0,   19412.2D 0 /
      DATA PA /  5.203363D 0,  9.537070D 0, 19.191264D 0, 30.068963D 0 /
      DATA PE /  0.048393D 0,  0.054151D 0,  0.047168D 0,  0.008586D 0 /
      DATA PJ /  0.022782D 0,  0.043362D 0,  0.013437D 0,  0.030878D 0 /
      DATA PO /  1.755036D 0,  1.984702D 0,  1.295556D 0,  2.298977D 0 / 
      DATA PW /  0.257503D 0,  1.613242D 0,  2.983889D 0,  0.784898D 0 /
      DATA PL /  0.600470D 0,  0.871693D 0,  5.466933D 0,  5.321160D 0 /
      DATA PN /  1.450138D-3,  5.841727D-4,  2.047497D-4,  1.043891D-4 /
     
      IF ( TLAST .LT. 1.D0 ) THEN
*         FIRST TIME COMPUTATIONS
*         MASS OF SUN PLUS FOUR INNER PLANETS
          TMASS = 1.D0 + 5.977D-6
          SE = DSIN ( OBL * PI / 180.D0 )
          CE = DCOS ( OBL * PI / 180.D0 )
          DO 15 I = 1, 4
              TMASS = TMASS + 1.D0 / PM(I)
*             COMPUTE SINE AND COSINE OF ORBITAL ANGLES
			  SI = DSIN ( PJ(I) )
			  CI = DCOS ( PJ(I) )
			  SN = DSIN ( PO(I) )
			  CN = DCOS ( PO(I) )
			  SW = DSIN ( PW(I) - PO(I) )
			  CW = DCOS ( PW(I) - PO(I) )
*             COMPUTE P AND Q VECTORS (SEE BROUWER & CLEMENCE (1961), 
*             METHODS OF CELESTIAL MECHANICS, PP. 35-36.)
			  P1 =    CW * CN - SW * SN * CI
			  P2 = (  CW * SN + SW * CN * CI ) * CE - SW * SI * SE
			  P3 = (  CW * SN + SW * CN * CI ) * SE + SW * SI * CE
			  Q1 =   -SW * CN - CW * SN * CI
			  Q2 = ( -SW * SN + CW * CN * CI ) * CE - CW * SI * SE
			  Q3 = ( -SW * SN + CW * CN * CI ) * SE + CW * SI * CE
              ROOTE = DSQRT ( 1.D0 - PE(I)**2 )
              A(1,I) = PA(I) * P1
              A(2,I) = PA(I) * P2
              A(3,I) = PA(I) * P3
              B(1,I) = PA(I) * ROOTE * Q1
              B(2,I) = PA(I) * ROOTE * Q2
              B(3,I) = PA(I) * ROOTE * Q3
  15      CONTINUE
          TLAST = 1.D0
      END IF 
      
      IERR = 0
*     VALID DATES ARE WITHIN 3 CENTURIES OF J2000, ALTHOUGH RESULTS
*     DETERIORATE GRADUALLY      
      IF ( TJD .LT. 2340000.5D0 ) IERR = 1
      IF ( TJD .GT. 2560000.5D0 ) IERR = 2
      IF ( IERR .NE. 0 ) GO TO 110
      IF ( M .GE. 2 ) GO TO 30

*     FORM HELIOCENTRIC COORDINATES OF SUN
  20  DO 25 J=1,3
          POS(J) = 0.D0
          VEL(J) = 0.D0
  25  CONTINUE
      IF ( K .GE. 1 ) GO TO 110
      GO TO 90
      
*     FORM HELIOCENTRIC COORDINATES OF EARTH
*     VELOCITIES ARE OBTAINED FROM CRUDE NUMERICAL DIFFERENTIATION
  30  DO 35 I = 1, 3
          QJD = TJD + DFLOAT(I-2) * 0.1D0
C         SUBROUTINE SUN COMPUTES EARTH-SUN VECTOR
          CALL SUN ( QJD, EL, C )
          CALL PRECES ( QJD, C(11), T0, POS )
          P(I,1) = -POS(1)
          P(I,2) = -POS(2)
          P(I,3) = -POS(3)
  35  CONTINUE
      DO 40 J=1,3
          POS(J) =   P(2,J)
          VEL(J) = ( P(3,J) - P(1,J) ) / 0.2D0  
  40  CONTINUE
      IF ( K .GE. 1 ) GO TO 110

*     IF K=0, MOVE ORIGIN TO SOLAR SYSTEM BARYCENTER
*     SOLAR SYSTEM BARYCENTER COORDINATES COMPUTED FROM KEPLERIAN
*     APPROXIMATIONS OF THE COORDINATES OF THE FOUR LARGEST PLANETS
  90  IF ( DABS ( TJD - TLAST ) .LT. 1.D-6 ) GO TO 99
      DO 92 J = 1, 3
          PBARY(J) = 0.D0
          VBARY(J) = 0.D0
  92  CONTINUE 
*     THE FOLLOWING LOOP CYCLES ONCE FOR EACH OF THE FOUR LARGE PLANETS
      DO 98 I = 1, 4
*         COMPUTE MEAN LONGITUDE, MEAN ANOMALY, AND ECCENTRIC ANOMOLY
          E = PE(I)
          MLON = PL(I) + PN(I) * ( TJD - T0 )
          MA = DMOD ( MLON - PW(I), TWOPI )
          U = MA + E * DSIN ( MA ) + 0.5D0 * E * E * DSIN ( 2.D0 * MA )
          SINU = DSIN ( U )
          COSU = DCOS ( U )
*         COMPUTE VELOCITY FACTOR     
          ANR = PN(I) / ( 1.D0 - E * COSU )
*         COMPUTE PLANET'S POSITION AND VELOCITY WRT EQ & EQ J2000
          PPLAN(1) = A(1,I) * ( COSU - E ) + B(1,I) * SINU 
          PPLAN(2) = A(2,I) * ( COSU - E ) + B(2,I) * SINU
          PPLAN(3) = A(3,I) * ( COSU - E ) + B(3,I) * SINU
          VPLAN(1) = ANR * ( -A(1,I) * SINU + B(1,I) * COSU )
          VPLAN(2) = ANR * ( -A(2,I) * SINU + B(2,I) * COSU )
          VPLAN(3) = ANR * ( -A(3,I) * SINU + B(3,I) * COSU )
*         COMPUTE MASS FACTOR AND ADD IN TO TOTAL DISPLACEMENT      
          F = 1.D0 / ( PM(I) * TMASS )
          PBARY(1) = PBARY(1) + PPLAN(1) * F
          PBARY(2) = PBARY(2) + PPLAN(2) * F
          PBARY(3) = PBARY(3) + PPLAN(3) * F
          VBARY(1) = VBARY(1) + VPLAN(1) * F
          VBARY(2) = VBARY(2) + VPLAN(2) * F
          VBARY(3) = VBARY(3) + VPLAN(3) * F
  98  CONTINUE
      TLAST = TJD
  99  DO 100 J=1,3
          POS(J) = POS(J) - PBARY(J)
          VEL(J) = VEL(J) - VBARY(J)
 100  CONTINUE 

 110  RETURN

      END



      SUBROUTINE SUN (DJ,EL,C)
C
C     FOR USE WITH SUBROUTINE SOLSYS VERSION 3.
C     THIS SUBROUTINE COMPUTES THE COORDINATES OF THE EARTH-SUN
C     POSITION VECTOR WITH RESPECT TO THE ECLIPTIC AND EQUATOR
C     OF DATE.  A MODIFIED FORM OF NEWCOMB'S THEORY ('TABLES OF THE
C     SUN', 1898) IS USED.  ONLY THE LARGEST PERIODIC PERTURBATIONS
C     ARE EVALUATED, AND VAN FLANDERN'S EXPRESSIONS FOR THE FUNDAMENTAL
C     ARGUMENTS ('IMPROVED MEAN ELEMENTS FOR THE EARTH AND MOON', 1981)
C     ARE USED.  THE ABSOLUTE ACCURACY IS NO WORSE THAN 1 ARCSECOND
C     (AVERAGE ERROR ABOUT 0.2 ARCSECOND) OVER 1800-2200.
C     (ADAPTED FROM SUBROUTINE IAUSUN BY P. M. JANICZEK, USNO.)
C
C          DJ   = TDB JULIAN DATE OF DESIRED EPOCH (IN)
C          EL   = ARRAY OF ORBITAL ELEMENTS (SEE BELOW) FOR
C                 EPOCH DJ (OUT)
C          C    = ARRAY OF COORDINATES (SEE BELOW) FOR
C                 EPOCH DJ (OUT)
C
C
      DOUBLE PRECISION DJ,EL,C,T,TP,T20,RO,GV,GM,GJ,GS,DL,DR,DB,DG,
     1 DBLARG,D,TWOPI,STR,RTD,R,TR,
     2 SINO,COSO,SINL,COSL,SINB,COSB,
     3 DSIN,DCOS,DMOD
C
      DIMENSION EL(21)
C
C     EL( 1)= SEMI-MAJOR AXIS, AU
C     EL( 2)= ORBITAL ECCENTRICITY
C     EL( 5)= LONGITUDE OF PERIGEE, RADIANS
C     EL( 9)= UNPERTURBED MEAN LONGITUDE, RADIANS
C     EL(10)= MEAN ANOMALY, AFFECTED BY LONG-PD PERTURBATIONS, RADIANS
C     EL(11)= UNPERTURBED RADIUS, AU
C     EL(12)= EQUATION OF THE CENTER, RADIANS
C     EL(13)= MEAN OBLIQUITY OF ECLIPTIC, RADIANS
C     EL(14)= MEAN LONGITUDE OF MOON, RADIANS
C     EL(15)= MEAN ANOMALY OF MOON, RADIANS
C     EL(16)= LUNAR MEAN ARGUMENT OF LATITUDE, RADIANS
C     EL(17)= MEAN LONGITUDE OF LUNAR ASCENDING NODE, RADIANS
C     EL(21)= MEAN LONGITUDE OF MOON'S PERIGEE, RADIANS
C             (REMAINING ELEMENTS OF ARRAY EL NOT USED)
C
      DIMENSION C(13)
C
C     C( 1) = PERTURBED RADIUS VECTOR, AU
C     C( 2) = SAME AS C(4), DEGREES
C     C( 3) = SAME AS C(5), DEGREES
C     C( 4) = ECLIPTIC LONGITUDE WRT MEAN ECL & EQUX OF DATE, RADIANS
C     C( 5) = ECLIPTIC LATITUDE  WRT MEAN ECL        OF DATE, RADIANS
C     C(11) = EQUATORIAL X WRT MEAN EQU & EQUX OF DATE, AU
C     C(12) = EQUATORIAL Y WRT MEAN EQU & EQUX OF DATE, AU
C     C(13) = EQUATORIAL Z WRT MEAN EQU & EQUX OF DATE, AU
C             (REMAINING ELEMENTS OF ARRAY C NOT USED)
C
C
C***********************************************************************
C
C     PART I    TABLES OF THE PERTURBATIONS
C
      DIMENSION X(8,46), X1(80), X2(80), X3(80), X4(80), X5(48)
      EQUIVALENCE (X(1, 1),X1(1))
      EQUIVALENCE (X(1,11),X2(1))
      EQUIVALENCE (X(1,21),X3(1))
      EQUIVALENCE (X(1,31),X4(1))
      EQUIVALENCE (X(1,41),X5(1))
C
C     PERTURBATIONS BY VENUS
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA X1 /  - 1.,  0., +  33.,-  67., -  85.,-  39., +  24.,-  17.,
     2           - 1.,+ 1., +2353.,-4228., -2062.,-1146., -   4.,+   3.,
     3           - 1.,+ 2., -  65.,-  34., +  68.,-  14., +   6.,-  92.,
     4           - 2.,+ 1., -  99.,+  60., +  84.,+ 136., +  23.,-   3.,
     5           - 2.,+ 2., -4702.,+2903., +3593.,+5822., +  10.,-   6.,
     6           - 2.,+ 3., +1795.,-1737., - 596.,- 632., +  37.,-  56.,
     7           - 3.,+ 3., - 666.,+  27., +  44.,+1044., +   8.,+   1.,
     8           - 3.,+ 4., +1508.,- 397., - 381.,-1448., + 185.,- 100.,
     9           - 3.,+ 5., + 763.,- 684., + 126.,+ 148., +   6.,-   3.,
     *           - 4.,+ 4., - 188.,-  93., - 166.,+ 337.,     0.,    0./
      DATA X2 /  - 4.,+ 5., - 139.,-  38., -  51.,+ 189., -  31.,-   1.,
     2           - 4.,+ 6., + 146.,-  42., -  25.,-  91., +  12.,    0.,
     3           - 5.,+ 5., -  47.,-  69., - 134.,+  93.,     0.,    0.,
     4           - 5.,+ 7., - 119.,-  33., -  37.,+ 136., -  18.,-   6.,
     5           - 5.,+ 8., + 154.,    0.,     0.,-  26.,     0.,    0.,
     6           - 6.,+ 6., -   4.,-  38., -  80.,+   8.,     0.,    0.,
C
C     PERTURBATIONS BY MARS
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
     7           + 1.,- 1., - 216.,- 167., -  92.,+ 119.,     0.,    0.,
     8           + 2.,- 2., +1963.,- 567., - 573.,-1976.,     0.,-   8.,
     9           + 2.,- 1., -1659.,- 617., +  64.,- 137.,     0.,    0.,
     *           + 3.,- 3., +  53.,- 118., - 154.,-  67.,     0.,    0./
      DATA X3 /  + 3.,- 2., + 396.,- 153., -  77.,- 201.,     0.,    0.,
     2           + 4.,- 3., - 131.,+ 483., + 461.,+ 125., +   7.,+   1.,
     3           + 4.,- 2., + 526.,- 256., +  43.,+  96.,     0.,    0.,
     4           + 5.,- 4., +  49.,+  69., +  87.,-  62.,     0.,    0.,
     5           + 5.,- 3., -  38.,+ 200., +  87.,+  17.,     0.,    0.,
     6           + 6.,- 4., - 104.,- 113., - 102.,+  94.,     0.,    0.,
     7           + 6.,- 3., -  11.,+ 100., -  27.,-   4.,     0.,    0.,
     8           + 7.,- 4., -  78.,-  72., -  26.,+  28.,     0.,    0.,
     9           + 9.,- 5., +  60.,-  15., -   4.,-  17.,     0.,    0.,
     *           +15.,- 8., + 200.,-  30., -   1.,-   6.,     0.,    0./
C
C     PERTURBATIONS BY JUPITER
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA X4 /  + 1.,- 2., - 155.,-  52., -  78.,+ 193., +   7.,    0.,
     2           + 1.,- 1., -7208.,+  59., +  56.,+7067., -   1.,+  17.,
     3           + 1.,  0., - 307.,-2582., + 227.,-  89., +  16.,    0.,
     4           + 1.,+ 1., +   8.,-  73., +  79.,+   9., +   1.,+  23.,
     5           + 2.,- 3., +  11.,+  68., + 102.,-  17.,     0.,    0.,
     6           + 2.,- 2., + 136.,+2728., +4021.,- 203.,     0.,    0.,
     7           + 2.,- 1., - 537.,+1518., +1376.,+ 486., +  13.,+ 166.,
     8           + 3.,- 3., - 162.,+  27., +  43.,+ 278.,     0.,    0.,
     9           + 3.,- 2., +  71.,+ 551., + 796.,- 104., +   6.,-   1.,
     *           + 3.,- 1., -  31.,+ 208., + 172.,+  26., +   1.,+  18./
      DATA X5 /  + 4.,- 3., -  43.,+   9., +  13.,+  73.,     0.,    0.,
     2           + 4.,- 2., +  17.,+  78., + 110.,-  24.,     0.,    0.,
C
C     PERTURBATIONS BY SATURN
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
     3           + 1.,- 1., -  77.,+ 412., + 422.,+  79., +   1.,+   6.,
     4           + 1.,  0., -   3.,- 320., +   8.,-   1.,     0.,    0.,
     5           + 2.,- 2., +  38.,- 101., - 152.,-  57.,     0.,    0.,
     6           + 2.,- 1., +  45.,- 103., - 103.,-  44.,     0.,    0./
C
C
C***********************************************************************
C
C     PART II   NECESSARY PRELIMINARIES
C
      DATA TWOPI /6.283185307179586D0/
      DATA STR   /206264806.2470964D0/
      DATA RTD   /57.295779513082321D0/
      DATA R     /1296000.0D0/
      TR = 1000.0D0 / STR
C
C     T  = TIME IN JULIAN CENTURIES FROM 1900 JANUARY 0
      T  = (DJ - 2415020.D0)/36525.D0
C
C     TP = TIME IN JULIAN YEARS     FROM 1850 JANUARY 0
      TP = (DJ - 2396758.D0)/365.25D0
C
C     T20= TIME IN JULIAN CENTURIES FROM J2000.0
      T20= (DJ - 2451545.D0)/36525.D0
C
C
C***********************************************************************
C
C     PART III  COMPUTATION OF ELLIPTIC ELEMENTS AND SECULAR TERMS
C
C     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
      EL( 1) = 1.00000030007166D0
      EL( 2) = 0.016708320D0 + (-0.42229D-04 - 0.126D-06 * T20) * T20
      EL( 5) = 1018578.046D0 + (6190.046D0 +
     1                (1.666D0 + 0.012D0 * T20) * T20) * T20
      EL( 5) = EL( 5) * TR
      EL( 9) = 1009677.850D0 + (100.0D0 * R + 2771.27D0 +
     1                1.089D0 * T20) * T20
      EL( 9) = DMOD (EL( 9) * TR, TWOPI)
      EL(10) = 1287099.804D0 + (99.0D0 * R + 1292581.224D0 +
     1                (-0.577D0 - 0.012D0 * T20) * T20) * T20
      EL(10) = DMOD (EL(10) * TR, TWOPI)
      
C     EXPRESSION FOR OBLIQUITY FROM P03 (IAU 2006) PRECESSION       
      EL(13) = 84381.406D0 + (-46.836769D0 +
     1               (-0.0001831D0 + 0.00200340D0 * T20) * T20) * T20
      EL(13) = EL(13) * TR

C     KAPLAN CORRECTION TO SUN'S MEAN LONGITUDE TO FIT DE405 OVER
C     INTERVAL 1800-2200, USING P03 (IAU 2006) PRECESSION
      EL(9) = EL(9)
     1      + ( 0.1320D0 - 0.1355D0 * T20 ) * TR

C
C***********************************************************************
C
C     PART IV   LUNAR TERMS
C
C     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
      EL(14) = 785939.157D0 + (1336.0D0 * R + 1108372.598D0
     1                + (-5.802D0 + 0.019D0 * T20) * T20) * T20
      EL(14) = DMOD (EL(14) * TR, TWOPI)
      EL(17) = 450160.280D0 + (-5.0D0 * R - 482890.539D0 +
     1                (7.455D0 + 0.008D0 * T20) * T20) * T20
      EL(17) = DMOD (EL(17) * TR, TWOPI)
      EL(21) = 300072.424D0 + (11.0D0 * R + 392449.965D0 +
     1                (-37.112D0 - 0.045D0 * T20) * T20) * T20
      EL(21) = DMOD (EL(21) * TR, TWOPI)
C
C     DERIVED ARGUMENTS
      EL(15) = EL(14) - EL(21)
      EL(16) = EL(14) - EL(17)
      EL(15) = DMOD (EL(15),TWOPI)
      EL(16) = DMOD (EL(16),TWOPI)
C     MEAN ELONGATION
      D      = EL(14) - EL(9)
C
C     COMBINATIONS OF ARGUMENTS AND THE PERTURBATIONS
      D = DMOD (D,TWOPI)
      ARG = D
      DL =    +  6469.*SIN(ARG) +  13.*SIN(3.*ARG)
      DR =    + 13390.*COS(ARG) +  30.*COS(3.*ARG)
C
      DBLARG = D + EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +  177.*SIN(ARG)
      DR = DR +  370.*COS(ARG)
C
      DBLARG = D - EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL -  424.*SIN(ARG)
      DR = DR - 1330.*COS(ARG)
C
      DBLARG = 3.D0*D - EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +   39.*SIN(ARG)
      DR = DR +   80.*COS(ARG)
C
      DBLARG = D + EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL -   64.*SIN(ARG)
      DR = DR -  140.*COS(ARG)
C
      DBLARG = D - EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +  172.*SIN(ARG)
      DR = DR +  360.*COS(ARG)
C
      EL(16) = DMOD (EL(16),TWOPI)
      ARG = EL(16)
      DB =    + 576.*SIN(ARG)
C
C
C***********************************************************************
C
C     PART V    COMPUTATION OF PERIODIC PERTURBATIONS
C
C     THE PERTURBING MEAN ANOMALIES
C
      GV  = 0.19984020D+01 + .1021322923D+02*TP
      GM  = 0.19173489D+01 + .3340556174D+01*TP
      GJ  = 0.25836283D+01 + .5296346478D+00*TP
      GS  = 0.49692316D+01 + .2132432808D+00*TP
      GV  = DMOD (GV,TWOPI)
      GM  = DMOD (GM,TWOPI)
      GJ  = DMOD (GJ,TWOPI)
      GS  = DMOD (GS,TWOPI)
C
C
C     MODIFICATION OF FUNDAMENTAL ARGUMENTS
C
C     APPLICATION OF THE JUPITER-SATURN GREAT INEQUALITY
C     TO JUPITER'S MEAN ANOMALY
C
      GJ = GJ + 0.579904067D-02 * DSIN (5.D0*GS - 2.D0*GJ
     1                 + 1.1719644977D0 - 0.397401726D-03*TP)
      GJ = DMOD (GJ,TWOPI)
C
C     LONG PERIOD PERTURBATIONS OF MEAN ANOMALY
C
      ST = T
C                ARGUMENT IS ( 4 MARS - 7 EARTH + 3 VENUS )
      DG = 266.* SIN (0.555015 + 2.076942*ST)
C                ARGUMENT IS ( 3 JUPITER - 8 MARS + 4 EARTH )
     1    + 6400.* SIN (4.035027 + 0.3525565*ST)
C                ARGUMENT IS ( 13 EARTH - 8 VENUS )
     2    + (1882.-16.*ST) * SIN (0.9990265 + 2.622706*ST)
C
C
C     COMPUTATION OF THE EQUATION OF THE CENTER
C
C     FORM PERTURBED MEAN ANOMALY
      EL(10) = DG/STR + EL(10)
      EL(10) = DMOD (EL(10),TWOPI)
      EL(12) =   DSIN(     EL(10)) * (6910057.D0 -(17240.D0+52.D0*T)*T)
     1         + DSIN(2.D0*EL(10)) * (  72338.D0 -    361.D0*T)
     2         + DSIN(3.D0*EL(10)) * (   1054.D0 -      1.D0*T)
C
C     THE UNPERTURBED RADIUS VECTOR
      RO     =                          30570.D0 -    150.D0*T
     1         - DCOS(     EL(10)) * (7274120.D0 - (18140.D0+50.D0*T)*T)
     2         - DCOS(2.D0*EL(10)) * (  91380.D0 -    460.D0*T)
     3         - DCOS(3.D0*EL(10)) * (   1450.D0 -     10.D0*T)
      EL(11) = 10.D0**(RO*1.D-09)
C
C
C     SELECTED PLANETARY PERTURBATIONS FROM NEWCOMB'S THEORY FOLLOW
C
C     PERTURBATIONS BY VENUS
      DO 20 K=1,16
C     ARGUMENT J * VENUS +   I * EARTH
      DBLARG = X(1,K)*GV + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   20 CONTINUE
C
C     PERTURBATIONS BY MARS
      DO 30 K=17,30
C     ARGUMENT  J * MARS +   I * EARTH
      DBLARG = X(1,K)*GM + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   30 CONTINUE
C
C     PERTURBATIONS BY JUPITER
      DO 40 K=31,42
C     ARGUMENT J*JUPITER +   I * EARTH
      DBLARG = X(1,K)*GJ + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   40 CONTINUE
C
C     PERTURBATIONS BY SATURN
      DO 50 K=43,46
C     ARGUMENT J*SATURN  +   I * EARTH
      DBLARG = X(1,K)*GS + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   50 CONTINUE
C
C
C***********************************************************************
C
C     PART VI   COMPUTATION OF ECLIPTIC AND EQUATORIAL COORDINATES
C
      C(1) = EL(11)*10.D0**(DR*1.D-09)
      C(4) = (DL + DG + EL(12))/STR + EL(9)
      C(4) = DMOD (C(4),TWOPI)
      C(5) = DB/STR
      C(2) = C(4)*RTD
      C(3) = C(5)*RTD
      SINO = DSIN (EL(13))
      COSO = DCOS (EL(13))
      SINL = DSIN (C(4))
      COSL = DCOS (C(4))
      SINB = DSIN (C(5))
      COSB = DCOS (C(5))
      C(11) = C(1) * (COSB * COSL)
      C(12) = C(1) * (COSB * SINL * COSO - SINB * SINO)
      C(13) = C(1) * (COSB * SINL * SINO + SINB * COSO)
C
C
C***********************************************************************
C
C
      RETURN
C
      END



      INTEGER FUNCTION IDSS ( NAME )
*
*     THIS FUNCTION RETURNS THE ID NUMBER OF A SOLAR SYSTEM BODY
*     FOR THE VERSION OF SOLSYS (OR SOLSYS-AUXPOS COMBINATION) IN USE.
*
*         NAME   = NAME OF BODY WHOSE ID NUMBER IS DESIRED, E.G.,
*                  'SUN', 'MOON, 'MERCURY', ETC., EXPRESSED AS ALL
*                  UPPER-CASE LETTERS (IN)
*         IDSS   = ID NUMBER OF BODY, FOR USE IN CALLS TO SOLSYS
*                  (FUNCTION VALUE RETURNED)
*
*     NOTE 1: IN THIS VERSION, ONLY THE FIRST THREE LETTERS OF THE
*     BODY'S NAME ARE USED FOR IDENTIFICATION.  ALTERNATIVE VERSIONS
*     MIGHT USE MORE LETTERS.
*
*     NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=1, SINCE SOLSYS 
*     VERSION 3 DOES NOT PROCESS SPLIT JULIAN DATES.    
*
*     NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
*     THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.
*
*
      CHARACTER NAME*(*), NAMEIN*3, NAMES*3
      DIMENSION NAMES(35), IDS(35)

      DATA NAMES / 'SUN', 'EAR', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---'  /
      DATA IDS   /     0,     3,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0  /
      DATA NUM   / 2 /

   3  FORMAT ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', A )

      IDSS = -9999
      NAMEIN = NAME

*     LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
      DO 20 I = 1, NUM
          IF ( NAMEIN .EQ. NAMES(I) ) THEN
              IDSS = IDS(I)
              GO TO 30
          END IF
  20  CONTINUE
  
*     IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES   
      IF ( NAMEIN .EQ. 'JD ' ) THEN
*         IN THIS CASE, SET IDSS=2 IF SOLSYS PROCESSES SPLIT
*         JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE 
          IDSS = 1
          GO TO 30
      END IF    

      WRITE ( *, 3 ) NAME

  30  RETURN

      END
