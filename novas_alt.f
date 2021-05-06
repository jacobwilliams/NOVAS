*  NOVAS FORTRAN VERS F2.0 (1 NOV 98)
*  ALTERNATIVE VERSIONS OF SOME SUBROUTINES

************************************************************************
*                                                                      *
*                              N O V A S                               *
*           NAVAL OBSERVATORY VECTOR ASTROMETRY SUBROUTINES            *
*                                                                      *
*                            G. H. KAPLAN                              *
*                        U.S. NAVAL OBSERVATORY                        *
*                                                                      *
************************************************************************



      SUBROUTINE SUNFLD (POS1,PE,POS2)
C
C     SUBROUTINE SUNFLD VERSION 2.
C     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR THE DEFLECTION
C     OF LIGHT IN THE GRAVITATIONAL FIELD OF THE SUN.
C     THIS VERSION IS A DUMMY.  NO CORRECTION IS APPLIED.
C
C          POS1 = POSITION VECTOR, REFERRED TO ORIGIN AT CENTER OF MASS
C                 OF THE EARTH, COMPONENTS IN AU (IN)
C          PE   = POSITION VECTOR OF CENTER OF MASS OF THE EARTH,
C                 REFERRED TO ORIGIN AT CENTER OF MASS OF
C                 THE SUN, COMPONENTS IN AU (IN)
C          POS2 = POSITION VECTOR, REFERRED TO ORIGIN AT CENTER OF MASS
C                 OF THE EARTH, CORRECTED FOR GRAVITATIONAL DEFLEC-
C                 TION, COMPONENTS IN AU (OUT)
C
C
      DOUBLE PRECISION POS1,PE,POS2
      DIMENSION POS1(3), PE(3), POS2(3)
C
      DO 20 J=1,3
   20 POS2(J) = POS1(J)
C
      RETURN
C
      END



      SUBROUTINE SUNFLD (POS1,PE,POS2)
C
C     SUBROUTINE SUNFLD VERSION 3.
C     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR THE DEFLECTION
C     OF LIGHT IN THE GRAVITATIONAL FIELD OF THE SUN.  SEE
C     MURRAY (1981) MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.  THIS
C     SUBROUTINE VALID FOR BODIES WITHIN THE SOLAR SYSTEM AS WELL AS
C     FOR STARS.
C
C          POS1 = POSITION VECTOR, REFERRED TO ORIGIN AT CENTER OF MASS
C                 OF THE EARTH, COMPONENTS IN AU (IN)
C          PE   = POSITION VECTOR OF CENTER OF MASS OF THE EARTH,
C                 REFERRED TO ORIGIN AT CENTER OF MASS OF
C                 THE SUN, COMPONENTS IN AU (IN)
C          POS2 = POSITION VECTOR, REFERRED TO ORIGIN AT CENTER OF MASS
C                 OF THE EARTH, CORRECTED FOR GRAVITATIONAL DEFLEC-
C                 TION, COMPONENTS IN AU (OUT)
C
C
      DOUBLE PRECISION POS1,PE,POS2,PQ,PMAG,EMAG,QMAG,PHAT,EHAT,QHAT,
     .     MAU,GS,C,PDOTQ,EDOTP,QDOTE,FAC1,FAC2,P2J,
     .     DSQRT
      DIMENSION POS1(3), PE(3), POS2(3), PQ(3),
     .     PHAT(3), EHAT(3), QHAT(3)
C
      DATA MAU / 1.49597870D11 /
C     MAU = NUMBER OF METERS PER AU
      DATA GS / 1.32712438D20 /
C     GS = HELIOCENTRIC GRAVITATIONAL CONSTANT
      DATA C / 299792458.0D0 /
C     C = SEED OF LIGHT
C
C     CONSTRUCT VECTOR PQ BETWEEN SUN AND BODY
      DO 20 J=1,3
   20 PQ(J) = PE(J) + POS1(J)
C
C     COMPUTE VECTOR MAGNITUDES AND UNIT VECTORS
      PMAG = DSQRT (POS1(1)**2 + POS1(2)**2 + POS1(3)**2)
      EMAG = DSQRT (  PE(1)**2 +   PE(2)**2 +   PE(3)**2)
      QMAG = DSQRT (  PQ(1)**2 +   PQ(2)**2 +   PQ(3)**2)
      DO 30 J=1,3
      PHAT(J) = POS1(J) / PMAG
      EHAT(J) =   PE(J) / EMAG
   30 QHAT(J) =   PQ(J) / QMAG
C
C     COMPUTE DOT PRODUCTS OF VECTORS
      PDOTQ = PHAT(1)*QHAT(1) + PHAT(2)*QHAT(2) + PHAT(3)*QHAT(3)
      EDOTP = EHAT(1)*PHAT(1) + EHAT(2)*PHAT(2) + EHAT(3)*PHAT(3)
      QDOTE = QHAT(1)*EHAT(1) + QHAT(2)*EHAT(2) + QHAT(3)*EHAT(3)
C
C     COMPUTE SCALAR FACTORS
      FAC1 = 2.0D0 * GS / (C * C * EMAG * MAU)
      FAC2 = 1.0D0 + QDOTE
C
C     CONSTRUCT CORRECTED POSITION VECTOR POS2
      DO 50 J=1,3
      P2J = PHAT(J) + FAC1 * (PDOTQ*EHAT(J) - EDOTP*QHAT(J)) / FAC2
   50 POS2(J) = P2J * PMAG
C
      RETURN
C
      END



      SUBROUTINE NOD (TIME, LONGIT,OBLIQ)
C
C-----------------------------------------------------------------------
C     SUBROUTINE NOD VERSION 1F.
C
C---PURPOSE: TO PROVIDE FAST EVALUATION OF THE 1980 IAU THEORY OF
C            NUTATION.
C
C---REFERENCES: SEIDELMANN (1982) CELESTIAL MECHANICS 27, 79-106
C                  (1980 IAU THEORY OF NUTATION).
C               COFFEY AND DEPRIT (1980) ASTRONOMY AND ASTROPHYSICS
C                  81, 310-315 (FAST EVALUATION OF FOURIER SERIES).
C
C---INPUT ARGUMENTS:    TIME = TDB TIME IN JULIAN CENTURIES SINCE
C                              J2000.0.
C
C---OUTPUT ARGUMENTS: LONGIT = NUTATION IN LONGITUDE IN ARCSECONDS.
C                      OBLIQ = NUTATION IN OBLIQUITY IN ARCSECONDS.
C
C---COMMON BLOCKS: NONE.
C
C---SUBROUTINES CALLED: NONE.
C
C---NOTES: 1. WARNING: THIS SUBROUTINE CONTAINS COMPUTER-GENERATED
C             CODE.  MODIFY AT YOUR OWN RISK.
C          2. CODE GENERATED ON 11/29/88 16:35:35 AT THE NATIONAL
C             INSTITUTES OF STANDARDS AND TECHNOLOGY (NIST).
C             BRUCE R. MILLER
C             NIST - CENTER FOR COMPUTING AND APPLIED MATHEMATICS
C             ADMIN. A302
C             GAITHERSBURG, MD 20899
C             CONTACT: MILLER@VAX.CAM.NBS.GOV OR (301) 975-2708
C          3. THIS IS NIST ROUTINE 'NUT80-S'.
C
C-----------------------------------------------------------------------
C
      REAL CLNG,CLNGX,COBL,COBLX
C
      INTEGER I,II,I1,I2,IOP,NAV1,NAV2,NAV,LLNG,LLNGX,LOBL,LOBLX
C
      DOUBLE PRECISION TIME,LONGIT,OBLIQ,A,ANGLE,CC,SS,CS,SC,C,S,LNG,
     . LNGX,OBL,OBLX
C
      DIMENSION A(5),C(106),S(106),NAV1(10),NAV2(10),NAV(183),CLNG(106
     . ),LLNG(106),CLNGX(14),LLNGX(14),COBL(64),LOBL(64),COBLX(8),
     . LOBLX(8)
C
      DATA NAV1/  1,  1,  2,  1,  3,  2,  4,  1,  5,  1/
      DATA NAV2/  1,  1,  1,  6,  2,  2,  4,  4,  5,  5/
      DATA NAV/  3,  1,  2,  2,  6,  3,  3,  1,  3,  2,  1,  4,  3,  6
     . ,  9,  2, 18,  9,  2, 19,  1,  3,  1,  9,  1,  2,  4,  3,  2,
     .   9,  1, 18,  2,  2, 16,  2,  3, 22,  2,  2,  3,  9,  3,  1,
     .  30,  2, 22,  3,  3,  2, 30,  3,  1, 10,  3,  6,  5,  3,  1,
     .   5,  1,  2, 10,  3,  2,  5,  1,  3, 10,  3,  3,  5,  2, 15,
     .  45,  3,  1, 46,  3,  6, 45,  3, 51,  1,  2, 37,  3,  3,  6,
     .  46,  2, 38,  3,  3,  2, 46,  3,  2, 45,  3, 54,  2,  3,  9,
     .   5,  2, 41,  4,  3, 18,  5,  3,  1, 65,  2, 40,  9,  3, 28,
     .   5,  2, 51, 19,  2, 22, 48,  3, 45,  4,  3, 45,  9,  3, 46,
     .   9,  2, 47,  9,  1, 68,  3,  2,  6, 75,  2,  1, 75,  3, 51,
     .   9,  2,  6, 79,  3, 18, 54,  3, 54,  9,  3,  1, 81,  3,  1,
     .  82,  1,  8, 80,  2,  8, 82,  3,  2, 82,  3, 25, 45,  2,  2,
     .  80,  3, 28, 45/
      DATA CLNG/     +1.,     +1.,     -1.,     -1.,     +1.,     -1.,
     .      -1.,     -1.,     -1.,     -1.,     -1.,     +1.,     -1.
     . ,     +1.,     -1.,     +1.,     +1.,     -1.,     -1.,
     .      +1.,     +1.,     -1.,     +1.,     -1.,     +1.,     -1.
     . ,     -1.,     -1.,     +1.,     -1.,     -1.,     +1.,
     .      -1.,     +1.,     +2.,     +2.,     +2.,     +2.,     +2.
     . ,     -2.,     +2.,     +2.,     +2.,     +3.,     -3.,
     .      -3.,     +3.,     -3.,     +3.,     -3.,     +3.,     +4.
     . ,     +4.,     -4.,     -4.,     +4.,     -4.,     +5.,
     .      +5.,     +5.,     -5.,     +6.,     +6.,     +6.,     -6.
     . ,     +6.,     -7.,     +7.,     +7.,     -7.,     -8.,
     .     +10.,    +11.,    +12.,    -13.,    -15.,    -16.,    -16.
     . ,    +17.,    -21.,    -22.,    +26.,    +29.,    +29.,
     .     -31.,    -38.,    -46.,    +48.,    -51.,    +58.,    +59.
     . ,    +63.,    +63.,   -123.,   +129.,   -158.,   -217.,
     .    -301.,   -386.,   -517.,   +712.,  +1426.,  +2062.,  -2274.
     . , -13187.,-171996./
      DATA LLNG/ 58, 26, 83, 35, 42, 67, 34, 37, 20, 89, 19,105, 94,
     .  85, 48, 29, 84, 87, 70, 76, 90, 31, 59, 74, 47, 78, 24, 33,
     .  60, 73, 32, 17, 75, 23, 99, 39, 63, 97, 38, 36,  7, 77, 86,
     .  52, 27, 11, 14, 64,106, 53,103, 68,100, 16, 25, 15,  4,101,
     .  66, 12, 56, 69, 21, 88, 65, 96, 28, 61, 62, 81, 92, 95, 13,
     .  44, 72, 43, 98, 71,  8, 50, 30,  3,  6, 93, 51, 79, 57, 18,
     .  49, 41, 91,  9, 40, 55, 82, 22,104, 54, 46,102,  1,  2, 10,
     .  45, 80,  5/
      DATA CLNGX/    +0.1,    -0.1,    +0.1,    +0.1,    +0.1,    +0.1
     . ,    +0.2,    -0.2,    -0.4,    +0.5,    +1.2,    -1.6,
     .     -3.4,  -174.2/
      DATA LLNGX/ 82,  8, 98,  1, 40, 41, 10, 45, 46,104,102, 80,  2,
     .   5/
      DATA COBL/     +1.,     +1.,     +1.,     -1.,     -1.,     -1.,
     .      +1.,     +1.,     +1.,     +1.,     +1.,     -1.,     +1.
     . ,     -1.,     +1.,     -1.,     -1.,     -1.,     +1.,
     .      -1.,     +1.,     +1.,     -1.,     -2.,     -2.,     -2.
     . ,     +3.,     +3.,     -3.,     +3.,     +3.,     -3.,
     .      +3.,     +3.,     -3.,     +3.,     +3.,     +5.,     +6.
     . ,     +7.,     -7.,     +7.,     -8.,     +9.,    -10.,
     .     -12.,    +13.,    +16.,    -24.,    +26.,    +27.,    +32.
     . ,    -33.,    -53.,    +54.,    -70.,    -95.,   +129.,
     .    +200.,   +224.,   -895.,   +977.,  +5736., +92025./
      DATA LOBL/ 52, 99, 18, 22,  6,  3, 64,106, 39, 53,103, 63, 97,
     .  38, 36, 77, 37, 89, 86,105, 94, 85, 84, 68,100,  9, 69,101,
     .  61, 62, 92, 88, 65, 81, 96, 66, 56, 95, 44, 98,  1, 72, 71,
     .  43, 50, 93, 51, 79, 57, 91, 49, 41, 40, 55,  2, 82,104, 54,
     .  46,102, 10, 45, 80,  5/
      DATA COBLX/    -0.1,    -0.1,    +0.3,    +0.5,    -0.5,    -0.6
     . ,    -3.1,    +8.9/
      DATA LOBLX/ 54,  2,104, 10, 45,102, 80,  5/
C
C     Calculation of Angles: (L LP F D W)
      CALL FUNARG (TIME,A(1),A(2),A(3),A(4),A(5))
C
C     Call to Library routines for Initial Sine and Cosine
      I = 1
      DO 1 II = 1,10,2
        ANGLE = A(NAV1(II))*DFLOAT(NAV1(1+II))
        C(I) = DCOS(ANGLE)
        S(I) = DSIN(ANGLE)
        I = I+1
   1  CONTINUE
C     Loop to compute trigs using operation ADD
      I = 6
      DO 2 II = 1,10,2
        I1 = NAV2(II)
        I2 = NAV2(1+II)
C       Addition of Sines and Cosines
        C(I) = C(I1)*C(I2)-S(I1)*S(I2)
        S(I) = S(I1)*C(I2)+C(I1)*S(I2)
        I = I+1
   2  CONTINUE
C     Loop to compute trigs using operations (ADD SUB +/-)
      I = 11
      DO 3 II = 1,183,3
        IOP = NAV(II)
        I1 = NAV(1+II)
        I2 = NAV(2+II)
        GO TO (4,5,6,7)IOP
   4    CONTINUE
C         Addition of Sines and Cosines
          C(I) = C(I1)*C(I2)-S(I1)*S(I2)
          S(I) = S(I1)*C(I2)+C(I1)*S(I2)
          I = I+1
        GO TO 7
   5    CONTINUE
C         Subtraction of Sines and Cosines
          C(I) = C(I1)*C(I2)+S(I1)*S(I2)
          S(I) = S(I1)*C(I2)-C(I1)*S(I2)
          I = I+1
        GO TO 7
   6    CONTINUE
C         Addition and Subtraction of Sines and Cosines
          CC = C(I1)*C(I2)
          SS = S(I1)*S(I2)
          SC = S(I1)*C(I2)
          CS = C(I1)*S(I2)
          C(I) = CC-SS
          S(I) = SC+CS
          I = I+1
          C(I) = CC+SS
          S(I) = SC-CS
          I = I+1
   7    CONTINUE
   3  CONTINUE
      LNG = 0.0D0
      DO 8 I = 1,106,1
        LNG = LNG+DBLE(CLNG(I))*S(LLNG(I))
   8  CONTINUE
      LNGX = 0.0D0
      DO 9 I = 1,14,1
        LNGX = LNGX+DBLE(CLNGX(I))*S(LLNGX(I))
   9  CONTINUE
      OBL = 0.0D0
      DO 10 I = 1,64,1
        OBL = OBL+DBLE(COBL(I))*C(LOBL(I))
  10  CONTINUE
      OBLX = 0.0D0
      DO 11 I = 1,8,1
        OBLX = OBLX+DBLE(COBLX(I))*C(LOBLX(I))
  11  CONTINUE
      LONGIT = (LNG+TIME*LNGX)/10000.0D0
      OBLIQ = (OBL+TIME*OBLX)/10000.0D0
C
C     End of Subroutine NOD
      RETURN
C
      END



      SUBROUTINE NOD (T,DPSI,DEPS)
C
C     SUBROUTINE NOD VERSION 2.
C     THIS SUBROUTINE EVALUATES A TRUNCATED NUTATION SERIES AND
C     RETURNS APPROXIMATE VALUES FOR NUTATION IN LONGITUDE AND
C     NUTATION IN OBLIQUITY.
C     LARGEST FOUR TERMS OF 1980 IAU THEORY OF NUTATION.
C
C          T    = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
C          DPSI = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
C          DEPS = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
C
C
      DOUBLE PRECISION T,DPSI,DEPS,L,LP,F,D,OM,ARG,DBLE,DSIN,DCOS
      DIMENSION X(9,4)
C
      DATA X / 0.,  0.,  0.,  0.,  1., -171996., -174.2,  92025.,  8.9,
     .         0.,  0.,  2., -2.,  2.,  -13187.,   -1.6,   5736., -3.1,
     .         0.,  0.,  2.,  0.,  2.,   -2274.,   -0.2,    977., -0.5,
     .         0.,  0.,  0.,  0.,  2.,    2062.,    0.2,   -895.,  0.5/
C
C     COMPUTATION OF FUNDAMENTAL ARGUMENTS
      CALL FUNARG (T,L,LP,F,D,OM)
C
      DPSI = 0.D0
      DEPS = 0.D0
C
C     SUM NUTATION SERIES TERMS
      DO 10 J=1,4
      I = 5 - J
      ARG = DBLE(X(3,I)) * F
     .    + DBLE(X(4,I)) * D
     .    + DBLE(X(5,I)) * OM
      DPSI = (DBLE(X(6,I)) + DBLE(X(7,I))*T) * DSIN(ARG) + DPSI
      DEPS = (DBLE(X(8,I)) + DBLE(X(9,I))*T) * DCOS(ARG) + DEPS
   10 CONTINUE
C
      DPSI = DPSI * 1.0D-4
      DEPS = DEPS * 1.0D-4
C
      RETURN
C
      END



      SUBROUTINE NOD (T,DPSI,DEPS)
C
C     SUBROUTINE NOD VERSION 3.
C     THIS SUBROUTINE EVALUATES THE NUTATION SERIES AND RETURNS THE
C     VALUES FOR NUTATION IN LONGITUDE AND NUTATION IN OBLIQUITY.
C     IERS (1996) THEORY OF PRECESSION/NUTATION (OR ITS SUCCESSORS)
C     AS DESCRIBED IN IERS CONVENTIONS (1996), PP. 25 FF.
C
C          T    = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
C          DPSI = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
C          DEPS = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
C
C
      DOUBLE PRECISION T,DPSI,DEPS,
     .     T0,TJD,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10
C
      DATA T0 / 2451545.00000000D0 /
C     T0 = TDB JULIAN DATE OF EPOCH J2000.0
C
      TJD = T0 + T * 36525.D0
C
C     CALL IERS SUBROUTINE
C     OUTPUT TO BE USED WITH IAU (1976) PRECESSION
      CALL KSV_1996_3 ( TJD, X1,X2,X3,X4,X5,X6,X7,X8,X9,X10 )
C
      DPSI = X9  * 1.D-3
      DEPS = X10 * 1.D-3
C
      RETURN
C
      END



*** SOLSYS VERSION 1 PACKAGE: SOLSYS, BLOCK DATA, FILDEF ***

      SUBROUTINE SOLSYS (TJD,M,K,POS,VEL,IERR)
C
C     SUBROUTINE SOLSYS VERSION 1.
C     THIS SUBROUTINE READS A COORDINATE FILE CONTAINING HELIOCENTRIC
C     POSITIONS OF SOLAR SYSTEM BODIES AT DAILY INTERVALS AND PROVIDES
C     THE POSITION AND VELOCITY OF BODY M AT EPOCH TJD.
C
C          TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
C          M    = BODY IDENTIFICATION NUMBER (IN)
C          K    = ORIGIN SELECTION CODE (IN)
C                 SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
C                 SET K=1 FOR ORIGIN AT CENTER OF SUN
C          POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
C                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
C                 OF J2000.0, COMPONENTS IN AU (OUT)
C          VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
C                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
C                 OF J2000.0, COMPONENTS IN AU/DAY (OUT)
C          IERR = ERROR INDICATOR (OUT)
C                 IERR=0 MEANS EVERYTHING OK
C                 IERR=1 MEANS TJD BEFORE FIRST USABLE DATE IN FILE
C                 IERR=2 MEANS TJD AFTER LAST USABLE DATE IN FILE
C                 IERR=3 MEANS BAD VALUE OF M
C
C
      DOUBLE PRECISION TJD,POS,VEL,RMASS,XJD,XYZ,BPOS,BVEL,
     .     DS,DT,DF,DL,DSLAST,TMASS,ASTART,SSBARY,POSN,T,AK,AI,P,
     .     DABS,DFLOAT
      DIMENSION POS(3), VEL(3), RMASS(20), XJD(11), XYZ(11,20,3),
     .     BPOS(11,3), BVEL(11,3)
      COMMON /SSFILE/ LU,N,RMASS
C     COMMON BLOCK SSFILE CONTAINS INFORMATION ON THE COORDINATE FILE
      SAVE
C
      DATA DF,DL,DSLAST,MLAST,KLAST / 0.0D0,1.0D10,0.0D0,0,0 /
C
      IERR = 0
      IF (DSLAST.GT.0.0D0) GO TO 11
      REWIND LU
      READ (LU) DT
      REWIND LU
      DSLAST = 1.0D0
      DF = DT + 15.0D0
      TMASS = 1.0D0
      DO 10 L=1,N
   10 TMASS = TMASS + 1.0D0 / RMASS(L)
      INTPTS = 5
      LMIDDL = 6
      LSTART = LMIDDL - INTPTS/2 - 1
      ASTART = LSTART
C
C     LOGIC TO DETERMINE BEST WAY TO SEARCH COORDINATE FILE
   11 IF (M.LT.0.OR.M.GT.N) GO TO 75
      IF (TJD.LE.DF) GO TO 76
      IF (TJD.GE.DL) GO TO 77
      IF (DABS(TJD-DSLAST).LT.0.8D0) GO TO 12
      IS = TJD
      DS = DFLOAT(IS) + 0.5D0
      IF (DS-DSLAST) 14,12,15
   12 IF ( M- MLAST) 25,13,25
   13 IF ( K- KLAST) 25,60,25
   14 REWIND LU
      DSLAST = 1.0D0
   15 IF (DS-DSLAST.LT.20.0D0) GO TO 20
C
C     COURSE SEARCH THROUGH COORDINATE FILE
   16 READ (LU,END=77) DT
      IF (DS-DT    .GT.20.0D0) GO TO 16
      DO 18 I=1,11
   18 READ (LU,END=77) XJD(I), ((XYZ(I,L,J), J=1,3), L=1,N)
C
C     FINE SEARCH THROUGH COORDINATE FILE
   20 DO 22 I=1,10
      IOLD = I + 1
      XJD(I) = XJD(IOLD)
      DO 22 L=1,N
      DO 22 J=1,3
   22 XYZ(I,L,J) = XYZ(IOLD,L,J)
      READ (LU,END=77) XJD(11), ((XYZ(11,L,J), J=1,3), L=1,N)
      IF (DABS(DS-XJD(6)).GT.1.0D-6) GO TO 20
      DSLAST = DS
C
C     FILL ARRAY BPOS WITH DAILY POSITIONS OF BODY M
C     IF K=0, MOVE ORIGIN TO SOLAR SYSTEM BARYCENTER
   25 DO 40 J=1,3
      DO 40 I=1,11
      SSBARY = 0.0D0
      IF (K.EQ.1) GO TO 35
      DO 30 L=1,N
   30 SSBARY = SSBARY + XYZ(I,L,J) / (TMASS * RMASS(L))
   35 POSN = 0.0D0
      IF (M.GT.0) POSN = XYZ(I,M,J)
      BPOS(I,J) = POSN - SSBARY
   40 CONTINUE
C
C     FILL ARRAY BVEL WITH DAILY VELOCITIES OF BODY M
C     COMPUTED FROM NUMERICAL DIFFERENTIATION OF POSITIONS IN ARRAY BPOS
      DO 50 J=1,3
      DO 50 I=1,11
      BVEL(I,J) = 0.0D0
      IF (I.LT.4.OR.I.GT.8) GO TO 50
      BVEL(I,J) = (           BPOS(I+3,J) -  9.0D0 * BPOS(I+2,J)
     .             + 45.0D0 * BPOS(I+1,J) - 45.0D0 * BPOS(I-1,J)
     .             +  9.0D0 * BPOS(I-2,J) -          BPOS(I-3,J))
     .             / 60.0D0
   50 CONTINUE
C
      MLAST = M
      KLAST = K
C
C     PERFORM LAGRANGIAN INTERPOLATION FOR POSITION AND VELOCITY AT
C     EPOCH TJD
   60 T = TJD - XJD(6) + 6.0D0
      DO 73 J=1,3
      POS(J) = 0.0D0
      VEL(J) = 0.0D0
      DO 72 L=1,INTPTS
      AK = ASTART + DFLOAT(L)
      P = 1.0D0
      DO 71 I=1,INTPTS
      IF (I.EQ.L) GO TO 71
      AI = ASTART + DFLOAT(I)
      P = P * (T-AI) / (AK-AI)
   71 CONTINUE
      POS(J) = POS(J) + P * BPOS(LSTART+L,J)
      VEL(J) = VEL(J) + P * BVEL(LSTART+L,J)
   72 CONTINUE
   73 CONTINUE
      RETURN
C
   75 IERR = 3
      GO TO 88
   76 IERR = 1
      GO TO 80
   77 IERR = 2
      IF (TJD.LT.DL) DL = TJD
   80 REWIND LU
      DSLAST = 1.0D0
   88 RETURN
C
      END



      BLOCK DATA
C
C     FOR USE WITH SUBROUTINE SOLSYS VERSION 1.
C     COMMON BLOCK SSFILE CONTAINS INFORMATION ON THE COORDINATE FILE
C     USED BY SUBROUTINE SOLSYS.  THIS BLOCK DATA SEGMENT SETS UP THE
C     DEFAULT VALUES FOR THE PARAMETERS IN SSFILE.  THESE DEFAULTS CAN
C     BE ALTERED BY EXECUTABLE STATEMENTS IN THE MAIN PROGRAM OR ANY
C     SUBROUTINE, OR BY A CALL TO SUBROUTINE FILDEF.
C
C          LU    = FORTRAN LOGICAL UNIT NUMBER OF COORDINATE FILE
C          N     = NUMBER OF BODIES WITH COORDINATES IN FILE
C          RMASS = ARRAY OF RECIPROCAL MASSES OF BODIES WITH
C                  COORDINATES IN FILE, IN SOLAR MASS UNITS
C
C
      DOUBLE PRECISION RMASS
      COMMON /SSFILE/ LU,N,RMASS(20)
      DATA LU,N / 20,9 /
C     RECIPROCAL PLANETARY MASSES FROM JPL DE-405
      DATA RMASS / 6023600.0D0,408523.71D0,328900.56D0,3098708.0D0,
     .     1047.3486D0,3497.898D0,22902.98D0,19412.24D0,135200000.0D0,
     .     11*1.0D50 /
      END



      SUBROUTINE FILDEF (LLU,NN,RRMASS)
C
C     FOR USE WITH SUBROUTINE SOLSYS VERSION 1.
C     THIS SUBROUTINE MAY BE CALLED TO CHANGE THE VALUES IN
C     COMMON BLOCK SSFILE, WHICH CONTAINS INFORMATION ON THE
C     COORDINATE FILE USED BY SUBROUTINE SOLSYS.
C
C          LLU    = FORTRAN LOGICAL UNIT NUMBER TO BE USED FOR
C                   COORDINATE FILE (IN)
C          NN     = NUMBER OF BODIES WITH COORDINATES IN FILE (IN)
C          RRMASS = ARRAY OF RECIPROCAL MASSES OF BODIES WITH
C                   COORDINATES IN FILE, IN SOLAR MASS UNITS (IN)
C
C
      DOUBLE PRECISION RRMASS,RMASS
      DIMENSION RRMASS(20)
      COMMON /SSFILE/ LU,N,RMASS(20)
      SAVE
C
   10 IF (LLU.LT.1.OR.LLU.GT.99) GO TO 20
      LU = LLU
   20 IF (NN.LT.1.OR.NN.GT.20) GO TO 50
      N = NN
   30 IF (RRMASS(1).LT.0.01D0) GO TO 50
      DO 35 L=1,NN
      IF (RRMASS(L).LT.0.01D0) GO TO 35
      RMASS(L) = RRMASS(L)
   35 CONTINUE
C
   50 RETURN
C
      END



*** SOLSYS VERSION 3 PACKAGE: SOLSYS, SUN ***

      SUBROUTINE SOLSYS (TJD,M,K,POS,VEL,IERR)
C
C     SUBROUTINE SOLSYS VERSION 3.
C     THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF THE
C     EARTH AT EPOCH TJD BY EVALUATING A CLOSED-FORM THEORY WITHOUT
C     REFERENCE TO AN EXTERNAL FILE.  THIS ROUTINE CAN ALSO PROVIDE
C     THE POSITION AND VELOCITY OF THE SUN.
C
C          TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
C          M    = BODY IDENTIFICATION NUMBER (IN)
C                 SET M=0 OR M=1 OR M=10 FOR THE SUN
C                 SET M=2 OR M=3 FOR THE EARTH
C          K    = ORIGIN SELECTION CODE (IN)
C                 SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
C                 SET K=1 FOR ORIGIN AT CENTER OF SUN
C          POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
C                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
C                 OF J2000.0, COMPONENTS IN AU (OUT)
C          VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
C                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
C                 OF J2000.0, COMPONENTS IN AU/DAY (OUT)
C          IERR = ERROR INDICATOR (OUT)
C                 IERR=0 MEANS EVERYTHING OK
C                 IERR=1 MEANS TJD BEFORE FIRST ALLOWED DATE
C                 IERR=2 MEANS TJD AFTER LAST ALLOWED DATE
C
C
      DOUBLE PRECISION TJD,POS,VEL,
     .     PM,PA,PL,PN,TWOPI,TLAST,T0,OBL,SINE,COSE,TMASS,
     .     QJD,EL,C,P,
     .     F,PBARY,VBARY,DLON,SINL,COSL,X,Y,Z,XDOT,YDOT,ZDOT,
     .     DFLOAT,DABS,DMOD,DSIN,DCOS
      DIMENSION POS(3), VEL(3), EL(21), C(13), P(3,3),
     .     PM(4), PA(4), PL(4), PN(4),
     .     PBARY(3), VBARY(3)
      SAVE
C
      DATA EL,C / 34*0.0D0 /
C
C     ARRAYS BELOW CONTAIN DATA ON THE FOUR LARGEST PLANETS (SEE
C     EXPLANATORY SUPPLEMENT, P. 316)
C     THIS DATA USED FOR BARYCENTER COMPUTATIONS ONLY
C                 JUPITER        SATURN        URANUS       NEPTUNE
      DATA PM /  1047.349D 0,  3497.898D 0,   22903.0D 0,   19412.2D 0 /
      DATA PA /  5.203363D 0,  9.537070D 0, 19.191264D 0, 30.068963D 0 /
      DATA PL /  0.600470D 0,  0.871693D 0,  5.466933D 0,  5.321160D 0 /
      DATA PN /  1.450138D-3,  5.841727D-4,  2.047497D-4,  1.043891D-4 /
C
      DATA TWOPI / 6.283185307179586D0 /,     TLAST / 0.0D0 /
      DATA T0 / 2451545.00000000D0 /,     OBL / 23.43929111D0 /
C     T0 = TDB JULIAN DATE OF EPOCH J2000.0
C     OBL = OBLIQUITY OF ECLIPTIC AT EPOCH J2000.0
C
      IF (TLAST.GT.0.0D0) GO TO 12
      SINE = DSIN (OBL * TWOPI/360.0D0)
      COSE = DCOS (OBL * TWOPI/360.0D0)
      TMASS = 1.0D0
      DO 10 I=1,4
   10 TMASS = TMASS + 1.0D0 / PM(I)
      TLAST = 1.0D0
   12 IERR = 0
      IF (TJD.LT.2340000.5D0) IERR = 1
      IF (TJD.GT.2560000.5D0) IERR = 2
      IF (IERR.NE.0) GO TO 110
      IF (M.GE.10) GO TO 20
      IF (M.GE. 2) GO TO 30
C
C     FORM HELIOCENTRIC COORDINATES OF SUN
   20 DO 25 J=1,3
      POS(J) = 0.0D0
   25 VEL(J) = 0.0D0
      IF (K) 90,90,110
C
C     FORM HELIOCENTRIC COORDINATES OF EARTH
   30 DO 35 I=1,3
      QJD = TJD + DFLOAT(I-2) * 0.10D0
C     SUBROUTINE SUN COMPUTES EARTH-SUN VECTOR
      CALL SUN (QJD,EL,C)
      CALL PRECES (QJD,C(11),T0,POS)
      P(I,1) = -POS(1)
      P(I,2) = -POS(2)
      P(I,3) = -POS(3)
   35 CONTINUE
      DO 40 J=1,3
      POS(J) =  P(2,J)
      VEL(J) = (P(3,J) - P(1,J)) / 0.20D0
   40 CONTINUE
      IF (K) 90,90,110
C
C     IF K=0, MOVE ORIGIN TO SOLAR SYSTEM BARYCENTER
C     SOLAR SYSTEM BARYCENTER COORDINATES COMPUTED FROM ROUGH
C     APPROXIMATIONS OF THE COORDINATES OF THE FOUR LARGEST PLANETS
   90 IF (DABS(TJD-TLAST).LT.1.0D-6) GO TO 99
      DO 92 J=1,3
      PBARY(J) = 0.0D0
   92 VBARY(J) = 0.0D0
C     THE FOLLOWING LOOP CYCLES ONCE FOR EACH OF THE FOUR PLANETS
      DO 98 I=1,4
      DLON = PL(I) + PN(I) * (TJD - T0)
      DLON = DMOD(DLON,TWOPI)
      SINL = DSIN(DLON)
      COSL = DCOS(DLON)
C     SINL AND COSL ARE THE SINE AND COSINE OF PLANET'S MEAN LONGITUDE
      X    =  PA(I) * COSL
      Y    =  PA(I) * SINL * COSE
      Z    =  PA(I) * SINL * SINE
      XDOT = -PA(I) * PN(I) * SINL
      YDOT =  PA(I) * PN(I) * COSL * COSE
      ZDOT =  PA(I) * PN(I) * COSL * SINE
      F = 1.0D0 / (PM(I) * TMASS)
      PBARY(1) = PBARY(1) + X * F
      PBARY(2) = PBARY(2) + Y * F
      PBARY(3) = PBARY(3) + Z * F
      VBARY(1) = VBARY(1) + XDOT * F
      VBARY(2) = VBARY(2) + YDOT * F
      VBARY(3) = VBARY(3) + ZDOT * F
   98 CONTINUE
      TLAST = TJD
   99 DO 100 J=1,3
      POS(J) = POS(J) - PBARY(J)
  100 VEL(J) = VEL(J) - VBARY(J)
C
  110 RETURN
C
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
C     ARE USED.  THE ABSOLUTE ACCURACY IS OF ORDER 1 ARCSECOND
C     FOR EPOCHS NEAR THE YEAR 2000.
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
      EL(13) = 84381.448D0 + (-46.8150D0 +
     1               (-0.00059D0 + 0.001813D0 * T20) * T20) * T20
      EL(13) = EL(13) * TR
C
C     KAPLAN CORRECTION TO SUN'S MEAN LONGITUDE TO FIT DE405 OVER
C     INTERVAL 1800-2200
      EL(9) = EL(9)
     1      + ( 0.1402D0 + 0.1645D0 * T20 ) * TR
C
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
