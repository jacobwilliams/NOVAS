*  NOVAS FORTRAN VERS F3.1
*  ALTERNATIVE VERSIONS OF SOME SUBROUTINES

************************************************************************
*                                                                      *
*                              N O V A S                               *
*           NAVAL OBSERVATORY VECTOR ASTROMETRY SOFTWARE               *
*                                                                      *
*                            G. H. KAPLAN                              *
*                        U.S. NAVAL OBSERVATORY                        *
*                                                                      *
************************************************************************



      SUBROUTINE GRVDEF (TJD,LOC,POS1,POBS,POS2)
*
*     SUBROUTINE GRVDEF VERSION 2.
*     THIS SUBROUTINE COMPUTES THE TOTAL GRAVITATIONAL DEFLECTION OF
*     LIGHT FOR THE OBSERVED OBJECT DUE TO THE MAJOR GRAVITATING BODIES
*     IN THE SOLAR SYSTEM.
*     THIS VERSION IS A DUMMY.  NO CORRECTION IS APPLIED.
*
*          TJD    = (NOT USED)
*          LOC    = (NOT USED)
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
*                   TO ICRS AXES, COMPONENTS IN AU (IN)
*          POBS   = (NOT USED)
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
*                   TO ICRS AXES, CORRECTED FOR GRAVITATIONAL
*                   DEFLECTION, COMPONENTS IN AU (OUT)
*
*     NOTE:  IN THIS VERSION, POS2 = POS1.  THE GRAVITATIONAL DEFLECTION
*     THAT IS NEGLECTED HERE CAN REACH 1.8 ARCSECONDS AT THE LIMB OF 
*     THE SUN, BUT IS LESS THAN 10 MILLIARCSECONDS OVER THE AREA OF SKY  
*     MORE THAN 45 DEGREES FROM THE SUN.  SEE TABLE 3.26.1 ON PAGE 138
*     OF THE EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC (1992).
*
*
      DOUBLE PRECISION TJD,POS1,POBS,POS2
      DIMENSION POS1(3), POBS(3), POS2(3)

      DO 20 J = 1, 3
  20  POS2(J) = POS1(J)

      RETURN

      END



      SUBROUTINE NOD (T,DPSI,DEPS)
*
*     SUBROUTINE NOD VERSION 2.
*     IN LOW-ACCURACY MODE, THIS SUBROUTINE EVALUATES A SHORT
*     NUTATION SERIES AND RETURNS APPROXIMATE VALUES FOR NUTATION IN
*     LONGITUDE AND NUTATION IN OBLIQUITY FOR A GIVEN TDB JULIAN DATE.
*     IN THIS MODE, ONLY THE LARGEST 13 TERMS OF THE IAU 2000A NUTATION
*     SERIES ARE EVALUATED.  IN HIGH-ACCURACY MODE, THE STANDARD IERS
*     SUBROUINE IS CALLED TO EVALUATE THE FULL IAU 2000A NUTATION
*     SERIES.
*
*          T    = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          DPSI = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
*          DEPS = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
*
*     NOTE:  IN LOW-ACCURACY MODE, MAX ERROR IN DPSI < 0.05 ARCSEC,
*     MAX ERROR IN DEPS < 0.02 ARCSEC, AVERAGE ERROR ABOUT 1/4 OF MAX.
*
*
      DOUBLE PRECISION T,DPSI,DEPS,PI,SECCON,T0,T1,DP,DE,
     .     X,EL,ELP,F,D,OM,ARG,DSIN,DCOS
      DIMENSION X(9,13)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

*     LARGEST 13 TERMS OF IAU 2000A NUTATION SERIES, WITH PRECISION
*     OF COEFFICIENTS TRUNCATED
      DATA X / 0., 0., 0., 0., 1., -17.2064,-0.01747, 9.2052, 0.00091,
     .         0., 0., 2.,-2., 2.,  -1.3171,-0.00017, 0.5730,-0.00030,
     .         0., 0., 2., 0., 2.,  -0.2276,-0.00002, 0.0978,-0.00005,
     .         0., 0., 0., 0., 2.,   0.2075, 0.00002,-0.0897, 0.00005,
     .         0., 1., 0., 0., 0.,   0.1476,-0.00036, 0.0074,-0.00002,
     .         0., 1., 2.,-2., 2.,  -0.0517, 0.00012, 0.0224,-0.00007,
     .         1., 0., 0., 0., 0.,   0.0711, 0.00001,-0.0007, 0.00000,
     .         0., 0., 2., 0., 1.,  -0.0387,-0.00004, 0.0201, 0.00000,
     .         1., 0., 2., 0., 2.,  -0.0301, 0.00000, 0.0129,-0.00001,
     .         0.,-1., 2.,-2., 2.,   0.0216,-0.00005,-0.0096, 0.00003,
     .         0., 0., 2.,-2., 1.,   0.0128, 0.00001,-0.0069,-0.00000,
     .        -1., 0., 2., 0., 2.,   0.0123, 0.00000,-0.0053, 0.00000,
     .        -1., 0., 0., 2., 0.,   0.0157, 0.00000,-0.0001, 0.00000 / 
*     REMAINING TERMS ALL HAVE AMPLITUDES < 0.01 ARCSECOND

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN

*         HIGH ACCURACY MODE -- CALL IERS SUBROUTINE

          T1 = T * 36525.D0
          CALL NU2000A ( T0, T1, DP, DE )
          DPSI = DP * SECCON
          DEPS = DE * SECCON

      ELSE 

*         LOW ACCURACY MODE -- EVALUATE SHORT NUTATION SERIES ABOVE

*         COMPUTATION OF FUNDAMENTAL ARGUMENTS
          CALL FUNARG ( T,   EL, ELP, F, D, OM )

          DPSI = 0.D0
          DEPS = 0.D0

*         SUM NUTATION SERIES TERMS
          DO 10 I = 13, 1, -1
              ARG = X(1,I) * EL
     .            + X(2,I) * ELP
     .            + X(3,I) * F
     .            + X(4,I) * D
     .            + X(5,I) * OM
              DPSI = ( X(6,I) + X(7,I) * T ) * DSIN ( ARG ) + DPSI
              DEPS = ( X(8,I) + X(9,I) * T ) * DCOS ( ARG ) + DEPS
  10      CONTINUE
  
*         ADD IN OUT-OF-PHASE COMPONENT OF PRINCIPAL (18.6-YEAR) TERM
*         (TO AVOID SMALL BUT LONG-TERM BIAS IN RESULTS)
          DPSI = DPSI + 0.0033D0 * DCOS ( OM )
          DEPS = DEPS + 0.0015D0 * DSIN ( OM )
          
      END IF    
  
      RETURN

      END
