*
*     NOVAS SAMPLE CALCULATIONS
*
*     SEE CHAPTER 3 OF USER'S GUIDE FOR EXPLANATION
*     REQUIRES SOLSYS VERSION 2 
*     FOR USE WITH SOLSYS VERSION 1, SEE COMENTS BELOW
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      IMPLICIT INTEGER ( I-N )
      
      DIMENSION STAR(6), OBSERV(6), SKYPOS(7),
     .     POS(3), VEL(3), POSE(3), VTER(3), VCEL(3),
     .     SSS(3), VALUES(500) 
     
      CHARACTER NAMES(500)*6, TXTEPH*158, FILNAM*80

      INTEGER DENUM

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           ) 

      DATA STAR, OBSERV / 12 * 0.D0 /

      DATA IYEAR, MONTH, IDAY, HOUR, LEAPS, UT1UTC, XP, YP /
     .     2008, 4, 24, 10.605D0, 33, -0.387845D0, -0.002D0, +0.529D0 /
     
      DATA GLON, GLAT, HT / -70.D0, 42.D0, 0.D0 /
      
      DATA LU, FILNAM / 20, 'SS_EPHEM.TXT' /


*     FORMAT STATEMENTS FOR OUTPUT

 9010 FORMAT ( 1X, 'JPL ephemeris DE', I3, ' open. Start JD = ', F10.2, 
     .         '  End JD = ', F10.2 ) 
 9020 FORMAT ( A158 )
 9030 FORMAT ( 1X, 3 ( F15.10, 8X ) )
 9040 FORMAT ( 1X, 2 ( F15.6, 8X ), 1F16.11 )
 9050 FORMAT ( 1X, 2 ( F15.10, 8X ), 1F15.12 )
 9060 FORMAT ( 1X, 2 ( F16.11, 8X ), 1F15.10 )
 
*     CHECK WHICH JPL EPHEMERIS IS AVAILABLE
*     REMOVE THIS BLOCK FOR USE WITH SOLSYS 1
      CALL CONST ( NAMES, VALUES, SSS, N )
      DENUM = VALUES ( 1 )
      EPHBEG = SSS ( 1 )
      EPHEND = SSS ( 2 )
      WRITE ( *, 9010 )  DENUM, EPHBEG, EPHEND  
      
*     VERIFY THAT TEXT EPHEMERIS FILE IS AVAILABLE
*     UNCOMMENT THIS BLOCK FOR USE WITH SOLSYS 1
*     OPEN ( UNIT=LU, FILE=FILNAM, STATUS='UNKNOWN' )
*     READ ( LU, 9020 ) TXTEPH
*     WRITE ( *, 9020 ) TXTEPH
*     CLOSE ( LU )

      WRITE ( *, * )
      WRITE ( *, * ) 'NOVAS Sample Calculations'
      WRITE ( *, * ) '-------------------------'
  
*     SETUP CALLS
*     HIGH ACCURACY AND EQUINOX MODE ARE NOVAS DEFAULTS.
      CALL HIACC
      CALL EQINOX
      IEARTH = IDSS ( 'EARTH' )
    
*     WRITE OUT ASSUMED LONGITUDE, LATITUDE, HEIGHT (ITRS = WGS-84)
      WRITE ( *, * )
      WRITE ( *, * ) 'Geodetic location:'
      WRITE ( *, 9030 ) GLON, GLAT, HT
    
*     ESTABLISH TIME ARGUMENTS
      CALL JULDAT  ( IYEAR, MONTH, IDAY, HOUR,   UTCJD )
      TTJD = UTCJD + ( LEAPS + 32.184D0 ) / 86400.D0
      UT1JD = UTCJD + UT1UTC / 86400.D0
      DELTAT = 32.184D0 + LEAPS - UT1UTC
      CALL SETDT ( DELTAT )
      WRITE ( *, * )
      WRITE ( *, * ) 'TT and UT1 Julian Dates and Delta-T:'
      WRITE ( *, 9040 ) TTJD, UT1JD, DELTAT
      
*     APPARENT AND TOPOCENTRIC PLACE OF STAR FK6 1307 = GRB 1830
      RA2000 = 11.88299133D0
      DC2000 = 37.71867646D0
      PMRA   =  4003.27D0
      PMDEC  = -5815.07D0
      PARX   =   109.21D0
      RV     =   -98.8D0
      CALL APSTAR ( TTJD, IEARTH, RA2000, DC2000, PMRA, PMDEC, PARX, RV,
     .     RA, DEC )
      CALL TPSTAR ( UT1JD, GLON, GLAT, HT,   RAT, DECT )
      WRITE ( *, * )
      WRITE ( *, * ) 'FK6 1307 geocentric and topocentric positions:'
      WRITE ( *, 9030 ) RA, DEC
      WRITE ( *, 9030 ) RAT, DECT

*     APPARENT AND TOPOCENTRIC PLACE OF THE MOON
      MOON = IDSS ( 'MOON' )
      CALL APPLAN ( TTJD, MOON, IEARTH,   RA, DEC, DIS )
      CALL TPPLAN ( UT1JD, GLON, GLAT, HT,   RAT, DECT, DIST )
      WRITE ( *, * )
      WRITE ( *, * ) 'Moon geocentric and topocentric positions:'
      WRITE ( *, 9050 ) RA, DEC, DIS
      WRITE ( *, 9050 ) RAT, DECT, DIST
      
*     TOPOCENTRIC PLACE OF THE MOON -- SHOULD BE SAME AS ABOVE
      OBSERV(1) = GLON
      OBSERV(2) = GLAT
      OBSERV(3) = HT
C      OBSERV(4) = ( TTJD - UT1JD ) * 86400.D0 
      CALL PLACE ( TTJD, 'MOON', 1, 1, STAR, OBSERV, SKYPOS )
      RAT  = SKYPOS(4)
      DECT = SKYPOS(5)
      DIST = SKYPOS(6)
      WRITE ( *, 9050 ) RAT, DECT, DIST

*     POSITION OF THE MOON IN LOCAL HORIZON COORDINATES
*     (POLAR MOTION IGNORED HERE)
      CALL ZDAZ ( UT1JD, 0.D0, 0.D0, GLON, GLAT, HT, RAT, DECT, 1,
     .              ZD, AZ, RAR, DECR )
      WRITE ( *, * )
      WRITE ( *, * ) 'Moon zenith distance and azimuth:'
      WRITE ( *, 9030 ) ZD, AZ

*     GREENWICH APPARENT SIDEREAL TIME AND EARTH ROTATION ANGLE
      CALL SIDTIM ( UT1JD, 0.D0, 1, GAST )
      ST = GAST + GLON / 15.D0
      IF ( ST .GE. 24.D0 ) ST = ST - 24.D0
      IF ( ST .LT.  0.D0 ) ST = ST + 24.D0
      CALL EROT ( UT1JD, 0.D0,   ERA )
      WRITE ( *, * )
      WRITE ( *, * ) 'Greenwich and local sidereal time and ',
     .               'Earth Rotation Angle:'
      WRITE ( *, 9060 ) GAST, ST, ERA
      
*     HELIOCENTRIC POSITION OF MARS IN BCRS      
      TDBJD = TTJD
*     ABOVE IS AN APPROXIMATION GOOD TO 0.0017 SECONDS 
*     APPROXIMATION COULD LEAD TO ERROR OF ~50 M IN POSITION OF MARS
      MARS = IDSS ( 'MARS' )
      CALL SOLSYS ( TDBJD, MARS, 1,   POS, VEL, IERR )
      IF ( IERR .EQ. 0 ) THEN
          CALL EQEC ( 0.D0, 0, POS,   POSE )
          CALL ANGLES ( POSE,   ELON, ELAT )
          WRITE ( *, * )
          WRITE ( *, * ) 'Mars heliocentric ecliptic longitude and ',
     .                   'latitude and radius vector:'
          R = DSQRT (POSE(1)**2 + POSE(2)**2 +POSE(3)**2)
          WRITE ( *, 9050 ) ELON * 15.D0, ELAT, R 
      ELSE
          WRITE ( *, * ) 'For Mars, SOLSYS returns error ', IERR
      END IF
      
*     TERRESTRIAL TO CELESTIAL TRANSFORMATION
      SINLON = DSIN ( GLON * DEGRAD ) 
      COSLON = DCOS ( GLON * DEGRAD ) 
      SINLAT = DSIN ( GLAT * DEGRAD ) 
      COSLAT = DCOS ( GLAT * DEGRAD ) 
*     FORM VECTOR TOWARD LOCAL ZENITH (ORTHOGONAL TO ELLIPSOID) IN ITRS
      VTER(1) = COSLAT * COSLON
      VTER(2) = COSLAT * SINLON
      VTER(3) = SINLAT
*     TRANSFORM VECTOR TO GCRS
      CALL TERCEL ( UT1JD, 0.D0, XP, YP, VTER,   VCEL )
      CALL ANGLES ( VCEL,   RA, DEC )
      WRITE ( *, * ) 
      WRITE ( *, * ) 'Direction of zenith vector (RA & Dec) in GCRS:'
      WRITE ( *, 9030 ) RA, DEC
      
      WRITE ( *, * )
      STOP
      END
      
      