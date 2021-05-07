*** SOLSYS VERSION 2 PACKAGE: SOLSYS, AUXPOS, IDSS ***
   
      SUBROUTINE SOLSYS (TJD,BODY,ORIGIN, POS,VEL,IERR)
*
*     SUBROUTINE SOLSYS VERSION 2.
*
*----------------------------------------------------------------------
*
*---PURPOSE: THIS IS SOLSYS VERSION 2.  IT IS INTENDED TO PROVIDE
*            AN INTERFACE BETWEEN THE JPL BINARY DIRECT-ACCESS SOLAR
*            SYSTEM EPHEMERIDES AND THE 'NOVAS' ASTROMETRIC SUBROUTINE
*            LIBRARY.
*
*---REFERENCE:  JPL. 2007, JPL Planetary and Lunar Ephemerides: Export 
*               Information, (Pasadena, CA: JPL) 
*               http://ssd.jpl.nasa.gov/?planet_eph_export.
*
*---INPUT ARGUMENTS:     TJD = JULIAN DATE OF THE DESIRED TIME,
*                              OR FRACTION OF A DAY (SEE NOTE 4),
*                              ON THE T_EPH OR TDB TIME SCALE
*                              (DOUBLE PRECISION).
*                       BODY = BODY IDENTIFICATION NUMBER FOR THE
*                              SOLAR SYSTEM OBJECT OF INTEREST;
*                              MERCURY= 1,...,PLUTO= 9, SUN= 10,
*                              MOON= 11 (INTEGER).
*                     ORIGIN = ORIGIN CODE; SOLAR SYSTEM BARYCENTER= 0,
*                              CENTER OF MASS OF THE SUN= 1 (INTEGER).
*
*---OUTPUT ARGUMENTS:    POS = POSITION VECTOR OF 'BODY' AT TJD;
*                              EQUATORIAL RECTANGULAR COORDINATES,
*                              REFERRED TO ICRS AXES, COMPONENTS IN
*                              AU (DOUBLE PRECISION).
*                        VEL = VELOCITY VECTOR OF 'BODY' AT TJD;
*                              EQUATORIAL RECTANGULAR COORDINATES,
*                              REFERRED TO ICRS AXES, COMPONENTS IN
*                              AU/DAY (DOUBLE PRECISION).
*                       IERR = 0 ... EVERYTHING OK
*                            = 1 ... 'TJD' BEFORE FIRST EPHEMERIS DATE
*                            = 2 ... 'TJD' AFTER LAST EPHEMERIS DATE
*                            = 3 ... INVALID VALUE OF 'BODY' OR
*                                    'ORIGIN' (INTEGER).
*
*---COMMON BLOCKS: NONE.
*
*---SUBROUTINES CALLED: SUBROUTINE CONST   (JPL)
*                       SUBROUTINE PLEPH   (JPL)
*                       SUBROUTINE AUXPOS  (SUPPLIED)
*
*---VERSION/DATE/PROGRAMMER: V1/02-90/JAB
*                            V2/07-91/GHK
*                            V3/05-98/GHK
*                            V4/02-04/GHK
*
*---NOTES: 1. SUBROUTINE PLEPH IS A JPL-SUPPLIED ROUTINE THAT CALLS
*             A VARIETY OF OTHER JPL SUBROUTINES.
*          2. THIS ROUTINE IS FOR USE WITH 1997 VERSION OF JPL
*             EPHEMERIS SOFTWARE, E.G., AS DISTRIBUTED ON CD-ROM
*             'JPL PLANETARY AND LUNAR EPHEMERIDES' (C)1997 BY JPL.
*             IT MAY NOT WORK PROPERLY WITH PREVIOUS VERSIONS OF THE
*             JPL SOFTWARE.
*          3. FOR BODY IDENTIFICATION NUMBERS OUTSIDE THE RANGE 1-11,
*             SUBROUTINE 'AUXPOS' WILL BE CALLED TO SUPPLY POSITIONS
*             AND VELOCITIES FROM SOURCES EXTERNAL TO THE JPL
*             EPHEMERIDES.  A DUMMY VERSION OF THIS ROUTINE IS
*             PROVIDED, WHICH CAN BE REPLACED BY THE USER.
*          4. IF TJD IS BETWEEN -1.D0 AND +1.D0, IT IS ASSUMED TO
*             REPRESENT A FRACTION OF A DAY, WITH THE WHOLE DAYS PART
*             OF THE JULIAN DATE TAKEN FROM A PREVIOUS CALL.
*
*----------------------------------------------------------------------

      INTEGER BODY,ORIGIN,IERR,TARG,CENT,I,N

      DOUBLE PRECISION TJD,POS(3),VEL(3),POSVEL(6),VALUES(500),SSS(3),
     .     BEGJD,ENDJD,TJD1,TJD2(2),TLAST,DINT

      CHARACTER NAMES(500)*6

      LOGICAL FIRST

      SAVE FIRST, BEGJD, ENDJD, TLAST

      DATA FIRST / .TRUE. /

    3 FORMAT ( ' SOLSYS: ERROR ',I1,' AT JD ', F10.1, ', BODY ID ',
     .     I2 )

*---ON FIRST CALL, CALL JPL ROUTINE 'CONST' TO OBTAIN BEGINNING
*   AND ENDING JULIAN DATES OF EPHEMERIS
      IF ( FIRST ) THEN
          CALL CONST ( NAMES, VALUES, SSS, N )
          BEGJD = SSS(1)
          ENDJD = SSS(2)
          TLAST = 0.D0
          FIRST = .FALSE.
      END IF

*---INITIALIZE OUTPUT ARGUMENTS
      POS(1) =  0.D0
      POS(2) =  0.D0
      POS(3) = 99.D0
      VEL(1) =  0.D0
      VEL(2) =  0.D0
      VEL(3) =  0.D0
      IERR = 0
      
*---SET UP SPLIT JULIAN DATE      
      IF ( DABS ( TJD ) .LE. 1.D0 ) THEN
          TJD2(1) = TLAST
          TJD2(2) = TJD
      ELSE
          TLAST = DINT ( TJD )
          TJD2(1) = TLAST
          TJD2(2) = TJD - TLAST
      END IF
      TJD1 = TJD2(1) + TJD2(2) 

*---PERFORM SANITY CHECKS ON THE INPUT BODY AND ORIGIN.
      IF ( ( ORIGIN .LT. 0 ) .OR. ( ORIGIN .GT. 1 ) ) THEN
          IERR = 3
          GO TO 99
      ELSE IF ( ( BODY .LT. 1 ) .OR. ( BODY .GT. 11 ) ) THEN
*         CALL AUXPOS FOR AUXILIARY BODIES (IF ANY)
          CALL AUXPOS ( TJD1, BODY, ORIGIN,   POS, VEL, JERR )
          IERR = JERR
          GO TO 99
      ENDIF

*---CHECK THAT REQUESTED JULIAN DATE IS WITHIN RANGE OF EPHEMERIS.
      IF ( TJD1 .LT. BEGJD ) THEN
          IERR = 1
          GO TO 99
      ELSE IF ( TJD1 .GT. ENDJD ) THEN
          IERR = 2
          GO TO 99
      ENDIF

*---SELECT 'TARG' ACCORDING TO VALUE OF 'BODY'.
      IF ( BODY .EQ. 10 ) THEN
          TARG = 11
      ELSE IF ( BODY .EQ. 11 ) THEN
          TARG = 10
      ELSE
          TARG = BODY
      ENDIF

*---SELECT 'CENT' ACCORDING TO THE VALUE OF 'ORIGIN'.
      IF ( ORIGIN .EQ. 0 ) THEN
          CENT = 12
      ELSE IF ( ORIGIN .EQ. 1 ) THEN
          CENT = 11
      ENDIF

*---CALL JPL ROUTINE 'DPLEPH' TO OBTAIN POSITION AND VELOCITY ARRAY
*   'POSVEL'.
      
      CALL DPLEPH ( TJD2, TARG, CENT,   POSVEL)
      
*     ABOVE IS EQUIVALENT TO CALL PLEPH ( TJD1, TARG, CENT,   POSVEL )
*     IF JULIAN DATE IS NOT SPLIT

*---DECOMPOSE 'POSVEL' INTO POSITION 'POS' AND VELOCITY 'VEL'.
      DO 10 I = 1, 3
          POS(I) = POSVEL(I)
          VEL(I) = POSVEL(I+3)
   10 CONTINUE

*---ROTATE POSITION AND VELOCITY FROM DYNAMICAL TO ICRS FRAME
*   (NECESSARY ONLY FOR JPL EPHEMERIDES PRIOR TO DE405/LE405)
*     CALL FRAME (POSVEL(1),-1,POS)
*     CALL FRAME (POSVEL(4),-1,VEL)

   99 IF ( IERR .NE. 0 ) WRITE ( *, 3 ) IERR, TJD1, BODY

      RETURN

      END



      SUBROUTINE AUXPOS (TJD,M,K, POS,VEL,JERR)
*
*     FOR USE WITH SUBROUTINE SOLSYS VERSION 2.
*     THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF
*     AN AUXILIARY SOLAR SYSTEM BODY AT TIME TJD.  IT IS CALLED
*     FROM SOLSYS VERSION 2 (JPL EPHEMERIS ACCESS) WHEN THE BODY
*     IDENTIFICATION NUMBER IS OUTSIDE THE RANGE 1-11.  ITS
*     INTENDED USE IS TO SUPPLY EPHEMERIS DATA FROM NON-STANDARD
*     SOURCES.
*
*          TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
*          M    = BODY IDENTIFICATION NUMBER (IN)
*          K    = ORIGIN SELECTION CODE (IN)
*                 SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
*                 SET K=1 FOR ORIGIN AT CENTER OF SUN
*          POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO ICRS AXES, COMPONENTS
*                 IN AU (OUT)
*          VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO ICRS AXES, COMPONENTS
*                 IN AU/DAY (OUT)
*          JERR = ERROR INDICATOR (OUT)
*                 JERR=0 MEANS EVERYTHING OK
*                 JERR=1 MEANS TJD BEFORE FIRST VALID DATE
*                 JERR=2 MEANS TJD AFTER LAST VALID DATE
*                 JERR=3 MEANS INVALID VALUE OF M OR K
*
*
*----------------------------------------------------------------------

*     THIS IS A DUMMY VERSION OF SUBROUTINE AUXPOS.  IT SIMPLY
*     RETURNS AN ERROR CODE, SINCE IT CANNOT SUPPLY POSITIONS
*     OR VELOCITIES.  FOR NORMAL (CORRECT) ACCESS TO THE JPL
*     EPHEMERIDES VIA SOLSYS VERSION 2, THIS ROUTINE WILL NEVER
*     BE CALLED.
*
*     A WORKING VERSION MUST BE SUPPLIED BY THE USER ONLY IF THE
*     POSITIONS/VELOCITIES OF AUXILIARY SOLAR SYSTEM BODIES (E.G.,
*     ASTEROIDS) ARE OF INTEREST.  SUCH DATA COULD BE OBTAINED FROM
*     EPHEMERIS FILES OR CLOSED-FORM SERIES APPROXIMATIONS.  THE
*     BODY IDENTIFICATION NUMBERS USED FOR SUCH OBJECTS MUST BE
*     OUTSIDE THE RANGE 1-11 USED FOR MAJOR SOLAR SYSTEM BODIES.
*
*     GENERALLY, IN SUCH CASES IT WOULD BE NECESSARY FOR THIS ROUTINE
*     TO PROVIDE ONLY BARYCENTRIC POSITIONS FOR THE INPUT JD.
*     VELOCITIES ARE NEEDED ONLY FOR THE EARTH, AND HELIOCENTRIC
*     POSITIONS ARE NOT USED IN NOVAS.  ALSO, DO NOT USE FORTRAN
*     I/O UNIT NUMBER 12, WHICH IS USED BY THE JPL ROUTINES.

*----------------------------------------------------------------------

      DOUBLE PRECISION TJD,POS,VEL
      DIMENSION POS(3), VEL(3)

      POS(1) =  0.D0
      POS(2) =  0.D0
      POS(3) = 99.D0
      VEL(1) =  0.D0
      VEL(2) =  0.D0
      VEL(3) =  0.D0
      JERR = 3
      RETURN

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
*     NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=2 IF SOLSYS PROCESSES
*     SPLIT JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE    
*
*     NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
*     THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.
*
*
      CHARACTER NAME*(*), NAMEIN*3, NAMES*3
      DIMENSION NAMES(35), IDS(35)

      DATA NAMES / 'SUN', 'MOO', 'EAR', 'MER', 'VEN', 'MAR', 'JUP',
     .             'SAT', 'URA', 'NEP', 'PLU', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---'  /
      DATA IDS   /    10,    11,     3,     1,     2,     4,     5,
     .                 6,     7,     8,     9,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0  /
      DATA NUM   / 11 /

   3  FORMAT ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', A )

      IDSS = -9999
      NAMEIN = NAME

C     LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
      DO 20 I = 1, NUM
          IF ( NAMEIN .EQ. NAMES(I) ) THEN
              IDSS = IDS(I)
              GO TO 30
          END IF
  20  CONTINUE
  
C     IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES   
      IF ( NAMEIN .EQ. 'JD ' ) THEN
C         IN THIS CASE, SET IDSS=2 IF SOLSYS PROCESSES SPLIT
C         JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE 
          IDSS = 2
          GO TO 30
      END IF    

      WRITE ( *, 3 ) NAME

  30  RETURN

      END

