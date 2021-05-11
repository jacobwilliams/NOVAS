!** SOLSYS VERSION 2 PACKAGE: SOLSYS, AUXPOS, IDSS ***

!***********************************************************************
!>
!  SUBROUTINE SOLSYS VERSION 2.
!
!----------------------------------------------------------------------
!
!---PURPOSE: THIS IS SOLSYS VERSION 2.  IT IS INTENDED TO PROVIDE
!            AN INTERFACE BETWEEN THE JPL BINARY DIRECT-ACCESS SOLAR
!            SYSTEM EPHEMERIDES AND THE 'NOVAS' ASTROMETRIC SUBROUTINE
!            LIBRARY.
!
!---REFERENCE:  JPL. 2007, JPL Planetary and Lunar Ephemerides: Export
!               Information, (Pasadena, CA: JPL)
!               http://ssd.jpl.nasa.gov/?planet_eph_export.
!
!---INPUT ARGUMENTS:     TJD = JULIAN DATE OF THE DESIRED TIME,
!                              OR FRACTION OF A DAY (SEE NOTE 4),
!                              ON THE T_EPH OR TDB TIME SCALE
!                              (DOUBLE PRECISION).
!                       BODY = BODY IDENTIFICATION NUMBER FOR THE
!                              SOLAR SYSTEM OBJECT OF INTEREST;
!                              MERCURY= 1,...,PLUTO= 9, SUN= 10,
!                              MOON= 11 (INTEGER).
!                     ORIGIN = ORIGIN CODE; SOLAR SYSTEM BARYCENTER= 0,
!                              CENTER OF MASS OF THE SUN= 1 (INTEGER).
!
!---OUTPUT ARGUMENTS:    POS = POSITION VECTOR OF 'BODY' AT TJD;
!                              EQUATORIAL RECTANGULAR COORDINATES,
!                              REFERRED TO ICRS AXES, COMPONENTS IN
!                              AU (DOUBLE PRECISION).
!                        VEL = VELOCITY VECTOR OF 'BODY' AT TJD;
!                              EQUATORIAL RECTANGULAR COORDINATES,
!                              REFERRED TO ICRS AXES, COMPONENTS IN
!                              AU/DAY (DOUBLE PRECISION).
!                       IERR = 0 ... EVERYTHING OK
!                            = 1 ... 'TJD' BEFORE FIRST EPHEMERIS DATE
!                            = 2 ... 'TJD' AFTER LAST EPHEMERIS DATE
!                            = 3 ... INVALID VALUE OF 'BODY' OR
!                                    'ORIGIN' (INTEGER).
!
!---COMMON BLOCKS: NONE.
!
!---SUBROUTINES CALLED: SUBROUTINE CONST   (JPL)
!                       SUBROUTINE PLEPH   (JPL)
!                       SUBROUTINE AUXPOS  (SUPPLIED)
!
!---VERSION/DATE/PROGRAMMER: V1/02-90/JAB
!                            V2/07-91/GHK
!                            V3/05-98/GHK
!                            V4/02-04/GHK
!
!---NOTES: 1. SUBROUTINE PLEPH IS A JPL-SUPPLIED ROUTINE THAT CALLS
!             A VARIETY OF OTHER JPL SUBROUTINES.
!          2. THIS ROUTINE IS FOR USE WITH 1997 VERSION OF JPL
!             EPHEMERIS SOFTWARE, E.G., AS DISTRIBUTED ON CD-ROM
!             'JPL PLANETARY AND LUNAR EPHEMERIDES' (C)1997 BY JPL.
!             IT MAY NOT WORK PROPERLY WITH PREVIOUS VERSIONS OF THE
!             JPL SOFTWARE.
!          3. FOR BODY IDENTIFICATION NUMBERS OUTSIDE THE RANGE 1-11,
!             SUBROUTINE 'AUXPOS' WILL BE CALLED TO SUPPLY POSITIONS
!             AND VELOCITIES FROM SOURCES EXTERNAL TO THE JPL
!             EPHEMERIDES.  A DUMMY VERSION OF THIS ROUTINE IS
!             PROVIDED, WHICH CAN BE REPLACED BY THE USER.
!          4. IF TJD IS BETWEEN -1.D0 AND +1.D0, IT IS ASSUMED TO
!             REPRESENT A FRACTION OF A DAY, WITH THE WHOLE DAYS PART
!             OF THE JULIAN DATE TAKEN FROM A PREVIOUS CALL.

subroutine solsys (tjd,body,origin, pos,vel,ierr)

integer body,origin,ierr,targ,cent,i,n

double precision tjd,pos(3),vel(3),posvel(6),values(500),sss(3), &
     begjd,endjd,tjd1,tjd2(2),tlast,dint

character names(500)*6

logical first

save first, begjd, endjd, tlast

data first / .true. /

3 format ( ' SOLSYS: ERROR ',i1,' AT JD ', f10.1, ', BODY ID ', &
     i2 )

!---ON FIRST CALL, CALL JPL ROUTINE 'CONST' TO OBTAIN BEGINNING
!   AND ENDING JULIAN DATES OF EPHEMERIS
if ( first ) then
    call const ( names, values, sss, n )
    begjd = sss(1)
    endjd = sss(2)
    tlast = 0.d0
    first = .false.
end if

!---INITIALIZE OUTPUT ARGUMENTS
pos(1) =  0.d0
pos(2) =  0.d0
pos(3) = 99.d0
vel(1) =  0.d0
vel(2) =  0.d0
vel(3) =  0.d0
ierr = 0

!---SET UP SPLIT JULIAN DATE
if ( dabs ( tjd ) <= 1.d0 ) then
    tjd2(1) = tlast
    tjd2(2) = tjd
else
    tlast = dint ( tjd )
    tjd2(1) = tlast
    tjd2(2) = tjd - tlast
end if
tjd1 = tjd2(1) + tjd2(2)

!---PERFORM SANITY CHECKS ON THE INPUT BODY AND ORIGIN.
if ( ( origin < 0 ) .or. ( origin > 1 ) ) then
    ierr = 3
    go to 99
else if ( ( body < 1 ) .or. ( body > 11 ) ) then
!         CALL AUXPOS FOR AUXILIARY BODIES (IF ANY)
    call auxpos ( tjd1, body, origin,   pos, vel, jerr )
    ierr = jerr
    go to 99
endif

!---CHECK THAT REQUESTED JULIAN DATE IS WITHIN RANGE OF EPHEMERIS.
if ( tjd1 < begjd ) then
    ierr = 1
    go to 99
else if ( tjd1 > endjd ) then
    ierr = 2
    go to 99
endif

!---SELECT 'TARG' ACCORDING TO VALUE OF 'BODY'.
if ( body == 10 ) then
    targ = 11
else if ( body == 11 ) then
    targ = 10
else
    targ = body
endif

!---SELECT 'CENT' ACCORDING TO THE VALUE OF 'ORIGIN'.
if ( origin == 0 ) then
    cent = 12
else if ( origin == 1 ) then
    cent = 11
endif

!---CALL JPL ROUTINE 'DPLEPH' TO OBTAIN POSITION AND VELOCITY ARRAY
!   'POSVEL'.

call dpleph ( tjd2, targ, cent,   posvel)

!     ABOVE IS EQUIVALENT TO CALL PLEPH ( TJD1, TARG, CENT,   POSVEL )
!     IF JULIAN DATE IS NOT SPLIT

!---DECOMPOSE 'POSVEL' INTO POSITION 'POS' AND VELOCITY 'VEL'.
do i = 1, 3
    pos(i) = posvel(i)
    vel(i) = posvel(i+3)
end do

!---ROTATE POSITION AND VELOCITY FROM DYNAMICAL TO ICRS FRAME
!   (NECESSARY ONLY FOR JPL EPHEMERIDES PRIOR TO DE405/LE405)
!     CALL FRAME (POSVEL(1),-1,POS)
!     CALL FRAME (POSVEL(4),-1,VEL)

99 if ( ierr /= 0 ) write ( *, 3 ) ierr, tjd1, body

end subroutine solsys
!***********************************************************************

!***********************************************************************
!>
!  FOR USE WITH SUBROUTINE SOLSYS VERSION 2.
!  THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF
!  AN AUXILIARY SOLAR SYSTEM BODY AT TIME TJD.  IT IS CALLED
!  FROM SOLSYS VERSION 2 (JPL EPHEMERIS ACCESS) WHEN THE BODY
!  IDENTIFICATION NUMBER IS OUTSIDE THE RANGE 1-11.  ITS
!  INTENDED USE IS TO SUPPLY EPHEMERIS DATA FROM NON-STANDARD
!  SOURCES.
!
!       TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
!       M    = BODY IDENTIFICATION NUMBER (IN)
!       K    = ORIGIN SELECTION CODE (IN)
!              SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
!              SET K=1 FOR ORIGIN AT CENTER OF SUN
!       POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
!              COORDINATES, REFERRED TO ICRS AXES, COMPONENTS
!              IN AU (OUT)
!       VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
!              COORDINATES, REFERRED TO ICRS AXES, COMPONENTS
!              IN AU/DAY (OUT)
!       JERR = ERROR INDICATOR (OUT)
!              JERR=0 MEANS EVERYTHING OK
!              JERR=1 MEANS TJD BEFORE FIRST VALID DATE
!              JERR=2 MEANS TJD AFTER LAST VALID DATE
!              JERR=3 MEANS INVALID VALUE OF M OR K
!
!----------------------------------------------------------------------
!
!  THIS IS A DUMMY VERSION OF SUBROUTINE AUXPOS.  IT SIMPLY
!  RETURNS AN ERROR CODE, SINCE IT CANNOT SUPPLY POSITIONS
!  OR VELOCITIES.  FOR NORMAL (CORRECT) ACCESS TO THE JPL
!  EPHEMERIDES VIA SOLSYS VERSION 2, THIS ROUTINE WILL NEVER
!  BE CALLED.
!
!  A WORKING VERSION MUST BE SUPPLIED BY THE USER ONLY IF THE
!  POSITIONS/VELOCITIES OF AUXILIARY SOLAR SYSTEM BODIES (E.G.,
!  ASTEROIDS) ARE OF INTEREST.  SUCH DATA COULD BE OBTAINED FROM
!  EPHEMERIS FILES OR CLOSED-FORM SERIES APPROXIMATIONS.  THE
!  BODY IDENTIFICATION NUMBERS USED FOR SUCH OBJECTS MUST BE
!  OUTSIDE THE RANGE 1-11 USED FOR MAJOR SOLAR SYSTEM BODIES.
!
!  GENERALLY, IN SUCH CASES IT WOULD BE NECESSARY FOR THIS ROUTINE
!  TO PROVIDE ONLY BARYCENTRIC POSITIONS FOR THE INPUT JD.
!  VELOCITIES ARE NEEDED ONLY FOR THE EARTH, AND HELIOCENTRIC
!  POSITIONS ARE NOT USED IN NOVAS.  ALSO, DO NOT USE FORTRAN
!  I/O UNIT NUMBER 12, WHICH IS USED BY THE JPL ROUTINES.

subroutine auxpos (tjd,m,k, pos,vel,jerr)

double precision tjd,pos,vel
dimension pos(3), vel(3)

pos(1) =  0.d0
pos(2) =  0.d0
pos(3) = 99.d0
vel(1) =  0.d0
vel(2) =  0.d0
vel(3) =  0.d0
jerr = 3

end subroutine auxpos
!***********************************************************************

!***********************************************************************
!>
!  THIS FUNCTION RETURNS THE ID NUMBER OF A SOLAR SYSTEM BODY
!  FOR THE VERSION OF SOLSYS (OR SOLSYS-AUXPOS COMBINATION) IN USE.
!
!      NAME   = NAME OF BODY WHOSE ID NUMBER IS DESIRED, E.G.,
!               'SUN', 'MOON, 'MERCURY', ETC., EXPRESSED AS ALL
!               UPPER-CASE LETTERS (IN)
!      IDSS   = ID NUMBER OF BODY, FOR USE IN CALLS TO SOLSYS
!               (FUNCTION VALUE RETURNED)
!
!  NOTE 1: IN THIS VERSION, ONLY THE FIRST THREE LETTERS OF THE
!  BODY'S NAME ARE USED FOR IDENTIFICATION.  ALTERNATIVE VERSIONS
!  MIGHT USE MORE LETTERS.
!
!  NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=2 IF SOLSYS PROCESSES
!  SPLIT JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE
!
!  NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
!  THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.

integer function idss ( name )

character name*(*), namein*3, names*3
dimension names(35), ids(35)

data names / 'SUN', 'MOO', 'EAR', 'MER', 'VEN', 'MAR', 'JUP', &
             'SAT', 'URA', 'NEP', 'PLU', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---'  /
data ids   /    10,    11,     3,     1,     2,     4,     5, &
                 6,     7,     8,     9,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0  /
data num   / 11 /

3 format ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', a )

idss = -9999
namein = name

! LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
do i = 1, num
    if ( namein == names(i) ) then
        idss = ids(i)
        return
    end if
end do

    ! IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES
if ( namein == 'JD ' ) then
    ! IN THIS CASE, SET IDSS=2 IF SOLSYS PROCESSES SPLIT
    ! JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE
    idss = 2
    return
end if

write ( *, 3 ) name

end function idss
!***********************************************************************