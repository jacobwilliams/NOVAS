!** SOLSYS VERSION 1 PACKAGE: SOLSYS, FILDEF, IDSS, BLOCK DATA ***

subroutine solsys (tjd,m,k,pos,vel,ierr)
!
!     SUBROUTINE SOLSYS VERSION 1.
!     THIS SUBROUTINE READS A COORDINATE FILE CONTAINING BARYCENTRIC
!     POSITIONS OF SOLAR SYSTEM BODIES AT DAILY INTERVALS AND PROVIDES
!     THE POSITION AND VELOCITY OF BODY M AT EPOCH TJD.
!
!     TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
!     M    = BODY IDENTIFICATION NUMBER (IN)
!     K    = ORIGIN SELECTION CODE (IN)
!            SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
!            SET K=1 FOR ORIGIN AT CENTER OF SUN
!     POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
!            COORDINATES, REFERRED TO ICRS, COMPONENTS IN AU (OUT)
!     VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
!            COORDINATES, REFERRED TO ICRS, COMPONENTS IN AU/DAY (OUT)
!     IERR = ERROR INDICATOR (OUT)
!            IERR=0 MEANS EVERYTHING OK
!            IERR=1 MEANS TJD BEFORE FIRST USABLE DATE IN FILE
!            IERR=2 MEANS TJD AFTER LAST USABLE DATE IN FILE
!            IERR=3 MEANS BAD VALUE OF M
!            IERR=4 MEANS PROBLEM OPENING FILE
!
!     NOTE 1:  INFORMATION ON THE COORDINATE FILE READ IN:
!        - PATH/NAME OF FILE SPECIFIED IN COMMON /SSFILE/
!        - ALL RECORDS ASCII (FORMATTED)
!        - FIRST RECORD:  HEADER OR IDENTIFYING INFORMATION, IGNORED 
!             HERE
!        - OTHER RECORDS:  TDB JULIAN DATE, X,Y,Z COORDINATES OF SUN, 
!             X,Y,Z COORDINATES OF MERCURY, X,Y,Z COORDINATES OF VENUS,
!             X,Y,Z COORDINATES OF EARTH, X,Y,Z COORDINATES OF MARS, ...
!             READ IN AS PER FORMAT IN COMMON /SSFILE/
!        - RECORDS AT FIXED INTERVALS OF +1 DAY OF TDB
!        - X,Y,Z VALUES IN AU WITH RESPECT TO BCRS (ICRS AXES)
!        - FILE MUST CONTAIN AT LEAST THE COORDINATES OF THE SUN AND
!             THE EARTH 
!        - EARTH REFERS TO GEOCENTER, NOT EARTH-MOON BARYCENTER
!     MANY ASPECTS OF THE FILE CAN BE CONTROLLED AT EXECUTION TIME VIA
!     SUBROUTINE FILDEF.  DEFAULTS: FILE PATH/NAME 'SS_EPHEM.TXT'
!     READ ON LOGICAL UNIT 20, CONTAINING, IN EACH RECORD, A TDB 
!     JULIAN DATE AND COORDINATES OF 11 BODIES (SUN, MERCURY, VENUS, 
!     EARTH, ..., PLUTO, MOON), READ IN WITH FORMAT (F10.2,11(3F16.12)).
!
!     NOTE 2:  IN SUCCESSIVE CALLS TO THIS SUBROUTINE, INPUT
!     JULIAN DATES (TJD) SHOULD GENERALLY BE IN ASCENDING ORDER
!     TO AVOID MULTIPLE SEARCHES STARTING AT BEGINNING OF FILE.
!
!     NOTE 3:  THIS SUBROUTINE USES A 7-POINT LAGRANGIAN INTERPOLATION
!     SCHEME ON FIXED INTERVAL DATA, WHICH PRODUCES INTERPOLATION ERRORS
!     THAT VARY WITH BODY. 
!
!
double precision tjd,pos,vel,xjd,xyz,bpos,bvel,t, &
     tbeg,tend,tlast,astart,origin,ti,ak,ai,p, &
     dabs,dfloat
character filnam*80, formt*80
dimension pos(3), vel(3), xjd(13), xyz(13,50,3), &
     bpos(13,3), bvel(13,3)

!     COMMON BLOCK SSFILE CONTAINS INFORMATION ON THE COORDINATE FILE
common /ssfile/ lu,n,filnam,formt

save

data tbeg, tend, tlast, mlast, klast / 0.d0, 1.d10, 0.d0, 0, 0 /

1 format ( a )
3 format ( ' SOLSYS: ERROR ',i1,' AT JD ', f10.1, ', BODY ID ', &
     i2 )
4 format ( ' SOLSYS: ERROR 4 TRYING TO OPEN FILE ', a, ' ON UNIT ', &
     i2 ) 

if ( tlast <= 0.d0 ) then
    open ( unit=lu, file=filnam, status='UNKNOWN', err=84 )
    read ( lu, 1 )
    read ( lu, formt ) t
    rewind lu
    read ( lu, 1 )
    msun = idss ( 'SUN' )
    tlast = 1.d0
    tbeg = t + 8.d0
    npts = 13 
    intpts = 7
    lmiddl = npts / 2 + 1
    lstart = lmiddl - intpts/2 - 1
    astart = lstart
end if

!     CHECK FOR COMMON ERROR CONDITIONS 
ierr = 0
if ( m < 0 .or. m > n - 1 ) go to 73
if ( tjd < tbeg ) go to 71
if ( tjd > tend ) go to 78

!     LOGIC TO DETERMINE BEST WAY TO SEARCH COORDINATE FILE

!     CHECK IF NEEDED DATA ALREADY IN ARRAYS      
if ( dabs ( tjd - tlast ) <= 0.8d0 ) then
    if ( m /= mlast .or. k /= klast ) go to 30 
    go to 60
end if

!     IF INPUT JD LESS THAN LAST JD, START FROM BEGINNING OF FILE
if ( tjd < tlast ) then
    rewind lu
    read ( lu, 1 ) 
    tlast = 1.d0
end if

!     DECIDE ON COURSE OR FINE SEARCH       
if ( tjd - tlast <= 15.d0 ) go to 20

!     COURSE SEARCH THROUGH COORDINATE FILE
15 read ( lu, formt, end=77 ) t
if ( tjd - t > 10.d0 ) go to 15
do 18 i = 1, npts
   read ( lu, formt, end=77 ) t, ((xyz(i,l,j), j=1,3), l=1,n)
   xjd(i) = t
18 continue

!     FINE SEARCH THROUGH COORDINATE FILE
20 do 22 i = 1, npts - 1
    iold = i + 1
    xjd(i) = xjd(iold)
    do 21 l = 1, n
    do 21 j = 1, 3
        xyz(i,l,j) = xyz(iold,l,j)
21     continue          
22 continue
read ( lu, formt, end=77 ) t, ((xyz(npts,l,j), j=1,3), l=1,n)
xjd(npts) = t
if ( dabs ( tjd - xjd(lmiddl) ) > 0.5d0 ) go to 20
tlast = xjd(lmiddl)

!     AT THIS POINT, THE FILE IS POSITIONED CORRECTLY AND ARRAYS
!     XJD AND XYZ ARE FILLED IN 
30 continue

!     FILL ARRAY BPOS WITH DAILY POSITIONS OF BODY M (WITH INDEX = M+1) 
!     IF K=1, MOVE ORIGIN TO SUN
do 40 i = 1, npts
do 40 j = 1, 3
    origin = 0.d0
    if ( k >= 1 ) origin = xyz(i,msun+1,j)
    bpos(i,j) = xyz(i,m+1,j) - origin
40 continue

!     FILL ARRAY BVEL WITH DAILY VELOCITIES OF BODY M
!     COMPUTED FROM NUMERICAL DIFFERENTIATION OF POSITIONS IN ARRAY BPOS
do 50 i = 1, npts
do 50 j = 1, 3
    bvel(i,j) = 0.d0
    if ( i >= 4 .and. i <= 10 ) &
        bvel(i,j) = (          bpos(i+3,j) -  9.d0 * bpos(i+2,j) &
                     + 45.d0 * bpos(i+1,j) - 45.d0 * bpos(i-1,j) &
                     +  9.d0 * bpos(i-2,j) -         bpos(i-3,j) ) &    !
                    / 60.d0
50 continue

mlast = m
klast = k

!     PERFORM LAGRANGIAN INTERPOLATION FOR POSITION AND VELOCITY AT
!     EPOCH TJD
60 ti = tjd - xjd(lmiddl) + lmiddl
do 63 j = 1, 3
    pos(j) = 0.d0
    vel(j) = 0.d0
    do 62 l = 1, intpts
        ak = astart + dfloat(l)
        p = 1.d0
        do 61 i = 1, intpts
            if ( i == l ) go to 61
            ai = astart + dfloat(i)
            p = p * (ti-ai) / (ak-ai)
61         continue
        pos(j) = pos(j) + p * bpos(lstart+l,j)
        vel(j) = vel(j) + p * bvel(lstart+l,j)
62     continue
63 continue
go to 99

71 ierr = 1
go to 80
73 ierr = 1
go to 80
77 tend = t - 6.d0
rewind lu
read ( lu, 1 ) 
tlast = 1.d0
78 ierr = 2
80 write ( *, 3 ) ierr, tjd, m
go to 99
84 ierr = 4
write ( *, 4 ) filnam, lu

99 return

end



subroutine fildef (lun,nbod,filnm,fmt)
!
!     FOR USE WITH SUBROUTINE SOLSYS VERSION 1.
!     THIS SUBROUTINE MAY BE CALLED TO CHANGE THE VALUES IN
!     COMMON BLOCK SSFILE, WHICH CONTAINS INFORMATION ON THE
!     COORDINATE FILE USED BY SUBROUTINE SOLSYS.
!
!     LUN    = FORTRAN LOGICAL UNIT NUMBER TO BE USED FOR
!              COORDINATE FILE (IN)
!     NBOD   = NUMBER OF BODIES WITH COORDINATES IN FILE (IN)
!     FILNM  = CHARACTER VARIABLE CONTAINING PATH AND FILE NAME
!              OF COORDINATE FILE (IN)
!     FMT    = CHARACTER VARIABLE CONTAINING FORMAT STATEMENT,
!              INCLUDING PARENTHESES AND EVERYTHING IN BETWEEN (IN)
!
!     NOTE:  IF LUN OR NBOD IS ZERO OR LESS, OR FILNM OR FMT IS BLANK,
!     THE CORRESPONDING DEFAULT VALUE IS NOT CHANGED.  DEFAULT VALUES
!     ARE SET IN BLOCK DATA FOR COMMON /SSFILE/.
!  
!  
character filnm*(*), fmt*(*), filnam*80, formt*80
common /ssfile/ lu,n,filnam,formt
save

if ( lun >= 1 ) lu = lun

if ( nbod >= 1 ) n = nbod

if ( filnm /= ' ' ) filnam = filnm

if ( fmt /= ' ' ) formt = fmt

return

end



integer function idss ( name )
!
!     FOR USE WITH SOLSYS VERSION 1.
!     THIS FUNCTION RETURNS THE ID NUMBER OF A SOLAR SYSTEM BODY
!     FOR THE VERSION OF SOLSYS (OR SOLSYS-AUXPOS COMBINATION) IN USE.
!     FOR SOLSYS VERSION 1, THE ID NUMBER OF A BODY REFERS TO ITS
!     ORDER WITHIN EACH RECORD OF THE COORDINATE FILE, WITH ID NUMBERS
!     BEGINNING AT 0 FOR THE FIRST BODY (NORMALLY THE SUN).
!
!         NAME   = NAME OF BODY WHOSE ID NUMBER IS DESIRED, E.G.,
!                  'SUN', 'MOON, 'MERCURY', ETC., EXPRESSED AS ALL
!                  UPPER-CASE LETTERS (IN)
!         IDSS   = ID NUMBER OF BODY, FOR USE IN CALLS TO SOLSYS
!                  (FUNCTION VALUE RETURNED)
!
!     NOTE 1: IN THIS VERSION, ONLY THE FIRST THREE LETTERS OF THE
!     BODY'S NAME ARE USED FOR IDENTIFICATION.  ALTERNATIVE VERSIONS
!     MIGHT USE MORE LETTERS.
!
!     NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=1, SINCE SOLSYS 
!     VERSION 1 DOES NOT PROCESS SPLIT JULIAN DATES.    
!
!     NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
!     THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.
!
!
character name*(*), namein*3, names*3
dimension names(50), ids(50)

data names / 'SUN', 'MER', 'VEN', 'EAR', 'MAR', 'JUP', 'SAT', &
             'URA', 'NEP', 'PLU', 'MOO', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---', '---', '---', '---', '---', '---', '---', &
             '---'/
data ids   /     0,     1,     2,     3,     4,     5,     6, &
                 7,     8,     9,    10,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0,     0,     0,     0,     0,     0,     0, &
                 0/
data num   / 11 /

3 format ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', a )

idss = -9999
namein = name

!     LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
do 20 i = 1, num
    if ( namein == names(i) ) then
        idss = ids(i)
        go to 30
    end if
20 continue

!     IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES   
if ( namein == 'JD ' ) then
!         IN THIS CASE, SET IDSS=1, SINCE SOLSYS VERSION 1 DOES NOT
!         PROCESSES SPLIT JULIAN DATES (IN SUCCESSIVE CALLS) 
    idss = 1
    go to 30
end if    

write ( *, 3 ) name

30 return

end



block data
!
!     FOR USE WITH SUBROUTINE SOLSYS VERSION 1.
!     COMMON BLOCK /SSFILE/ CONTAINS INFORMATION ON THE COORDINATE FILE
!     USED BY SUBROUTINE SOLSYS.  THIS BLOCK DATA SEGMENT SETS UP THE
!     DEFAULT VALUES FOR THE PARAMETERS IN /SSFILE/.  THESE DEFAULTS CAN
!     BE ALTERED BY EXECUTABLE STATEMENTS IN THE MAIN PROGRAM OR ANY
!     SUBROUTINE, OR BY A CALL TO SUBROUTINE FILDEF.
!
!     LU     = FORTRAN LOGICAL UNIT NUMBER OF COORDINATE FILE
!     N      = NUMBER OF BODIES WITH COORDINATES IN FILE
!     FILNAM = CHARACTER VARIABLE CONTAINING PATH AND FILE NAME
!              OF COORDINATE FILE
!     FORMT  = CHARACTER VARIABLE CONTAINING FORMAT STATEMENT,
!              INCLUDING PARENTHESES AND EVERYTHING BETWEEN 
!
!
character filnam*80, formt*80

common /ssfile/ lu,n,filnam,formt

data lu / 20 /

data n / 11 /

data filnam / 'SS_EPHEM.TXT' /
      
data formt / '(F10.2,33F16.12)' /

end

